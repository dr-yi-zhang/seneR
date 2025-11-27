###############################################
# app.R (works BOTH locally and on shinyapps.io)
###############################################

library(shiny)
library(reticulate)

# Python settings (must come first)

py_require(
  packages = c(
    "numpy==1.26.4",
    "joblib",
    "scikit-learn==1.5.0"
  ),
  python_version = "3.10"
)


library(DT)
library(ggplot2)
library(Seurat)
library(dplyr)
library(tidyr)
library(seneR) 

options(shiny.maxRequestSize = 150 * 1024^2)   # Safe upload limit
options(expressions = 5e5)                 # Prevent recursion issues


############################################################
# 0. DETECT ENVIRONMENT (LOCAL vs SHINYAPPS.IO)
############################################################

running_local <- !nzchar(Sys.getenv("SHINY_PORT"))
if (running_local) {
  message(">>> Running in LOCAL mode")
} else {
  message(">>> Running on shinyapps.io (DEPLOY mode)")
}


############################################################
# 1. Python Environment Handling
############################################################

if (running_local) {
  # LOCAL ONLY: Install Python dependencies if missing
  message("Setting up Python environment locally...")
  
  py_require(
    packages = c(
      "scikit-learn==1.5.0",
      "joblib",
      "numpy==1.26.4"
    ),
    python_version = "3.10"
  )
  
} else {
  # DEPLOY MODE: shinyapps.io provides Python env from requirements.txt
  message("Using shinyapps.io Python environment.")
}

use_python(Sys.getenv("RETICULATE_PYTHON"), required = FALSE)
py_config()


############################################################
# 2. Load ALL Python Models ONCE (critical to prevent OOM)
############################################################

py_model_env <- new.env(parent = emptyenv())

joblib_load <- function(path) {
  py_run_string(sprintf("
import warnings
from sklearn.exceptions import InconsistentVersionWarning
warnings.filterwarnings('ignore', category=InconsistentVersionWarning)
import joblib
obj = joblib.load(r'%s')
", path))
  py$obj
}

load_python_models <- function() {
  message("Loading Python models...")
  
  base_path <- system.file("model", package = "seneR")
  
  # Load SID models
  for (i in 1:6) {
    m  <- file.path(base_path, sprintf("SID%d.pkl", i))
    ml <- file.path(base_path, sprintf("SID%d_L.pkl", i))
    
    py_model_env[[paste0("model",  i)]]  <- joblib_load(m)
    py_model_env[[paste0("modelL", i)]]  <- joblib_load(ml)
  }
  
  # Load recommendation model
  rec_path <- file.path(base_path, "recommend_model.pkl")
  py_model_env$recommender <- joblib_load(rec_path)
  
  message("Python models loaded successfully.")
}



############################################################
# 3. Helper to access pre-loaded Python models
############################################################

GetModels <- function(sidnum) {
  list(
    model = py_model_env[[paste0("model", sidnum)]],
    model_L = py_model_env[[paste0("modelL", sidnum)]]
  )
}


############################################################
# 4. UI
############################################################

ui <- fluidPage(
  
  titlePanel("seneR: Senescence Analysis Portal"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("mode", "Analysis mode:",
                  c("Bulk RNA-seq", "Single-cell (Seurat)")),
      
      hr(),
      
      conditionalPanel(
        condition = "input.mode == 'Bulk RNA-seq'",
        fileInput("expr_file", "Upload expression matrix (CSV)", accept = ".csv"),
        fileInput("meta_file", "Upload metadata (CSV)", accept = ".csv"),
        checkboxInput("use_example", "Use demo dataset (GSE246425)", TRUE),
        actionButton("run_bulk", "Run Bulk Analysis")
      ),
      
      conditionalPanel(
        condition = "input.mode == 'Single-cell (Seurat)'",
        fileInput("seurat_file", "Upload Seurat RDS", accept = ".rds"),
        checkboxInput("use_pbmc", "Use small PBMC example", TRUE),
        actionButton("run_seurat", "Run Seurat Analysis")
      )
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Log", verbatimTextOutput("log")),
        tabPanel("Table", DTOutput("result_table")),
        tabPanel("Plot 1", plotOutput("plot1", height="500px")),
        tabPanel("Plot 2", plotOutput("plot2", height="500px"))
      )
    )
  )
)

# DEFINE Pred_for_shiny HERE
Pred_for_shiny <- function(cpm_zcol, sidnum, binarize = FALSE) {
  df_test <- t(cpm_zcol)
  rfeatures <- GetFeatures(sidnum)
  x_test <- SplitLabel(as.data.frame(df_test), rfeatures)
  
  py$x_test <- as.matrix(x_test)
  
  py_run_string("
model_decf = model.decision_function(x_test)
Ltest_prob = model_L.predict_proba(x_test)
model_pretest = model_L.predict(x_test)
")
  
  model_decf <- py$model_decf
  Ltest_prob <- py$Ltest_prob
  model_pretest <- py$model_pretest
  
  #### 🔥 MEMORY RELEASE 🔥 ####
  py$x_test <- NULL
  py$model_decf <- NULL
  py$Ltest_prob <- NULL
  py$model_pretest <- NULL
  gc()
  py_run_string("import gc; gc.collect()")
  
  if (binarize) {
    return(data.frame(
      SID_Score = Ltest_prob[, 2],
      Decision = model_decf,
      Binarization = model_pretest
    ))
  } else {
    return(data.frame(
      SID_Score = Ltest_prob[, 2],
      Decision = model_decf
    ))
  }
}

Recommend_stub <- function(cpm_zcol) {
  n <- nrow(t(cpm_zcol))
  rn <- rownames(t(cpm_zcol))
  data.frame(
    RecSID   = factor(rep(NA, n)),
    rec_SID1 = NA, rec_SID2 = NA, rec_SID3 = NA,
    rec_SID4 = NA, rec_SID5 = NA, rec_SID6 = NA,
    row.names = rn
  )
}

############################################################
# 5. SERVER
############################################################

server <- function(input, output, session) {
  
  try({
    load_python_models()
  }, silent = TRUE)
  
  Pred <- Pred_for_shiny
  Recommend  <- Recommend_stub
  
  output$log <- renderPrint(cat("Environment Ready.\n"))
  
  
  ######################
  # BULK MODE
  ######################
  observeEvent(input$run_bulk, {
    output$log <- renderPrint(cat("Running bulk RNA-seq analysis... please wait.\n"))
    
    
    withProgress(message = "Processing...", value = 0, {
      
      incProgress(0.1, detail = "Loading data...")
      # Load data
      if (input$use_example) {
        expr <- read.csv(system.file("demo/GSE246425_counts.csv", package="seneR"),
                         row.names=1, check.names=FALSE)
        meta <- read.csv(system.file("demo/GSE246425_meta.csv", package="seneR"),
                         row.names=1, check.names=FALSE)
        
      } else {
        req(input$expr_file, input$meta_file)
        expr <- read.csv(input$expr_file$datapath, row.names=1, check.names=FALSE)
        meta <- read.csv(input$meta_file$datapath, row.names=1, check.names=FALSE)
      }
      
      # Run prediction
      incProgress(0.2, detail = "Running SID prediction...")
      
      
      SID_res <- SenCID(expr)
      score_res <- SID_res$score_res
      
      output$plot1 <- renderPlot({
        plot_group(score_res, meta, "group", "SID_Score",
                   comparisons=list(c("Old","Young")))
      })
      
      # Safe GSVA subset
      incProgress(0.3, detail = "Running GSVA...")
      
      expr_small <- expr[1:min(200, nrow(expr)), ]
      gsva_res <- seneGSVA(expr_small)
      
      incProgress(0.3, detail = "Rendering plots...")
      output$plot2 <- renderPlot({
        plot_violin(gsva_res, meta, "group")
      })
      
      output$result_table <- renderDT({
        datatable(score_res, options=list(pageLength=20))
      })
      
      output$log <- renderPrint(cat(
        "Bulk RNA-seq analysis completed.\n",
        "SID and GSVA results are ready.\n"
      ))
    })
  })
  
  ######################
  # SEURAT MODE
  ######################
  observeEvent(input$run_seurat, {
    output$log <- renderPrint(cat("Running Seurat analysis...\n"))
    
    
    if (input$use_pbmc) {
      output$log <- renderPrint(cat("use_pbmc...\n"))
      demo_path <- system.file("demo/pbmc_small.rds", package = "seneR")
      pbmc <- readRDS(demo_path)
    } else {
      req(input$seurat_file)
      pbmc <- readRDS(input$seurat_file$datapath)
    }
    
    expr <- GetAssayData(pbmc, layer="counts")
    rm(pbmc)
    gc()
    
    SID_res <- SenCID(expr, binarize=TRUE)
    pbmc@meta.data <- cbind(pbmc@meta.data, SID_res$score_res)
    
    output$plot1 <- renderPlot({
      VlnPlot(pbmc, features="SID_Score", group.by="seurat_annotations")
    })
    
    output$plot2 <- renderPlot({
      DimPlot(pbmc, group.by="RecSID")
    })
    
    output$result_table <- renderDT({
      datatable(SID_res$score_res)
    })
    
    log_msg("Seurat analysis completed.")
  })
}


############################################################
# 6. Launch App
############################################################

shinyApp(ui, server)
