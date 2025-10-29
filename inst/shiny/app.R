library(shiny)
library(reticulate) 
# ... (other libraries) ...

# CRITICAL: This tells shinyapps.io to provision a Python environment 
# and install your required packages and versions.
py_require(
  packages = c("scikit-learn==1.5.0", "joblib", "NumPy==1.26.4"),
  python_version = "3.10" # Use a compatible, standard version
)

# inst/shiny/app.R
if (!requireNamespace("seneR", quietly = TRUE)) {
  install.packages("seneR_0.0.0.9000.tar.gz", repos = NULL, type = "source")
}
library(seneR) # Your package
library(Seurat)
library(ggplot2)
library(DT)
library(SeuratData) # For the pbmc3k example

# ----------------- PYTHON CONFIGURATION (CRITICAL FOR DEPLOYMENT) ------------------
# This block tells rsconnect what Python environment to create and packages to install.
# It MUST be in the global scope of app.R (or global.R).



ui <- fluidPage(
  titlePanel("seneR: Senescence Analysis Portal"),
  sidebarLayout(
    sidebarPanel(
      actionButton("init_env", "Initialize seneR environment"),
      selectInput("mode", "Analysis mode:",
                  choices = c("Bulk RNA-seq", "Single-cell (Seurat)")),
      hr(),
      
      conditionalPanel(
        condition = "input.mode == 'Bulk RNA-seq'",
        fileInput("expr_file", "Upload expression matrix (CSV)", accept = ".csv"),
        fileInput("meta_file", "Upload metadata (CSV)", accept = ".csv"),
        checkboxInput("use_example", "Use example dataset (GSE246425)", TRUE),
        actionButton("run_bulk", "Run bulk analysis")
      ),
      
      conditionalPanel(
        condition = "input.mode == 'Single-cell (Seurat)'",
        fileInput("seurat_file", "Upload Seurat object (RDS)", accept = ".rds"),
        checkboxInput("use_pbmc", "Use SeuratData PBMC3K sample", TRUE),
        actionButton("run_seurat", "Run Seurat analysis")
      )
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Summary / Messages", verbatimTextOutput("log")),
        tabPanel("Results Table", DTOutput("result_table")),
        tabPanel("Plot 1", plotOutput("plot1", height = "500px")),
        tabPanel("Plot 2", plotOutput("plot2", height = "500px"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Initialize environment
  observeEvent(input$init_env, {
    output$log <- renderPrint({
      cat("Setting up environment...\n")
      env_Load() # This will now execute the simplified function from A
      cat("Environment ready.")
    })
  })
  
  #### BULK RNA-SEQ MODE ####
  observeEvent(input$run_bulk, {
    output$log <- renderPrint({ cat("Running bulk RNA-seq analysis...\n") })
    
    if (isTRUE(input$use_example)) {
      expr <- read.csv(system.file("demo", "GSE246425_counts.csv", package = "seneR"), row.names = 1)
      meta <- read.csv(system.file("demo", "GSE246425_meta.csv", package = "seneR"), row.names = 1)
    } else {
      req(input$expr_file, input$meta_file)
      expr <- read.csv(input$expr_file$datapath, row.names = 1)
      meta <- read.csv(input$meta_file$datapath, row.names = 1)
    }
    
    # Run senescence prediction
    SID_res <- SenCID(expr)
    score_res <- SID_res$score_res
    
    # Plot group comparison
    output$plot1 <- renderPlot({
      plot_group(score_res, meta, "group", "SID_Score", comparisons = list(c("Old", "Young")))
    })
    
    # GSVA analysis
    gsva_res <- seneGSVA(expr)
    
    # Violin plot
    output$plot2 <- renderPlot({
      plot_violin(gsva_res, meta, "group")
    })
    
    # Results table
    output$result_table <- renderDT({
      datatable(score_res)
    })
    
    output$log <- renderPrint({
      cat("Bulk analysis complete.\nSID and GSVA results generated.")
    })
  })
  
  #### SINGLE-CELL (SEURAT) MODE ####
  observeEvent(input$run_seurat, {
    output$log <- renderPrint({ cat("Running Seurat analysis...\n") })
    
    if (isTRUE(input$use_pbmc)) {
      library(SeuratData)
      options(timeout = 600)
      SeuratData::InstallData("pbmc3k")
      pbmc3k <- LoadData("pbmc3k", type = "pbmc3k.final")
    } else {
      req(input$seurat_file)
      pbmc3k <- readRDS(input$seurat_file$datapath)
    }
    
    expr_matrix <- GetAssayData(pbmc3k, layer = "counts")
    SID_res <- SenCID(expr_matrix, binarize = TRUE)
    pbmc3k@meta.data <- cbind(pbmc3k@meta.data, SID_res$score_res)
    
    # UMAP/Feature plots
    output$plot1 <- renderPlot({
      VlnPlot(pbmc3k, features = "SID_Score", group.by = "seurat_annotations")
    })
    
    output$plot2 <- renderPlot({
      DimPlot(pbmc3k, group.by = "RecSID")
    })
    
    output$result_table <- renderDT({
      datatable(SID_res$score_res)
    })
    
    output$log <- renderPrint({
      cat("Seurat analysis complete.\nSID scores and plots generated.")
    })
  })
}

shinyApp(ui, server)
