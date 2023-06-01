#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)
library(progress)
library(tidyverse)
reportsDir <- "../data"
Sys.setenv("REPORTS_DIR"=reportsDir)
source("definitions.R")
source("functions_io.R")
source("functions_steps.R")
source("functions_methods.R")
source("functions_ggpubr-custom.R")
source("functions_plots.R")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Ceramide Ring Trial"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          # textInput("labId", "Enter an id for your lab", value = "", placeholder = "myLab"),
          # Input: Select a file ----
          fileInput("file1", "Choose XLSX File",
                    multiple = FALSE,
                    accept = c("application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                               ".xlsx")),
          
          # Horizontal line ----
          tags$hr(),
        ),

        # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel(
            tabPanel("Raw Report", tableOutput("rawReport")),
            tabPanel("Summary", verbatimTextOutput("summary")),
            tabPanel("Calibration Lines", tableOutput("calibrationLines")),
            tabPanel("QC Samples", tableOutput("qcSamples")),
            tabPanel("NIST SRM 1950", tableOutput("nistSrm")),
            tabPanel("NIST YAA", tableOutput("nistYaa")),
          )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  data <- reactiveValues()
  
  output$rawReport <- renderTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    # req(input$labId)
    req(input$file1)
    
    # when reading semicolon separated files,
    # having a comma separator causes `read.csv` to error
    tryCatch(
      {
        reportsDir <- file.path(tempdir(),"userReport")
        dir.create(reportsDir, showWarnings = FALSE)
        file.copy(input$file1$datapath, file.path(reportsDir, "mylab.xlsx"), overwrite = TRUE)
        df <- createAssayTable(
          reportsToInclude = c("mylab"), 
          ceramideColNames = ceramideColNames, 
          filePrefix = "", 
          protocol="User", 
          reportsDir=reportsDir,
          blankTypes = blankTypesTable,
          na = naValues
        )
        reactiveValues$data <- df
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    return(head(df))
  })
  
  output$calibrationLines <- renderTable({
    req(reactiveValues$data)
    
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
