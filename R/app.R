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
dataDir <- "../data"
Sys.setenv("DATA_DIR"=dataDir)
source("definitions.R")
source("functions_io.R")
source("functions_steps.R")
source("functions_methods.R")
source("functions_ggpubr-custom.R")
source("functions_plots.R")

dir.create(outputDir, showWarnings = F, recursive = T)

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
          selectInput("protocol", "Protocol", choices=c("Standard","Preferred"), selected = "Standard", multiple = FALSE)
        ),

        # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel(
            tabPanel("Summary", tableOutput("summary")),
            tabPanel("Comparison", tableOutput("comparison")),
            tabPanel("Calibration Lines", tableOutput("calibrationLines")),
            tabPanel("QC Samples", tableOutput("qcSamples")),
            tabPanel("NIST SRM 1950", tableOutput("nistSrm")),
            tabPanel("NIST YAA", tableOutput("nistYaa")),
            tabPanel("NIST hTAG", tableOutput("nistHtag")),
            tabPanel("NIST T1D", tableOutput("nistT1d")),
            tabPanel("Raw Report", tableOutput("rawReport"))
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
    req(input$protocol)
    
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
          protocol=input$protocol, 
          reportsDir=reportsDir,
          blankTypes = blankTypesTable,
          na = naValues
        ) |> mutate(SampleType = forcats::as_factor(as.character(SampleType)))
        data$table <- df
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    return(df)
  })
  
  output$summary <- renderTable({
    req(data$table)
    data$table |> group_by(LabId, SampleType, ceramideName, isotope) |> summarize(n=n(), min_area=min(area), max_area=max(area))
  })
  
  output$calibrationLines <- renderTable({
    req(data$table)
    data$table |>
      filter(SampleType %in% c("Calibration Line 1", "Calibration Line 2"))
  })
  
  # output$calibrationLinesPlot <- 
  
  output$qcSamples <- renderTable({
    req(data$table)
    data$table |>
      filter(grepl("QC", SampleType))
  })
  
  output$nistSrm <- renderTable({
    req(data$table)
    data$table |>
      filter(grepl("NIST SRM", SampleType))
  })
  
  output$nistYaa <- renderTable({
    req(data$table)
    data$table |>
      filter(grepl("NIST YAA", SampleType))
  })
  
  output$nistHtag <- renderTable({
    req(data$table)
    data$table |>
      filter(grepl("NIST hTAG", SampleType))
  })
  
  output$nistT1d <- renderTable({
    req(data$table)
    data$table |>
      filter(grepl("NIST T1D", SampleType))
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
