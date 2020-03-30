#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)
library(shinyjs)
library(tidyr)
library(dplyr)
source("data-generation-app-fun.R") #source helper function
file_size_max = 10^4
# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Synthetic DECL Data Generator"),
   
   fluidRow(
   column(4, 
          numericInput(inputId = "file_size", "Number of Reads in FASTA file", 10^3, min = 1, max = file_size_max)),
   column(4, 
          p("Status", style ="font-weight: bold"),
          verbatimTextOutput("file_path", placeholder = T)
          )
  ),
   
   fluidRow(
     column(4, style = "background-color:#e6eeff;",
            numericInput(inputId = "len_S", "Sequence length of selection identifiers", 6, min = 1, max = 100),
            numericInput(inputId = "n_S_list", "Number of Selection Identifiers", 2, min = 1, max = 100),
            numericInput(inputId = "len_S_list", "Length of selection identifiers list", 10, min = 1, max = 10)),
   
     column(4, style = "background-color:#ccdcff;",
            numericInput(inputId = "len_B", "Sequence length of building block identifiers", 6, min = 1, max = 100),
            numericInput(inputId = "n_B_list", "Number of building block Identifiers", 2, min = 1, max = 100),
            numericInput(inputId = "len_B_list", "Length of building block identifiers list", 10, min = 1, max = 10)),
    
     column(4, style = "background-color:#b3cbff;",
            numericInput(inputId = "len_C", "Sequence length of constant regions", 6, min = 1, max = 100),
            numericInput(inputId = "len_C_list", "Length of constant regions list", 1, min = 1, max = 10))
   ),
   actionButton("generate_button", label = "Genrate Data")

)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  set_code = c("A", "T", "G", "C")
  set_const = c("a", "t", "g", "c")
  
  x = observeEvent(input$file_size, {
  updateNumericInput(session, "file_size", value = min(file_size_max, input$file_size))}
  )
  
  observeEvent(input$generate_button, {
    output$file_path = renderText({"Generating data set...please wait"})
    
    dirr = strsplit(as.character(Sys.time()), " ")[[1]]
    dirr[2] = gsub(":", ".",dirr[2])
    dirr = paste(dirr[1], dirr[2], sep = "_")
    dir.create(dirr)
    setwd(dirr)
    
    n_C_list = input$n_S_list + input$n_B_list - 1 #number of const regions given by number of identifiers
    set_code[1]
    file = "test.fasta"

    tbl = create_fasta_file(file,
                            input$file_size,
                            set_code,
                            set_const,
                            input$len_S, input$len_S_list, input$n_S_list,
                            input$len_B, input$len_B_list, input$n_B_list,
                            input$len_C, input$len_C_list, n_C_list)

    # create_struct_file("structure.txt", file,
    #                    input$len_S, input$len_B, input$len_C,
    #                    input$n_S_list, input$n_B_list, n_C_list)
    setwd("..")
    output$file_path = renderText({"Generating data set complete"})
   })
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)

