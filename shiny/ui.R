library(shiny)
library(ggplot2)
library(Cairo)

fluidPage(
    titlePanel("Uploading Files"),
    fluidRow(
        column(width = 5,
               wellPanel(
                       fileInput('file1', 'Choose CSV File',
                                 accept=c('text/csv',
                                          'text/comma-separated-values,text/plain',
                                          '.csv')),
                       tags$hr(),
                       checkboxInput('header', 'Header', TRUE),
                       radioButtons('sep', 'Separator',
                                    c(Comma=',',
                                      Semicolon=';',
                                      Tab='\t'),
                                    ','),
                       radioButtons('quote', 'Quote',
                                    c(None='',
                                      'Double Quote'='"',
                                      'Single Quote'="'"),
                                    '"')
                   ),
                   mainPanel(
                       tableOutput('contents')
                   )
        ),
        column(width = 6,
               wellPanel(
               plotOutput("plot1", height = 350,
                          click = "plot1_click",
                          brush = brushOpts(
                              id = "plot1_brush"
                          )
               ),
               actionButton("exclude_toggle", "Toggle points"),
               actionButton("exclude_reset", "Reset")
               )
        )
    )
)
