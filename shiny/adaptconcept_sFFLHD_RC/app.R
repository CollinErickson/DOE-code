library(shiny)
#source("CopyOfsFFLHD.R")
library(mlegp)
if (F) { # Use for local runs. Then copy files and comment this section out before uploading
  #library(UGP)
  #library(contourfilled)
  #library(TestFunctions)
  #source("../../sFFLHD.R")
  #source("../../adaptconcept_helpers.R")
  #source("../../LHS.R")
  #source("../../adaptconcept_sFFLHD_RC.R")
} else { # ONLY USE FOR uploading to Shiny
  source("UGP.R")
  source("contourfilled.R")
  source("sFFLHD.R")
  source("adaptconcept_helpers.R")
  source("LHS.R")
  source("adaptconcept_sFFLHD_RC.R")
  source("TestFunctions.R")
}

ui <- fluidPage(
  titlePanel("Adaptive sFFLHD"),
  sidebarLayout(
    sidebarPanel(
      selectInput("func", label = "Function", 
                  choices = list("Gaussian" = "Gaussian", 
                                 "Sinumoid" = "Sinumoid", 
                                 "Diagonal plateau" = "Diagonal plateau",
                                 "Branin (2D)" = "Branin", 
                                 "Franke (2D)" = "Franke",
                                 "Currin (2D)" = "currin1991", 
                                 "Lim (2D)" = "lim2002", 
                                 "Zhou" = "zhou1998"
                                 ), 
                  selected = "Gaussian"),
      selectInput("obj", label = "Objective", 
                  choices = list("MSE" = "mse", 
                                 "Max Error" = "maxerr"
                  ), 
                  selected = "mse"),
      sliderInput("Dimension",
                  "Dimension",
                  min = 1,
                  max = 9,
                  value = 2
      ),
      sliderInput("Batch",
                  "Batch size",
                  min = 1,
                  max = 9,
                  value = 3
      ),
      sliderInput("g",
                  "Divisions",
                  min = 2,
                  max = 6,
                  value = 3
      ),
      tags$head(tags$script(src = "message-handler.js")),
      
      actionButton("addbatch", "Add batch"),
      fluidRow(
        column(7,actionButton("addbatches", "Add batches")),
        column(5,numericInput("numbatches", 
                     label = NULL, 
                     value = 2
                     #,width = "50px"
                     )    )), 
      actionButton("restart", "Restart")
      ,width = 2
    ),
    mainPanel(
      plotOutput("plot", width="975px", height="650px")
    )
  )
)


server <- function(input, output, session) {
  gaussian1 <- function(xx) exp(-sum((xx-.5)^2)/2/.1)
  sinumoid <- function(xx){sum(sin(2*pi*xx*3)) + 20/(1+exp(-80*(xx[[1]]-.5)))}
  diagonal_plateau <- function(xx){sum(xx) > length(xx) / 2}
  D <- reactive({input$Dimension})
  L <- reactive({input$Batch})
  g <- reactive({input$g})
  obj <- reactive({input$obj})
  funcr <- reactive({
    funcin <- input$func
    if (funcin == "Gaussian") {
      return(gaussian1)
    }
    if (funcin == "Sinumoid") {
      return(sinumoid)
    }
    if (funcin == "Diagonal plateau") {
      return(diagonal_plateau)
    }
    if (funcin == "Branin") {
      return(branin)
    }
    if (funcin == "Franke") {
      return(franke)
    }
    if (funcin == "currin1991") {
      return(currin1991)
    }
    if (funcin == "lim2002") {
      return(lim2002)
    }
    if (funcin == "zhou1998") {
      return(zhou1998)
    }
  })
  #s <- sFFLHD.seq(D = D(), L = L())
  #sr <- reactive({ sFFLHD.seq(D = input$Dimension, L = input$Batch) })
  #sr <- reactive({ sFFLHD.seq(D = D(), L = L()) })
  ar <- reactive({ adapt.concept.sFFLHD.RC(D=D(),L=L(),g=g(),func=funcr(), obj=obj()) })
  #D <- renderPrint({input$Dimension})
  #L <- renderPrint({input$Batch})
  #s <- sFFLHD.seq(D = D, L = 3)
  
  observeEvent(input$addbatch, {
    a <- ar() # for reactive
    output$plot <- renderPlot({
      a$run1()
    })
  })
  observeEvent(input$addbatches, {
    a <- ar() # for reactive
    output$plot <- renderPlot({
      nbat <- isolate(input$numbatches) # isolate -> won't run until button pressed
      a$run(nbat, plotlastonly = TRUE) 
    })
  })
  observeEvent(input$restart, {
    session$reload()
  })
}

shinyApp(ui, server)


