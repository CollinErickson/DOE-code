library(shiny)
source("CopyOfsFFLHD.R")

ui <- fluidPage(
  titlePanel("sFFLHD"),
  sidebarLayout(
    sidebarPanel(
    #wellPanel(
      sliderInput("Dimension",
                  "Dimension",
                  min = 1,
                  max = 9,
                  value = 2
                  #width = "400px"
                  ),
      sliderInput("Batch",
                  "Batch size",
                  min = 1,
                  max = 9,
                  value = 3
                  #width = "400px"
                  ),
      tags$head(tags$script(src = "message-handler.js")),
      actionButton("addbatch", "Add batch"),
      actionButton("restart", "Restart"),
      width = 2
    ),
    mainPanel(
      plotOutput("plot", width="1000px", height="1000px")
    )
  )
)


server <- function(input, output, session) {
  D <- reactive({input$Dimension})
  L <- reactive({input$Batch})
  #s <- sFFLHD.seq(D = D(), L = L())
  #sr <- reactive({ sFFLHD.seq(D = input$Dimension, L = input$Batch) })
  sr <- reactive({ sFFLHD.seq(D = D(), L = L()) })
  #D <- renderPrint({input$Dimension})
  #L <- renderPrint({input$Batch})
  #s <- sFFLHD.seq(D = D, L = 3)
  
  observeEvent(input$addbatch, {
    s <- sr() # for reactive
    #if (D != input$Dimension | L != input$Batch) {
    #  s <- sFFLHD.seq(D = D, L = 3)
    #}
    
    s$get.batch()
    output$plot <- renderPlot({
      ppch <- 21+sapply(1:((nrow(s$Xb)/s$L)),function(ii){rep(ii,s$L)})%%5
      pcol <- 1 + (sapply(1:((nrow(s$Xb)/s$L)),function(ii){rep(ii,s$L)})-1)%%8
      if (input$Dimension == 2) {
        ppch <- 21+sapply(1:((nrow(s$Xb)/s$L)),function(ii){rep(ii,s$L)})%%5
        pcol <- 1 + (sapply(1:((nrow(s$Xb)/s$L)),function(ii){rep(ii,s$L)})-1)%%8
        plot(s$Xb,xlim=0:1,ylim=0:1,col=pcol, pch=ppch, bg=pcol, cex=4, cex.axis=2, cex.lab=2,
        xlab="", ylab="")
        abline(h=(0:(s$Lb))/s$Lb,v=(0:(s$Lb))/s$Lb,col=1)
      } else {
        pairs_func <- function(x, y, ...) {
          abline(h=(0:(s$Lb))/s$Lb,v=(0:(s$Lb))/s$Lb)
          #abline(h=(0:(s$lb))/s$lb,v=(0:(s$lb))/s$lb)
          points(x, y, col=pcol, pch=ppch, bg=pcol, cex=4, cex.axis=2, cex.lab=2)
        }
        pairs(s$Xb, lower.panel=pairs_func, upper.panel=pairs_func, xlim=0:1, ylim=0:1)
      }
    })
  })
  observeEvent(input$restart, {
    session$reload()
  })
}

shinyApp(ui, server)


