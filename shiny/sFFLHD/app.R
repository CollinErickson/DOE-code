library(shiny)

ui <- fluidPage(
  titlePanel("sFFLHD"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("Dimension",
                  "Dimension",
                  min = 1,
                  max = 9,
                  value = 2),
      sliderInput("Batch",
                  "Batch size",
                  min = 1,
                  max = 9,
                  value = 3),
      tags$head(tags$script(src = "message-handler.js")),
      actionButton("do", "Add batch")
    ),
    mainPanel(
      plotOutput("plot", width="1000px", height="1000px")
    )
  )
)


server <- function(input, output, session) {
  #D <- reactive({input$Dimension})
  #s <- sFFLHD.seq(D = D(), L = 3)
  sr <- reactive({ sFFLHD.seq(D = input$Dimension, L = input$Batch) })
  observeEvent(input$do, {
    s <- sr()
    #session$sendCustomMessage(type = 'testmessage',
    #                          message = 'Thank you for clicking')
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
}

shinyApp(ui, server)


