library(ggplot2)
require(viridis)
ui <- basicPage(
  plotOutput("plot1", click = "plot_click"),
  verbatimTextOutput("info")
)

server <- function(input, output) {
  output$plot1 <- renderPlot({
    ggplot(df, aes(x, y, color = percent.mt, label = X)) +
      geom_point(size = 0.5, alpha=1) +
      scale_color_viridis(option = "D", direction = -1) +
      coord_fixed()
  })
  
  output$info <- renderPrint({
    # With ggplot2, no need to tell it what the x and y variables are.
    # threshold: set max distance, in pixels
    # maxpoints: maximum number of rows to return
    # addDist: add column with distance, in pixels
    nearPoints(df, input$plot_click, threshold = 10, maxpoints = 1,
               addDist = TRUE)
  })
}

shinyApp(ui, server)