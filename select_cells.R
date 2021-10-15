require(shiny)
require(plotly)
library(ggplot2)
require(viridis)
ui <- fluidPage(
  plotlyOutput("plot"),
  verbatimTextOutput("brush")
)
server <- function(input, output, session) {
  output$plot <- renderPlotly({
    # use the key aesthetic/argument to help uniquely identify selected observations
    key <- row.names(pd)
    p <- ggplot(pd, aes(x, y, color = mito)) +
      geom_point(size = 0.5, alpha=1) +
      scale_color_viridis(option = "D", direction = -1) +
      theme_void()
    ggplotly(p,height = 800, width = 900) %>% toWebGL() %>%
      highlight("plotly_selected", dynamic = F )
  })
  output$brush <- renderPrint({
    d <- event_data("plotly_selected")
    req(d)
    sel = dplyr::inner_join(d,pd)
    print(sel)
    write(sel$cell,
          paste0("data/OB",m,"_slideseq/Slide",s,"clean.txt"))
  })
}