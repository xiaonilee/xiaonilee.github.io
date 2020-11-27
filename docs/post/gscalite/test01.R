rm(list = ls())
options(stringsAsFactors = FALSE)


# Shiny app template
library(shiny)
ui <- fluidPage(h1('Hello, Shiny'), p(h3('Learn RShiny')))
server <- function(input, output) {}
shinyApp(ui=ui, server = server)
