library(shiny)
library(ggplot2)
library(dplyr)
bcl <- read.csv("bcl-data.csv", stringsAsFactors = FALSE)

ui <- fluidPage(
  titlePanel("BC Liquor Store prices"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("priceInput", "Price", min = 0, max = 100,
                  value = c(25, 40), pre = "$"),
      radioButtons("typeInput", "Product type",
                   choices = c("BEER", "REFRESHMENT", "SPIRITS", "WINE"),
                   selected = "WINE"),
      # selectInput("countryInput", "Country",
      #             choices = c("CANADA", "FRANCE", "ITALY")),
      uiOutput("subtypeOutput"),
      uiOutput("countryOutput")
      # textInput("txt", "num of results found:",
      #           value = 100)
    ),
    mainPanel(
      plotOutput("coolplot"),
      br(), br(),
      tableOutput("results")
      # br(), br(),
      # textOutput("text")
    )
  )
)
# print(ui)
server <- function(input, output) {
  # output$coolplot <- renderPlot({
  #   #print(input$priceInput)
  #   #plot(rnorm(input$priceInput[1]))
  #   ggplot(bcl, aes(Alcohol_Content)) +
  #   geom_histogram()
  # })
  filtered <- reactive({
    if (is.null(input$countryInput)) {
      return(NULL)
    }
    
    bcl %>%
      filter(Price >= input$priceInput[1],
             Price <= input$priceInput[2],
             Type == input$typeInput,
             Subtype == input$subtypeInput,
             Country == input$countryInput
      )
  })
  
  output$coolplot <- renderPlot({
    if (is.null(filtered())) {
      return()
    }
    
    ggplot(filtered(), aes(Alcohol_Content)) +
      geom_histogram()
  })
  output$results <- renderTable({
    filtered()
  })
  # output$text <- renderText({ input$txt })
  # print(input$priceInput)--will get error, need to use observe({})
#   observe({print(input$priceInput)})
#   observe({cat(input$x)})
#   priceDiff <- reactive({
#     diff(input$priceInput)
#   })
#   observe({ print(priceDiff()) }) #use priceDiff() rather than priceDiff
  output$subtypeOutput <- renderUI({
    selectInput("subtypeInput", "Subtype",
                sort(unique(bcl$Subtype)),
                selected = "ICE WINE WHITE")
  })
  
  output$countryOutput <- renderUI({
    selectInput("countryInput", "Country",
                sort(unique(bcl$Country)),
                selected = "CANADA")
  })
  
  
}
shinyApp(ui = ui, server = server)