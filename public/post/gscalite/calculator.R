
# calculator
library(shiny)
ui <- fluidPage(textInput('a', 'please input the first number:'),
                checkboxGroupInput('op', 'please choose a sign',
                                   c('+' = '+',
                                     '-' = '-',
                                     '*' = '*',
                                     '/' = '/',
                                     '%' = '%%')),
                textInput('b', 'please input the second number:'),
                submitButton('do', 'calculate'),
                textOutput('result')
)

server <- function(input, output) {
  output$result <- renderPrint({
    input$do
    
    print(eval(parse(text = paste0(input$a, input$op, input$b))))
  })
}
shinyApp(ui=ui, server = server)
