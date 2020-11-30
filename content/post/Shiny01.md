---
title: "Learn Shiny: Building a Shiny app"
date: 2020-11-27
lastmod: 2020-11-27
draft: false
tags: ["R", "Bioinformatics", "Shiny", "R Package"]
categories: ["R", "Bioinformatics", "R Package"]
author: "Xiaoni"

weight: 1

mathjax: true

# menu:
#   main:
#     parent: "docs"
#     weight: 1
---
 
Building a Shiny app - an interactive tutorial.

<!--more-->

## In summary

Learn to build a Shiny app step by step.

### 1. Prerequisites

- Install `shiny` package:

```r
Install.packages("shiny")
```

- To ensure successfully installed Shiny:

```r
library(shiny)
runExample(01_hello")
```

  ![hello](hello.png)

### 2. Shiny app basics

- ***UI***: a web page, showing the app to the user.
  - creating the layout of the app and telling Shiny exactly where things go.

- ***server***: a computer, powering the app.
  - local laptop.
  - server somewhere else.
  - for the logic of the app.
  - telling the web page what to show when the user interacts with the page.

### 3. Create an empty Shiny app

- All Shiny apps follow the same template:

```r
library(shiny)
ui <- fluidPage()
server <- function(input, output) {}
shinyApp(ui = ui, server = server)
```

#### 3.1 Alternate way to create a Shiny app

- separate UI file, ui.R.
- server file, server.R.
- Do not need to include the shinyApp(ui = ui, server = server) line.
- Put ui.R and server.R in the same folder.

#### 3.2 Fill out a Shiny app template using RStudio

- RStudio’s menu -> selecting File -> New File -> Shiny Web App
  
### 4. Load the dataset

- Dataset from [BC Liquor Store](bcl-data.csv), which has been cleaned up and simplified.

- Load and verify the dataset.
  - Place this line as the second line, just after `library(shiny)`.

```r
bcl <- read.csv("bcl-data.csv", stringsAsFactors = FALSE)

print(str(bcl))
```

### 5. Bulid the basic UI

- Add elements to the UI

#### 5.1 Add plain text to the UI

- render the text

```r
fluidPage("BC Liquor Store", "prices)
```

#### 5.2 Add formatted text and other HTML elements

- Try the `fluidPage()` function.

```r
ui <- fluidPage(
  h1("BC Liquor Store"), 
  "BC",
  "Liqqq",
  br(),
  "ssdse",
  strong("prices"), 
  div("this is blue", style = "color:blue;"),
  "Why it is boring")
```

  ![fig](fig1.png)

#### 5.3 Add a title

- Use `titlePanel()` function to set the "official" title of the web page.

```r
fluidPage(
  titlePanel("BC Liquor Store prices")
)
```

#### 5.4 Add a layout

- Use `sidebarLayout()` to add a simple structure.
  - a smaller sidebar and a larger main panel.

- Add the following code after the `titlePanel()`
  - Remember that all the arguments inside `fluidPage()` need to be separated by commas.

```r
sidebarLayout(
  sidebarPanel("our inputs will go here"),
  mainPanel("the results will go here")
)
```

#### 5.5 All UI functions are simply HTML wrappers

- Look at the output when printing the contents of the ui variable.

```r
print(ui)
# output
# <div class="container-fluid">
#   <h2>BC Liquor Store prices</h2>
#   <div class="row">
#     <div class="col-sm-4">
#       <form class="well">our inputs will go here</form>
#     </div>
#     <div class="col-sm-8">the results will go here</div>
#   </div>
# </div>
```

### 6. Add inputs to the UI

- Every input must have a unique inputId.

  - `textInput()` is used to let the user enter text, 
  - `numericInput()` lets the user select a number, 
  - `dateInput()` is for selecting a date, 
  - `selectInput()` is for creating a select box (aka a dropdown menu).
  
- The only way to find out what arguments you can use with a specific input function is to look at its help file.

#### 6.1 Input for price

```r
sliderInput("priceInput", "Price", min = 0, max = 100,
            value = c(25, 40), pre = "$")
```

- Place the code of `sliderInput()` inside `sidebarPanel()`, replace the text wrote earlier with this input("our inputs will go here").

#### 6.2 Input for product type

- Use radio buttons
  
```r
radioButtons("typeInput", "Product type",
            choices = c("BEER", "REFRESHMENT", "SPIRITS", "WINE"),
            selected = "WINE")
```

- Add this input code inside sidebarPanel(), after the previous input (separate them with a comma).

#### 6.3 Input for country

- Use select box.

```r
selectInput("countryInput", "Country",
            choices = c("CANADA", "FRANCE", "ITALY"))
```

- Continue to add this input code inside sidebarPanel(), after the previous input (separate them with a comma).

### 7. Add placeholders for outputs

#### 7.1 Output for a plot of the results

```r
plotOutput("coolplot")
```

#### 7.2 Output for a table summary of the results

- Create a UI element that will hold a table output:

```r
tableOutput("result")
```

### 8. Checkpoint: what our app looks like after implementing the UI

  ![fig2](fig2.png)

### 9. Implement server logic to create outputs

- For listening to changes to the inputs and creating outputs to show in the app.

- `input()`
  - a list you will read values ***from***. 
  - contain the values of all the different inputs at any given time
- `output()`
  - a list you will write values ***to***.
  - where you will save output objects (such as tables and plots) to display in your app.

#### 9.1 Building an output

- Three rules:
  - Save the output object into the output list (remember the app template - every server function has an output argument)
  - Build the object with a `render*` function, where `*` is the type of output
  - Access input values using the input list (every server function has an input argument)

- Create a plot and send it to the coolplot output.

```r
output$coolplot <- renderPlot({
    plot(rnorm(100))
  })
```

- Reload App and get output.
  
  ![fig3](fig3.png)

#### 9.2 Making an output react to an input

```r
output$coolplot <- renderPlot({
    plot(rnorm(input$priceInput[1]))
  })
```

- Reload App and get update.

  ![fig4](fig4.png)

#### 9.3 Building the plot output

- first, a histogram of the whole data, unfiltered.
  ![fig5](fig5.png)

- Use `dplyr()` function to filter the data and do histogram.
  ![fig6](fig6.png)
  ![fig7](fig7.png)
  ![fig8](fig8.png)

- At this point, get the code
  ![fig9](fig9.png)

#### 9.4 Building the table output

- Use `renderTable()`

### 10. Reactivity 101

```r
output$someoutput <- renderPlot({
  col <- input$mycolour
  num <- input$mynumber
  plot(rnorm(num), col = col)
})
```

The above render function accesses two different inputs: `input$mycolour` and `input$mynumber`. This means that this code block depends on both of these variables, so whenever either one of the two inputs is updated, the code gets re-executed with the new input values and `output$someoutput` is updated.

#### 10.1 Creating and accessing reactive variables

- `observe({})` and `reactive({})` functions in this section were just for learning purposes

#### 10.2 Using reactive variables to reduce code duplication

- As a reminder, Shiny creates ***a dependency tree*** with all the **reactive expressions** to know what value depends on what other value.

### 11. Using `uiOutput()` to create UI elements dynamically

- Any input that you normally create in the UI is created when the app starts, and it cannot be changed.

- In the server, use `uiOutput()` to be able to create an input dynamically.

#### 11.1 Basic example of `uiOutput()`

```r
library(shiny)
ui <- fluidPage(
  numericInput("num", "Maximum slider value", 5),
  uiOutput("slider")
)

server <- function(input, output) {
  output$slider <- renderUI({
    sliderInput("slider", "Slider", min = 0,
                max = input$num, value = 2.5)
  })
}

shinyApp(ui = ui, server = server)
```

  ![fig10](fig10.png)

#### 11.2 Use `uiOutput()` to populate the countries

  ![fig11](fig11.png)

#### 11.3 Errors showing up and quickly disappearing

when the app initializes, filtered is trying to access the country input, but the country input hasn’t been created yet. After Shiny finishes loading fully and the country input is generated, filtered tries accessing it again, this time it’s successful, and the error goes away.

#### 11.4 Use `uiOutput()` to add **subtypes**

- In the ui, add code:

```r
uiOutput("subtypeOutput"),
```

- In the server, modify `reactive({})` by adding code:

```r
Subtype == input$subtypeInput,
```

  ![fig12](fig12.png)

### 12. Final Shiny app code

- see [code](app.R).

### 13. Share apps with the world

#### 13.1 Host on ***shinyapps.io***

- [shinyapps.io](https://xiaonilee.shinyapps.io/shiny01/)

#### 13.2 Host on a ***Shiny Server***

- private Shiny serverc: DigitalOcean[(DO)](https://www.digitalocean.com/)

- DO droplet. droplet = your machine in the cloud

### 14. More Shiny features to check out

#### 14.1 Shiny in Rmarkdown

```r
---
output: html_document
runtime: shiny
---

```{r echo=FALSE}
sliderInput("num", "Choose a number",
            0, 100, 20)

renderPlot({
    plot(seq(input$num))
})
```

- Click `Knit` button

  ![figrmd](figrmd.png)

#### 14.2 Use `conditionalPanel()` to conditionally show UI elements

```r
library(shiny)
ui <- fluidPage(
  numericInput("num", "Number", 5, 1, 10),
  conditionalPanel(
    "input.num >=5",
    "Hello!"
  )
)
server <- function(input, output) {}
shinyApp(ui = ui, server = server)
```

  ![fig13](fig13.png)

#### 14.3 Use `navbarPage()` or `tabsetPanel()` to have multiple tabs in the UI

```r
library(shiny)
ui <- fluidPage(
  tabsetPanel(
    tabPanel("Tab 1", "Hello"),
    tabPanel("Tab 2", "there!")
  )
)
server <- function(input, output) {}
shinyApp(ui = ui, server = server)
```

  ![fig15](fig15.png)

#### 14.4 Use `DT` for beautiful, interactive tables

- In the `DT` package.
- Use `DT::dataTableOutput()` + `DT::renderDataTable()` to replace `tableOutput()` + `renderTable()`

#### 14.5 Use `isolate()` function to remove a dependency on a reactive variable

- To suppress some of reactive variables
- Cause a reactive variable to not be a dependency.
- With `isolate()` to wrap the code

#### 14.6 Use `update*Input()` functions to update input values programmatically

```r
library(shiny)
ui <- fluidPage(
  sliderInput("slider", "Move me", value = 5, 1, 10),
  numericInput("num", "Number", value = 5, 1, 10)
)
server <- function(input, output, session) {
  observe({
    updateNumericInput(session, "num", value = 5*(input$slider))
  })
}
shinyApp(ui = ui, server = server)
```

  ![fig16](fig16.png)

#### 14.7 Scoping rules in Shiny apps

#### 14.8 Use global.R to define objects available to both ui.R and server.R

#### 14.9 Add images

- Place an image under the "**www/**" folder
- Use the UI function `img(src = "image.png")`

#### 14.10 Add `JavaScript`/`CSS`

```r
library(shiny)
ui <- fluidPage(
  tags$head(tags$script("alert('Hello!');")),
  tags$head(tags$style("body{ color: blue; }")),
  "Hello"
)
server <- function(input, output) {
  
}
shinyApp(ui = ui, server = server)
```

### 15. Awesome add-on packages to Shiny

- shinyjs: Easily improve the user interaction and user experience in the Shiny apps in seconds.
- shinythemes: Easily alter app appearance.
- leaflet: Add interactive maps to apps.
- ggvis: Similar to ggplot2, but the plots are focused on being web-based and are more interactive.
- shinydashboard: Giving tools to create visual “dashboards”.
