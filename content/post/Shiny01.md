---
title: "Learn Shiny: Building Shiny apps"
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
 
Radiant - A shiny app for statistics and machine learning.

<!--more-->

## In summary

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

#### Alternate way to create a Shiny app

- separate UI file, ui.R.
- server file, server.R.
- Do not need to include the shinyApp(ui = ui, server = server) line.
- Put ui.R and server.R in the same folder.

#### Let RStudio fill out a Shiny app template for you

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

### 11. Using uiOutput() to create UI elements dynamically

- Any input that you normally create in the UI is created when the app starts, and it cannot be changed.

- In the server, use `uiOutput()` to be able to create an input dynamically.

#### 11.1 Basic example of uiOutput()
