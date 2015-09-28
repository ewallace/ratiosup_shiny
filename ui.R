library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
    
    # Application title
    titlePanel("Ratio in supernatant for yeast genes after 42C, 10min heat shock"),
    
    # Sidebar with a slider input for the number of bins
    sidebarLayout(
        sidebarPanel(
            textInput("ids",
                      "Enter gene identifiers separated by commas:",
                      value = "PGK1,DED1"),
            selectInput("idType", "identify by", c("gene","orf"), selected="gene"),
            selectInput("replicate", "replicate", choices=c(1,2), selected=1),
            selectInput("maxtime", "maximum time", choices=c(60,180), selected=60)
        ),
        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("plot") # ,plotOutput("tempPlot")
        )
    ),
    
    # Subtitle with paper reference
    mainPanel(
        p("Data from: 
          Reversible, Specific, Active Aggregates of Endogenous Proteins Assemble upon Heat Stress,
          Wallace et al., Cell 162 (6), 2015, http://drummondlab.org/endogenous-aggregates")
    )
))