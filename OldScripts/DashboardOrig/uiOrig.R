
library(shinydashboard)

sidebar <- dashboardSidebar(
  hr(),
  sidebarMenu(id="tabs",
    menuItem("Plot", tabName="plot", icon=icon("line-chart"), selected=TRUE),
    menuItem("Table", tabName = "table", icon=icon("table")),
    menuItem("Codes",  icon = icon("file-text-o"),
     menuSubItem("Mlxtran", tabName = "pkmodel", icon = icon("angle-right")),
     menuSubItem("ui.R", tabName = "ui", icon = icon("angle-right")),
     menuSubItem("server.R", tabName = "server", icon = icon("angle-right"))
    ),
    menuItem("ReadMe", tabName = "readme", icon=icon("mortar-board")),
    menuItem("About", tabName = "about", icon = icon("question"))
  ),
  hr(),
  conditionalPanel("input.tabs=='plot'",
   fluidRow(
     column(1),
     column(10,
      checkboxInput("first", "First order", TRUE),
      checkboxInput("zero", "Zero order", TRUE),
      checkboxInput("al", "alpha order", FALSE),
      checkboxInput("sequential", "Sequential (0-1)", FALSE),
      checkboxInput("mixed", "Simultaneous (0-1)", FALSE),
      checkboxInput("saturated", "Saturated", FALSE),
      checkboxInput("legend", "Legend", TRUE)
     )
   )
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "readme",
      withMathJax(), 
      includeMarkdown("readMe.Rmd")
    ),
    tabItem(tabName = "plot",
      fluidRow(
        column(width = 4, 
         tabBox( width = NULL,
           tabPanel(h5("parameters"),
            conditionalPanel(condition="input.sequential=='1' | input.mixed=='1' | input.first=='1' | input.alpha=='1' ", sliderInput("ka", "ka:", value = 0.5, min = 0.1, max = 3, step=0.1)),
            conditionalPanel(condition="input.sequential=='1' | input.mixed=='1' | input.zero=='1' ", sliderInput("Tk0", "Tk0:", value = 5, min = 0, max = 10, step=0.5)),
            conditionalPanel(condition="input.al=='1'", sliderInput("alpha", "alpha:", value = 0.5, min = 0, max = 2, step=0.1)),
            conditionalPanel(condition="input.sequential=='1' | input.mixed=='1' ", sliderInput("F0", "F0:", value = 0.5, min = 0, max = 1, step=0.1)),
            conditionalPanel(condition="input.saturated=='1'",
              sliderInput("Vm", "Vm:", value = 0.9, min = 0, max = 2, step=0.1),
              sliderInput("Km", "Km:", value = 0.2, min = 0, max = 1, step=0.1)),
            sliderInput("k", "k:", value = 0.1, min = 0, max = 2, step=0.05)),
           tabPanel(h5("dosage"),
              sliderInput("tfd", "Time of first dose:", value=0, min=0, max = 20, step=1),
              sliderInput("nd", "Number of doses:", value=1, min=0, max = 10, step=1),
              sliderInput("ii", "Interdose interval:", value = 9, min = 0.5, max = 15, step=0.5),
              sliderInput("amt", "Amount:", value = 5, min = 0, max = 20, step=1))
         )),
        column(width = 8,
         box(width = NULL, plotOutput("plot",height="500px"), collapsible = TRUE, title = "Plot", status = "primary", solidHeader = TRUE)))
    ),
    tabItem(tabName = "table",
      box( width = NULL, status = "primary", solidHeader = TRUE, title="Table",                
      downloadButton('downloadTable', 'Download'),
      br(),br(),
      tableOutput("table"))),
    tabItem(tabName = "pkmodel",
            box( width = NULL, status = "primary", solidHeader = TRUE, title="absorptionModel.txt",             downloadButton('downloadData1', 'Download'),
            br(),br(),
            pre(includeText("absorptionModel.txt"))
            )
    ),
    tabItem(tabName = "ui",
      box( width = NULL, status = "primary", solidHeader = TRUE, title="ui.R",
        downloadButton('downloadData2', 'Download'),
        br(),br(),
        pre(includeText("ui.R")))),
    tabItem(tabName = "server",
      box( width = NULL, status = "primary", solidHeader = TRUE, title="server.R",
        downloadButton('downloadData3', 'Download'),
        br(),br(),
        pre(includeText("server.R")))),
    tabItem(tabName = "about", includeMarkdown("../../about/about.Rmd")
    )
  )
)

dashboardPage(
  dashboardHeader(title = "Absorption processes"),
  sidebar,
  body
)