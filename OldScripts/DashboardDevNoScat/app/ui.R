library(shinydashboard)
library(shiny)
library(plotly)
library(hexbin)
library(htmlwidgets)
library(dplyr)
library(tidyr)

set.seed(1)
# Create data and subsets of data based on user selection of pairs
dat <- data.frame(ID = paste0("ID", 1:100), A.1 = c(1, abs(rnorm(99))), A.2 = c(1, abs(rnorm(99))), B.1 = c(1, abs(rnorm(99))), B.2 = c(1, abs(rnorm(99))), C.1 = c(1, abs(rnorm(99))), C.2 = c(1, abs(rnorm(99))), C.3 = c(1, abs(rnorm(99))), stringsAsFactors = FALSE)
datCol <- colnames(dat)[-which(colnames(dat) %in% "ID")]
myPairs <- unique(sapply(datCol, function(x) unlist(strsplit(x,"[.]"))[1]))
metrics <- list()
for (i in 1:(length(myPairs)-1)){
  for (j in (i+1):length(myPairs)){
    metrics[[paste0(myPairs[i],"vs",myPairs[j])]] <- data.frame(ID = paste0("ID", 1:100), pValAdj = runif(100, 0, 1), logFC = runif(100, 0, 6), AveExp = runif(100, 0, 60))
  }
}
myMetrics <- colnames(metrics[[1]])[-which(colnames(metrics[[1]]) %in% "ID")]
values <- reactiveValues(x=0, selPair=NULL, selMetric=NULL, selOrder=NULL)

sidebar <- dashboardSidebar(
  width = 180,
  hr(),
  sidebarMenu(id="tabs",
    menuItem("Binned scatterplot", tabName="hexPlot", selected=TRUE), # icon=icon("line-chart"),
    menuItem("Parallel coordinates", tabName = "boxPlot") # icon = icon("file-text-o")
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "hexPlot",
      fluidRow(
        column(width = 4, 
         box(width = NULL, status = "primary", title = "Plot metrics", solidHeader = TRUE,
              selectizeInput("selPair1", "Pairs:", choices = myPairs, multiple = TRUE, options = list(maxItems = 2)),
              selectInput("selMetric1", "Metric:", choices = myMetrics),
              selectInput("selOrder1", "Order:", choices = c("Increasing", "Decreasing")),
              numericInput("binSize", "Hexagon size:", value = 10),
              actionButton("goButton1", "Plot case!"))),
        column(width = 8,
          box(width = NULL, plotlyOutput("hexPlot"), collapsible = FALSE, background = "black", title = "Binned scatterplot", status = "primary", solidHeader = TRUE))),
      fluidRow(
        column(width = 8, offset = 4,
          box(width = NULL, verbatimTextOutput("info1"), collapsible = TRUE, title = "Observation metrics", status = "primary", solidHeader = TRUE)))),

    tabItem(tabName = "boxPlot",
      fluidRow(
        column(width = 4, 
           box(width = NULL, status = "primary", title = "Plot metrics", solidHeader = TRUE,
             selectizeInput("selPair2", "Pairs:", choices = myPairs, multiple = TRUE, options = list(maxItems = 2)),
             selectInput("selMetric2", "Metric:", choices = myMetrics),
             selectInput("selOrder2", "Order:", choices = c("Increasing", "Decreasing")),
             actionButton("goButton2", "Plot case!"))),
        column(width = 8,
           box(width = NULL, plotlyOutput("boxPlot"), collapsible = FALSE, background = "black", title = "Parallel coordinate plot", status = "primary", solidHeader = TRUE))),
      fluidRow(
        column(width = 8, offset = 4,
          box(width = NULL, verbatimTextOutput("info2"), collapsible = TRUE, title = "Observation metrics", status = "primary", solidHeader = TRUE))))
  )
)

dashboardPage(
  dashboardHeader(title = "Overlaying cases", titleWidth = 180),
  sidebar,
  body
)