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
#dat <- data.frame(ID = paste0("ID", 1:1010), A.1 = c(rep(0.5, 1000), abs(rnorm(10))), A.2 = c(rep(0.5, 1000), abs(rnorm(10))), B.1 = c(rep(0.5, 1000), abs(rnorm(10))), B.2 = c(rep(0.5, 1000), abs(rnorm(10))), C.1 = c(rep(0.5, 1000), abs(rnorm(10))), C.2 = c(rep(0.5, 1000), abs(rnorm(10))), C.3 = c(rep(0.5, 1000), abs(rnorm(10))), stringsAsFactors = FALSE)
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

shinyServer(function(input, output, session){
  observeEvent(c(input$goButton1, input$goButton2, input$goButton3), values$x <- values$x + 1)
  observeEvent(c(input$selPair1, input$selPair2, input$selPair3), values$x <- 0)
  observeEvent(c(input$selMetric1, input$selMetric2, input$selMetric3), values$x <- 0)
  observeEvent(c(input$selOrder1, input$selOrder2, input$selOrder3), values$x <- 0)
  observeEvent(input$binSize, values$x <- 0)
  
  observeEvent(input$selPair1, values$selPair <- input$selPair1)
  observeEvent(input$selPair2, values$selPair <- input$selPair2)
  observeEvent(input$selPair3, values$selPair <- input$selPair3)
  observeEvent(input$selMetric1, values$selMetric <- input$selMetric1)
  observeEvent(input$selMetric2, values$selMetric <- input$selMetric2)
  observeEvent(input$selMetric3, values$selMetric <- input$selMetric3)  
  observeEvent(input$selOrder1, values$selOrder <- input$selOrder1)
  observeEvent(input$selOrder2, values$selOrder <- input$selOrder2)
  observeEvent(input$selOrder3, values$selOrder <- input$selOrder3)
  
  observe({x <- values$selPair
    updateSelectizeInput(session, "selPair1", "Pairs:", choices = myPairs, options = list(maxItems = 2), selected = x)
    updateSelectizeInput(session, "selPair2", "Pairs:", choices = myPairs, options = list(maxItems = 2), selected = x)
    updateSelectizeInput(session, "selPair3", "Pairs:", choices = myPairs, options = list(maxItems = 2), selected = x)})

  observe({x <- values$selMetric
    updateSelectizeInput(session, "selMetric1", "Metric:", choices = myMetrics, selected = x)
    updateSelectizeInput(session, "selMetric2", "Metric:", choices = myMetrics, selected = x)
    updateSelectizeInput(session, "selMetric3", "Metric:", choices = myMetrics, selected = x)})
    
  observe({x <- values$selOrder
  updateSelectizeInput(session, "selOrder1", "Order:", choices = c("Increasing", "Decreasing"), selected = x)
  updateSelectizeInput(session, "selOrder2", "Order:", choices = c("Increasing", "Decreasing"), selected = x)
  updateSelectizeInput(session, "selOrder3", "Order:", choices = c("Increasing", "Decreasing"), selected = x)}) 
  
  # Create data subset based on two letters user chooses
  datSel <- eventReactive(values$selPair, {
    validate(need(length(values$selPair) == 2, "Select a pair."))
    sampleIndex <- reactive(which(sapply(colnames(dat), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(values$selPair[1], values$selPair[2])))
    dat[,c(1, sampleIndex())]
  }, ignoreNULL = FALSE)
  
  metricDF <- eventReactive(c(values$selPair, values$selMetric, values$selOrder), {
    metricDF <- metrics[[paste0(values$selPair[1], "vs", values$selPair[2])]]
    if (!is.null(metricDF[[values$selMetric]])){
      metricDF <- metricDF[order(metricDF[[values$selMetric]]),]
      if (values$selOrder == "Decreasing"){
        metricDF <- metricDF[order(-metricDF[[values$selMetric]]),]
      }
    }
    metricDF
  })
  
  currMetric <- eventReactive(values$x, {metricDF()[values$x, ]})
  currID <- eventReactive(currMetric(), {as.character(currMetric()$ID)})
  currGene <- eventReactive(currID(), {unname(unlist(datSel()[which(datSel()$ID == currID()), -1]))})
  output$info1 <- renderPrint({ currMetric() })
  output$info2 <- renderPrint({ currMetric() })
  output$info3 <- renderPrint({ currMetric() })
  
  output$hexPlot <- renderPlotly({
    
    sampleIndex1 <- which(sapply(colnames(datSel()), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(values$selPair[1]))
    sampleIndex2 <- which(sapply(colnames(datSel()), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(values$selPair[2]))
    
    minVal = min(datSel()[,-1])
    maxVal = max(datSel()[,-1])
    maxRange = c(minVal, maxVal)
    xbins= input$binSize
    buffer = (maxRange[2]-maxRange[1])/(xbins/2)
    x <- c()
    y <- c()
    for (i in 1:length(sampleIndex1)){
      for (j in 1:length(sampleIndex2)){
        x <- c(x, unlist(datSel()[,(sampleIndex1[i])]))
        y <- c(y, unlist(datSel()[,(sampleIndex2[j])]))
      }
    }
    
    h <- hexbin(x=x, y=y, xbins=xbins, shape=1, IDs=TRUE, xbnds=maxRange, ybnds=maxRange)
    hexdf <- data.frame (hcell2xy (h),  hexID = h@cell, counts = h@count)
    attr(hexdf, "cID") <- h@cID
    
    my_breaks <- c(2, 4, 6, 8, 20, 1000)    
    p <- reactive(ggplot(hexdf, aes(x=x, y=y, fill = counts, hexID=hexID)) + geom_hex(stat="identity") + geom_abline(intercept = 0, color = "red", size = 0.25) + labs(x = values$selPair[1], y = values$selPair[2]) + coord_fixed(xlim = c(-0.5, (maxRange[2]+buffer)), ylim = c(-0.5, (maxRange[2]+buffer))) + theme(aspect.ratio=1) + scale_fill_gradient(name = "count", trans = "log", breaks = my_breaks, labels = my_breaks, guide="legend"))
    
    plotlyHex <- reactive(ggplotly(p(), height = 400, width = 400)) #height = 400, width = 400
    
    # Use onRender() function to draw x and y values of selected row as orange point
    plotlyHex() %>% onRender("
     function(el, x, data) {
     noPoint = x.data.length;
     Shiny.addCustomMessageHandler('points', function(drawPoints) {
     if (x.data.length > noPoint){
     Plotly.deleteTraces(el.id, x.data.length-1);
     }
     var Traces = [];
     var trace = {
     x: drawPoints.geneX,
     y: drawPoints.geneY,
     mode: 'markers',
     marker: {
     color: 'orange',
     size: 7
     },
     hoverinfo: 'none'
     };
     Traces.push(trace);
     Plotly.addTraces(el.id, Traces);
     });}")
  })
  
  
  output$scatterPlot <- renderPlotly({
    
    sampleIndex1 <- which(sapply(colnames(datSel()), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(values$selPair[1]))
    sampleIndex2 <- which(sapply(colnames(datSel()), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(values$selPair[2]))
    
    minVal = min(datSel()[,-1])
    maxVal = max(datSel()[,-1])
    maxRange = c(minVal, maxVal)
    xbins= input$binSize
    buffer = (maxRange[2]-maxRange[1])/(xbins/2)
    x <- c()
    y <- c()
    for (i in 1:length(sampleIndex1)){
      for (j in 1:length(sampleIndex2)){
        x <- c(x, unlist(datSel()[,(sampleIndex1[i])]))
        y <- c(y, unlist(datSel()[,(sampleIndex2[j])]))
      }
    }
    tempDF <- data.frame(x=x, y=y)
    
    p2 <- reactive(ggplot(tempDF, aes(x=x, y=y)) + geom_point(alpha = input$alpha) + geom_abline(intercept = 0, color = "red", size = 0.25) + labs(x = values$selPair[1], y = values$selPair[2]) + coord_fixed(xlim = c(-0.5, (maxRange[2]+buffer)), ylim = c(-0.5, (maxRange[2]+buffer))) + theme(aspect.ratio=1))
    
    plotlyScatter <- reactive(ggplotly(p2(), height = 400, width = 400))
    
    # Use onRender() function to draw x and y values of selected row as orange point
    plotlyScatter() %>% onRender("
     function(el, x, data) {
     noPoint = x.data.length;
     Shiny.addCustomMessageHandler('points', function(drawPoints) {
     
     if (x.data.length > noPoint){
      Plotly.deleteTraces(el.id, noPoint);
     }

     var Traces = [];
     var trace = {
     x: drawPoints.geneX,
     y: drawPoints.geneY,
     mode: 'markers',
     marker: {
     color: 'orange',
     size: 7
     },
     hoverinfo: 'none'
     };
     Traces.push(trace);
     Plotly.addTraces(el.id, Traces);
     });}")
  })
  
  observe({
    #print(values$x)
    # Get x and y values of selected row
    sampleIndex1 <- which(sapply(colnames(datSel()), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(values$selPair[1]))
    sampleIndex2 <- which(sapply(colnames(datSel()), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(values$selPair[2]))
    
    geneX <- c()
    geneY <- c()
    for (i in 1:length(sampleIndex1)){
      for (j in 1:length(sampleIndex2)){
        geneX <- c(geneX, currGene()[sampleIndex1[i]-1])
        geneY <- c(geneY, currGene()[sampleIndex2[j]-1])
      }
    }
    # Send x and y values of selected row into onRender() function
    session$sendCustomMessage(type = "points", message=list(geneX=geneX, geneY=geneY))
  })
  
  output$boxPlot <- renderPlotly({
    nVar = reactive(ncol(datSel()))
    colNms <- reactive(colnames(datSel()[, c(2:nVar())]))
    
    boxDat <- eventReactive(datSel(), {
      boxDat <- datSel()[, c(1:nVar())] %>% gather(key, val, -c(ID))
      colnames(boxDat) <- c("ID", "Sample", "Count")
      boxDat
    })
    
    BP <- reactive(ggplot(boxDat(), aes(x = Sample, y = Count)) + geom_boxplot())
    ggBP <- reactive(ggplotly(BP(), width=600))
    
    observe({
      session$sendCustomMessage(type = "lines", currGene())
    })
    
    ggBP() %>% onRender("
      function(el, x, data) {
      
console.log(['x', x])
console.log(['x.data.length', x.data.length])
      noPoint = x.data.length;
      
      function range(start, stop, step){
      var a=[start], b=start;
      while(b<stop){b+=step;a.push(b)}
      return a;
      };
      
      Shiny.addCustomMessageHandler('lines',
      function(drawLines) {

console.log(['el.id', el.id])
console.log(['noPoint', noPoint])
console.log(['x', x])

i = x.data.length
if (i > 1){
while (i > 1){
Plotly.deleteTraces(el.id, (i-1));
i--;
}
      }
      
      var dLength = drawLines.length
      
      var Traces = [];
      var traceLine = {
      x: range(1, dLength, 1),
      y: drawLines,
      mode: 'lines',
      line: {
      color: 'orange',
      width: 1
      },
      opacity: 0.9,
      }
      Traces.push(traceLine);
      Plotly.addTraces(el.id, Traces);
      })
      }")})
})
