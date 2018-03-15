#' Plot interactive scatterplot matrices
#' 
#' Plot interactive scatterplot matrices.
#' 
##' @details There are four options:
##' \itemize{
##'  \item{"hexagon": }{Plot interactive scatterplot matrix with hexagon binning}
##'  \item{"foldChange": }{Plot interactive scatterplot matrix with fold change}
##'  \item{"orthogonal": }{Plot interactive scatterplot matrix with orthogonal distance}
##'  \item{"prediction": }{Plot interactive scatterplot matrix with prediction interval}
##' } 
#' 
#' @param data data frame containing read counts
#' @param outDir output directory to save all images (default current directory)
#' @param threshOrth threshold of orthogonal distance (default 3; used in "orthogonal")
#' @param xbins the number of bins partitioning the range of the plot (default 10; used in "hexagon")
#' @param option the type of plot (can choose from c("hexagon", "foldChange", "orthogonal", "prediction"); default "hexagon")
#' @importFrom plotly plotlyOutput ggplotly renderPlotly config
#' @importFrom ggplot2 ggplot aes_string aes xlim ylim geom_boxplot
#' @importFrom shiny verbatimTextOutput fluidPage reactive renderPrint shinyUI sliderInput shinyServer shinyApp HTML
#' @importFrom htmlwidgets onRender
#' @importFrom utils str
#' @importFrom tidyr gather
#' @importFrom stats qt lm coef
#' @importFrom hexbin hexbin hcell2xy
#' @importFrom stringr str_replace str_trim
#' @export
#' @examples
#' data(soybean_cn)
#' soybean_cn <- soybean_cn[,1:7]
#' plotScatterInteractive(soybean_cn)
plotScatterInteractive = function(data=data, outDir=getwd(), threshOrth=3, xbins=10, option="hexagon"){
  
  if (option=="foldChange"){
    scatMatFCPCP(data=data)
  }
  else if (option=="orthogonal"){
    scatMatOrthPCP(data=data, threshOrth=threshOrth)
  }
  else if (option=="prediction"){
    scatMatPIPCP(data=data)
  }  
  else{
    scatMatHexPCP(data=data, xbins=xbins)
  }
} 

##############################################################################
##############################################################################
##############################################################################
# HELPER FUNCTIONS

scatMatFCPCP = function(data=data){
  
  key <- val <- ID <- NULL
  
  ui <- shinyUI(fluidPage(
    sliderInput("threshold", "Fold Change:", min = 0, max = 30, value=15, step=0.1),
    plotlyOutput("myPlot", height = 620),
    plotlyOutput("boxPlot"),
    verbatimTextOutput("selectedValues")
  ))
  
  server <- shinyServer(function(input, output) {
    
    minVal = min(data[,-1])
    maxVal = max(data[,-1])
    # Designate end points of lines to be drawn
    minLine = minVal
    maxLine = maxVal
    maxRange = c(minVal, maxVal)
    buffer = maxRange[2]/xbins
    
    my_fn <- function(data, mapping, ...){
      x = data[,c(as.character(mapping$x))]
      y = data[,c(as.character(mapping$y))]
      p <- ggplot(data = data, aes(x=x, y=y)) + geom_point(alpha=0) + coord_cartesian(xlim = c(minVal, maxVal+buffer), ylim = c(minVal, maxVal+buffer))
      p
    }
    
    p <- ggpairs(data[1,-1], lower = list(continuous = my_fn))
    ggPS <- ggplotly(p)
    
    myLength <- length(ggPS[["x"]][["data"]])
    for (i in 1:myLength){
      item =ggPS[["x"]][["data"]][[i]]$text[1]
      if (!is.null(item))
        if (!startsWith(item, "co")){
          ggPS[["x"]][["data"]][[i]]$hoverinfo <- "none"
        }
    }
    
    output$myPlot <- renderPlotly(ggPS %>%
      onRender("
           function(el, x, data) {
           len = Math.sqrt(document.getElementsByClassName('cartesianlayer')[0].childNodes.length);
  console.log('Test Console')
           AxisNames = [];
           for (i = 1; i < (len+1); i++) {
           AxisNames.push(document.getElementsByClassName('infolayer')[0].childNodes[i].textContent);
           }
           noPoint = x.data.length;
           var SubPoints = [];
           var Traces = [];
           var i=0;
           var k=1;
           var inc = (data.maxLine-data.minLine)/100
           var xv = [];
           var uyv = [];
           var lyv = [];
           var a = data.minLine
           while (a < data.maxLine){
           var fract = a
           xv.push(a);
           uyv.push(a*(data.val+1));
           lyv.push(a/(data.val+1));
           a+=inc;
           }
           
           console.log(xv)
           console.log(uyv)
           console.log(lyv)
           
           while ((i*len+k)<=Math.pow((len-1),2)) {
           while ((i+k)<len){
           var selRows = [];
           data.dat.forEach(function(row){
           var fract = row[AxisNames[i]] / row[AxisNames[(len-k)]]
           if(fract > (data.val + 1) || fract < (1/(data.val+1))){
           selRows.push(row);
           }})
           var xArr = [];
           for (a=0; a<selRows.length; a++){
           xArr.push(selRows[a][AxisNames[i]])
           }
           var yArr = [];
           for (a=0; a<selRows.length; a++){
           yArr.push(selRows[a][AxisNames[(len-k)]])
           }
           var keepIndex = [];
           for (a=0; a<selRows.length; a++){
           keepIndex.push(selRows[a]['ID'])
           }
           SubPoints.push(keepIndex);
           var tracePoints = {
           x: xArr,
           y: yArr,
           hoverinfo: 'none',
           mode: 'markers',
           marker: {
           color: 'black',
           size: 4
           },
           xaxis: 'x' + (i+1),
           yaxis: 'y' + (i*len+k)
           };
           
           var hiLine = {
           x: xv,
           y: uyv,
           mode: 'lines',
           line: {
           color: 'gray',
           width: 1
           },
           xaxis: 'x' + (i+1),
           yaxis: 'y' + (i*len+k),
           opacity: 0.25,
           hoverinfo: 'none'
           };
           
           var lowLine = {
           x: xv,
           y: lyv,
           mode: 'lines',
           fill: 'tonexty',
           line: {
           color: 'gray',
           width: 1
           },
           xaxis: 'x' + (i+1),
           yaxis: 'y' + (i*len+k),
           opacity: 0.25,
           hoverinfo: 'none'
           };
           
           Traces.push(tracePoints);
           Traces.push(hiLine);
           Traces.push(lowLine);
           k++;
           }
           i++;
           k=1;
           }
           Plotly.addTraces(el.id, Traces);
                                         
           var idRows = []
           for (a=0; a<data.dat.length; a++){
           idRows.push(data.dat[a]['ID'])
           }
           
           
           el.on('plotly_selected', function(e) {
           numSel = e.points.length
           cN = e.points[0].curveNumber;
           
           var pointNumbers = [];
           for (a=0; a<numSel; a++){
           pointNumbers.push(e.points[a].pointNumber)
           }
           
           // Determine which subplot was selected
           subPlot = (cN - Math.pow(len,2))/3+1
           
           var selData = []
           for (a=0; a<pointNumbers.length; a++){
           selData.push(data.dat[idRows.indexOf(SubPoints[subPlot-1][pointNumbers[a]])])
           }
           
           var selData = []
           var selDots = []
           for (a=0; a<pointNumbers.length; a++){
           var selDot = SubPoints[subPlot-1][pointNumbers[a]]
           selData.push(data.dat[idRows.indexOf(selDot)])
           selDots.push(idRows.indexOf(selDot))
           }
           
           Shiny.onInputChange('selDots', selDots);
           
           var Traces = [];
           var i=0;
           var k=1;
           while ((i*len+k)<=Math.pow((len-1),2)) {
           var xArr = [];
           for (a=0; a<selData.length; a++){
           xArr.push(selData[a][AxisNames[i]])
           }
           while ((i+k)<len){
           var yArr = [];
           for (a=0; a<selData.length; a++){
           yArr.push(selData[a][AxisNames[(len-k)]])
           }
           var trace = {
           x: xArr,
           y: yArr,
           mode: 'markers',
           marker: {
           color: 'red',
           size: 4
           },
           xaxis: 'x' + (i+1),
           yaxis: 'y' + (i*len+k),
           hoverinfo: 'none'
           };
           Traces.push(trace);
           k++;
           }
           i++;
           k=1;
           }
           Plotly.addTraces(el.id, Traces);
           })}", data = list(dat=data, val = input$threshold, minLine=minLine, maxLine=maxLine)))
  
    selID <- reactive(input$selDots)
    
    pcpDat <- reactive(data[selID()+1, ])
    output$selectedValues <- renderPrint({str(pcpDat())})
    
    colNms <- colnames(data[, c(2:(ncol(data)))])
    
    boxDat <- data %>% gather(key, val, -c(ID))
    BP <- ggplot(boxDat, aes(x = key, y = val)) + geom_boxplot()
    ggBP <- ggplotly(BP)
    
    output$boxPlot <- renderPlotly({
      ggBP %>% onRender("
                        function(el, x, data) {
                        
                        var Traces = [];
                        
                        var dLength = data.pcpDat.length
                        var vLength = data.nVar
                        var cNames = data.colNms
                        
                        for (a=0; a<dLength; a++){
                        xArr = [];
                        yArr = [];
                        for (b=0; b<vLength; b++){
                        xArr.push(b+1)
                        yArr.push(data.pcpDat[a][cNames[b]]);
                        }
                        
                        var traceHiLine = {
                        x: xArr,
                        y: yArr,
                        mode: 'lines',
                        line: {
                        color: 'orange',
                        width: 1
                        },
                        opacity: 0.9,
                        }
                        Traces.push(traceHiLine);
                        }
                        Plotly.addTraces(el.id, Traces);
                        }", data = list(pcpDat = pcpDat(), nVar = p$nrow, colNms = colNms))})
  })

shinyApp(ui, server)
}



scatMatHexPCP = function(data=data, xbins=xbins){
  
  counts <- hexID <- key <- val <- ID <- NULL
  
  ui <- shinyUI(fluidPage(
    plotlyOutput("scatMatPlot", width = 700, height = 650),
    plotlyOutput("boxPlot", width = 700),
    HTML("<br><br><br>"),
    print("Selected Gene IDs:"),
    verbatimTextOutput("selectedValues")
  ))
  
  server <- shinyServer(function(input, output) {
    
    ################################ Prepare scatterplot matrix
    ###########################################################
    
    maxVal = max(abs(data[,-1]))
    maxRange = c(-1*maxVal, maxVal)

    my_fn <- function(data, mapping, ...){
      x = data[,c(as.character(mapping$x))]
      y = data[,c(as.character(mapping$y))]
      h <- hexbin(x=x, y=y, xbins=xbins, shape=1, IDs=TRUE, xbnds=maxRange, ybnds=maxRange)
      hexdf <- data.frame (hcell2xy (h),  hexID = h@cell, counts = h@count)
      attr(hexdf, "cID") <- h@cID
      p <- ggplot(hexdf, aes(x=x, y=y, fill = counts, hexID=hexID)) + geom_hex(stat="identity") + geom_abline(intercept = 0, color = "red", size = 0.25) + coord_cartesian(xlim = c(-0.5, maxRange[2]+0.5), ylim = c(-0.5, maxRange[2]+0.5)) + theme(legend.position="none")
      p
    }
    
    p <- ggpairs(data[,-1], lower = list(continuous = my_fn))
    pS <- p
    
    ggPS <- ggplotly(pS, width=700, height=600)
    
    myLength <- length(ggPS[["x"]][["data"]])
    for (i in 1:myLength){
      item =ggPS[["x"]][["data"]][[i]]$text[1]
      if (!is.null(item)){
        if (!startsWith(item, "co")){
          ggPS[["x"]][["data"]][[i]]$hoverinfo <- "none"
        }}
      hexHover = ggPS[["x"]][["data"]][[i]]$text
      if (!is.null(hexHover) && grepl("hexID", hexHover)){
        ggPS[["x"]][["data"]][[i]]$text <- strsplit(hexHover, "<")[[1]][1]
        ggPS[["x"]][["data"]][[i]]$t2 <- hexHover
        ggPS[["x"]][["data"]][[i]]$hoverinfo <- "text"
      }
    }
    
    for(i in 2:(p$nrow)) {
      for(j in 1:(p$nrow-1)) {
        data[[paste(i,j,sep="-")]] <- attr(pS[i,j]$data, "cID")
      }
    }
    
    output$scatMatPlot <- renderPlotly({
      ggPS %>% onRender("
        function(el, x, data) {
        
        function range(start, stop, step){
        var a=[start], b=start;
        while(b<stop){b+=step;a.push(b)}
        return a;
        };
        
        len = Math.sqrt(document.getElementsByClassName('cartesianlayer')[0].childNodes.length);
        AxisNames = [];
        for (i = 1; i < (len+1); i++) {
        AxisNames.push(Object.keys(data[0])[i]);
        }
        noPoint = x.data.length;
        
        el.on('plotly_click', function(e) {
        
        if (x.data.length > noPoint){
        Plotly.deleteTraces(el.id, range(noPoint, (noPoint+(len*(len-1)/2-1)), 1));
        }
        
        xVar = (e.points[0].xaxis._id).replace(/[^0-9]/g,'')
        if (xVar.length == 0) xVar = 1
        yVar = (e.points[0].yaxis._id).replace(/[^0-9]/g,'')
        if (yVar.length == 0) yVar = 1
        myX = len + 1 - (yVar - len * (xVar - 1))
        myY = xVar
        cN = e.points[0].curveNumber
        split1 = (x.data[cN].text).split(' ')
        hexID = (x.data[cN].t2).split(':')[2]
        counts = split1[1].split('<')[0]
        var selRows = [];
        data.forEach(function(row){
        if(row[myX+'-'+myY]==hexID) selRows.push(row);
        });
        selID = []
        for (a=0; a<selRows.length; a++){
        selID.push(selRows[a]['ID'])
        }
        // Save selected row IDs for PCP
        Shiny.onInputChange('selID', selID);
        
        var Traces = [];
        var i=0;
        var k=1;
        while ((i*len+k)<=Math.pow((len-1),2)) {
        var xArr = [];
        for (a=0; a<selRows.length; a++){
        xArr.push(selRows[a][AxisNames[i]])
        }
        while ((i+k)<len){
        var yArr = [];
        for (a=0; a<selRows.length; a++){
        yArr.push(selRows[a][AxisNames[(len-k)]])
        }
        //console.log(['xArr', xArr])
        //console.log(['yArr', yArr])
        var trace = {
        x: xArr,
        y: yArr,
        mode: 'markers',
        marker: {
        color: 'orange',
        size: 6
        },
        xaxis: 'x' + (i+1),
        yaxis: 'y' + (i*len+k),
        hoverinfo: 'none'
        };
        Traces.push(trace);
        k++;
        }
        i++;
        k=1;
        }
        Plotly.addTraces(el.id, Traces);
        })}
        ", data = data)
      
    })
    selID <- reactive(input$selID)
    
    pcpDat <- reactive(data[which(data$ID %in% selID()), c(1:(p$nrow+1))])
    output$selectedValues <- renderPrint({pcpDat()$ID})
    colNms <- colnames(data[, c(2:(p$nrow+1))])
    
    boxDat <- data[, c(1:(p$nrow+1))] %>% gather(key, val, -c(ID))
    colnames(boxDat) <- c("ID", "Sample", "Count")
    BP <- ggplot(boxDat, aes(x = Sample, y = Count)) + geom_boxplot()
    ggBP <- ggplotly(BP, width=700)
    
    output$boxPlot <- renderPlotly({
      ggBP %>% onRender("
      function(el, x, data) {
      
      var Traces = [];
      
      var dLength = data.pcpDat.length
      var vLength = data.nVar
      var cNames = data.colNms
      
      for (a=0; a<dLength; a++){
      xArr = [];
      yArr = [];
      for (b=0; b<vLength; b++){
      xArr.push(b+1)
      yArr.push(data.pcpDat[a][cNames[b]]);
      }
      
      var traceHiLine = {
      x: xArr,
      y: yArr,
      mode: 'lines',
      line: {
      color: 'orange',
      width: 1.5
      },
      opacity: 0.9,
      }
      Traces.push(traceHiLine);
      }
      Plotly.addTraces(el.id, Traces);
      
      }", data = list(pcpDat = pcpDat(), nVar = p$nrow, colNms = colNms))})
})
  
  shinyApp(ui, server)
  
    }


scatMatOrthPCP = function(data=data, threshOrth=threshOrth){
  
  key <- val <- ID <- NULL
  
  ui <- shinyUI(fluidPage(
    sliderInput("threshold", "Threshold:", min = 0, max = threshOrth, value=1, step=0.1),
    plotlyOutput("myPlot"),
    plotlyOutput("boxPlot"),
    verbatimTextOutput("selectedValues")
  ))
  
  server <- shinyServer(function(input, output) {
    dat <- data
    
    minVal = min(dat[,-1])
    maxVal = max(dat[,-1])
    # Designate end points of lines to be drawn
    minLine = minVal - 5*(maxVal-minVal)
    maxLine = maxVal + 5*(maxVal-minVal)
    cv = 1
    
    my_fn <- function(data, mapping, ...){
      x = data[,c(as.character(mapping$x))]
      y = data[,c(as.character(mapping$y))]
      p <- ggplot(data = data, aes(x=x, y=y)) + geom_point(alpha=0) + coord_cartesian(xlim = c(minVal, maxVal), ylim = c(minVal, maxVal))
      p
    }
    
    p <- ggpairs(dat[,-1], lower = list(continuous = my_fn))
    
    ggPS <- ggplotly(p)
    
    myLength <- length(ggPS[["x"]][["data"]])
    for (i in 1:myLength){
      item =ggPS[["x"]][["data"]][[i]]$text[1]
      if (!is.null(item))
        if (!startsWith(item, "co")){
          ggPS[["x"]][["data"]][[i]]$hoverinfo <- "none"
        }
    }
 
    output$myPlot <- renderPlotly(ggPS %>%    
#    output$myPlot <- renderPlotly(ggPS %>% config(displayModeBar = F) %>%
onRender("
       function(el, x, data) {
       
       len = Math.sqrt(document.getElementsByClassName('cartesianlayer')[0].childNodes.length);
       AxisNames = [];
       for (i = 1; i < (len+1); i++) {
       AxisNames.push(document.getElementsByClassName('infolayer')[0].childNodes[i].textContent);
       }
       
       noPoint = x.data.length;
       var SubPoints = [];
       var Traces = [];
       var i=0;
       var k=1;
       while ((i*len+k)<=Math.pow((len-1),2)) {
       while ((i+k)<len){
       var selRows = [];
       data.dat.forEach(function(row){
       if(Math.abs(row[AxisNames[i]]-row[AxisNames[(len-k)]]) > Math.sqrt(2)*data.val){
       selRows.push(row);
       }})
       var xArr = [];
       for (a=0; a<selRows.length; a++){
       xArr.push(selRows[a][AxisNames[i]])
       }
       var yArr = [];
       for (a=0; a<selRows.length; a++){
       yArr.push(selRows[a][AxisNames[(len-k)]])
       }
       var keepIndex = [];
       for (a=0; a<selRows.length; a++){
       keepIndex.push(selRows[a]['ID'])
       }
       SubPoints.push(keepIndex);
       
       var tracePoints = {
       x: xArr,
       y: yArr,
       hoverinfo: 'none',
       mode: 'markers',
       marker: {
       color: 'black',
       size: 4
       },
       xaxis: 'x' + (i+1),
       yaxis: 'y' + (i*len+k)
       };
       var traceHiLine = {
       x: [data.minLine, data.maxLine - Math.sqrt(2)*data.val],
       y: [data.minLine + Math.sqrt(2)*data.val, data.maxLine],
       mode: 'lines',
       line: {
       color: 'gray',
       width: 1
       },
       opacity: 0.25,
       xaxis: 'x' + (i+1),
       yaxis: 'y' + (i*len+k)
       }
       var traceLoLine = {
       x: [data.minLine + Math.sqrt(2)*data.val, data.maxLine],
       y: [data.minLine, data.maxLine - Math.sqrt(2)*data.val],
       mode: 'lines',
       fill: 'tonexty',
       line: {
       color: 'gray',
       width: 1
       },
       opacity: 0.25,
       xaxis: 'x' + (i+1),
       yaxis: 'y' + (i*len+k)
       }
       Traces.push(tracePoints);
       Traces.push(traceHiLine);
       Traces.push(traceLoLine);
       k++;
       }
       i++;
       k=1;
       }
       Plotly.addTraces(el.id, Traces);
       
       var idRows = []
       for (a=0; a<data.dat.length; a++){
       idRows.push(data.dat[a]['ID'])
       }
       
       noPoint = x.data.length;
       
       var nseltrace = 0;
       el.on('plotly_selected', function(e) {
       
       if (x.data.length > noPoint){ // had been crossed out//
       Plotly.deleteTraces(el.id, x.data.length-1); // had been crossed out//
       } // had been crossed out//
       numSel = e.points.length
       cN = e.points[0].curveNumber;
       
       var pointNumbers = [];
       for (a=0; a<numSel; a++){
       pointNumbers.push(e.points[a].pointNumber)
       }
       
       // Determine which subplot was selected
       subPlot = (cN - Math.pow(len,2))/3+1
       
       var selDots = []
       var selData = []
       for (a=0; a<pointNumbers.length; a++){
       var selDot = SubPoints[subPlot-1][pointNumbers[a]]
       //selData.push(myDat[selDot])
       selData.push(data.dat[idRows.indexOf(selDot)])
       selDots.push(selDot)
       }
       
       Shiny.onInputChange('selDots', selDots);
       
       
       if (nseltrace>0){
       Plotly.deleteTraces(el.id,range((-1*len*(len-1)/2),-1,1))
       }
       
       var Traces = [];
       var i=0;
       var k=1;
       while ((i*len+k)<=Math.pow((len-1),2)) {
       var xArr = [];
       for (a=0; a<selData.length; a++){
       xArr.push(selData[a][AxisNames[i]])
       }
       while ((i+k)<len){
       var yArr = [];
       for (a=0; a<selData.length; a++){
       yArr.push(selData[a][AxisNames[(len-k)]])
       }
       var trace = {
       x: xArr,
       y: yArr,
       mode: 'markers',
       marker: {
       color: 'red',
       size: 4
       },
       xaxis: 'x' + (i+1),
       yaxis: 'y' + (i*len+k),
       hoverinfo: 'none'
       };
       Traces.push(trace);
       k++;
       }
       i++;
       k=1;
       }
       
       nseltrace = nseltrace+1
       Plotly.addTraces(el.id, Traces);
       })
       
       }
       ", data = list(dat=dat, val = input$threshold, minLine=minLine, maxLine=maxLine)))
    
    selDots <- reactive(input$selDots)
    
    pcpDat <- reactive(dat[which(dat$ID %in% selDots()), ])
    output$selectedValues <- renderPrint({str(pcpDat())})
    colNms <- colnames(dat[, c(2:ncol(dat))])
    
    boxDat <- dat %>% gather(key, val, -c(ID))
    BP <- ggplot(boxDat, aes(x = key, y = val)) + geom_boxplot()
    ggBP <- ggplotly(BP)
    
    output$boxPlot <- renderPlotly({
      ggBP %>% onRender("
                        function(el, x, data) {
                        
                        var Traces = [];
                        
                        var dLength = data.pcpDat.length
                        var vLength = data.nVar
                        var cNames = data.colNms
                        
                        for (a=0; a<dLength; a++){
                        xArr = [];
                        yArr = [];
                        for (b=0; b<vLength; b++){
                        xArr.push(b+1)
                        yArr.push(data.pcpDat[a][cNames[b]]);
                        }
                        
                        var traceHiLine = {
                        x: xArr,
                        y: yArr,
                        mode: 'lines',
                        line: {
                        color: 'orange',
                        width: 1
                        },
                        opacity: 0.9,
                        }
                        Traces.push(traceHiLine);
                        }
                        Plotly.addTraces(el.id, Traces);
                        
                        }", data = list(pcpDat = pcpDat(), nVar = p$nrow, colNms = colNms))})    
})
  
  shinyApp(ui, server)
  
}

scatMatPIPCP = function(data=data){
  
  ID <- key <- val <- NULL
  
  ui <- shinyUI(fluidPage(
    sliderInput("cio", "Prediction interval:", min = 0.8, max = 0.995, value=0.995, step=0.005),
    plotlyOutput("myPlot", height = 620),
    plotlyOutput("boxPlot"),
    verbatimTextOutput("selectedValues")
  ))
  
  server <- shinyServer(function(input, output, session) {
    
    dat <- data
    nCol = ncol(dat)
    conf=seq(0,0.995,.005)
    st<- qt(1-(1-conf)/2,(nrow(dat)-2))
    
    b0 = c()
    b1 = c()
    sse = c()
    for (i in 2:(nCol-1)){
      j = nCol
      while (j >i){
        datXY <- as.data.frame(cbind(x = dat[,i], y = dat[,j]))
        datLm <- lm(y~x,data=datXY)
        b0 <- c(b0, coef(datLm)[1])
        b1 <- c(b1, coef(datLm)[2])
        sse <- c(sse, summary(datLm)[[6]])
        j = j-1
      }
    }
    b0<-as.vector(b0)
    b1<-as.vector(b1)
    sse<-as.vector(sse)
    
    minVal = 0
    maxVal = max(dat[,-1])
    
    my_fn <- function(data, mapping, ...){
      x = data[,c(as.character(mapping$x))]
      y = data[,c(as.character(mapping$y))]
      p <- ggplot(data = data, aes(x=x, y=y)) + geom_point(alpha=0) + coord_cartesian(xlim = c(minVal, maxVal), ylim = c(minVal, maxVal))
      p
    }
    
    p <- ggpairs(dat[,-1], lower = list(continuous = my_fn))
    ggPS <- ggplotly(p)
    
    myLength <- length(ggPS[["x"]][["data"]])
    for (i in 1:myLength){
      item =ggPS[["x"]][["data"]][[i]]$text[1]
      if (!is.null(item))
        if (!startsWith(item, "co")){
          ggPS[["x"]][["data"]][[i]]$hoverinfo <- "none"
        }
    }
    
    cio <- reactive(input$cio)
    
    shiny::observe({
      cio <- cio()
      # Send ci-original level into onRender() function
      session$sendCustomMessage(type = "cio", message=cio)
    })
    
    #output$myPlot <- renderPlotly(ggPS %>% config(displayModeBar = F) %>%
    output$myPlot <- renderPlotly(ggPS %>%
    onRender("
     function(el, x, data) {
     myDat = data.dat
     
     function range(start, stop, step){
     var a=[start], b=start;
     while(b<stop){b+=step;a.push(b)}
     return a;
     }
     
     len = Math.sqrt(document.getElementsByClassName('cartesianlayer')[0].childNodes.length);
     AxisNames = [];
     for (i = 1; i < (len+1); i++) {
     AxisNames.push(document.getElementsByClassName('infolayer')[0].childNodes[i].textContent);
     }
     
     Shiny.addCustomMessageHandler('cio',
     function(ci) {
     
     
     //Plotly.newPlot(el.id, data=[], layout=[]);
     
     stIndex = Math.round((1-0)/.01*ci)
     st = data.st[stIndex]
     
     var Traces = [];
     var i=0;
     var j=0;
     var k=1;
     var SubPoints = [];
     while ((i*len+k)<=Math.pow((len-1),2)) {
     while ((i+k)<len){
     var x = [];
     var y = [];
     var xTotal = 0;
     var ssx = 0;
     n = data.dat.length; // 100
     for (a=0; a<n; a++){
     xa = data.dat[a][AxisNames[i]]
     x.push(xa)
     y.push(data.dat[a][AxisNames[(len-k)]])
     xTotal+=xa
     }
     
     var xm = xTotal/n
     for (a=0; a<n; a++){
     ssx+=Math.pow((data.dat[a][AxisNames[i]] - xm),2)
     }
     
     var minX = -1
     var maxX = 2*Math.max.apply(null,x)
     
     var inc = (maxX-minX)/100
     var xv = [];
     var yv = [];
     var se = [];
     var ci = [];
     var uyv = [];
     var lyv = [];
     var a = minX
     while (a < maxX){
     xv.push(a);
     yva = data.b0[j]+data.b1[j]*a;
     // just changed this to have 1+1/n instead of just 1/n
     sea = data.sse[j] * Math.sqrt(1+1/n+Math.pow((a-xm),2)/ssx);
     yv.push(yva);
     se.push(sea);
     ci.push(st*sea);
     uyv.push(yva+st*sea);
     lyv.push(yva-st*sea);
     a+=inc;
     }
     
     var lwr = [];
     var upr = [];
     var ypred = [];
     var ssea = [];
     var outCI = [];
     var xPoints = [];
     var yPoints = [];
     var keepIndex = []
     for (a=0; a<n; a++){
     xa = data.dat[a][AxisNames[i]]
     // just changed this to have 1+1/n instead of just 1/n
     ssea.push(data.sse[j] * Math.sqrt(1+1/n+Math.pow((xa-xm),2)/ssx))
     ypred.push(data.b0[j]+data.b1[j]*xa)
     lwr.push(ypred[a] - ssea[a]*st)
     upr.push(ypred[a] + ssea[a]*st)
     if (!(y[a]>lwr[a] & y[a]<upr[a])){
     xPoints.push(xa)
     yPoints.push(data.dat[a][AxisNames[(len-k)]])
     keepIndex.push(a)
     }
     }
     SubPoints.push(keepIndex);
     
     var tracePoints = {
     x: xPoints,
     y: yPoints,
     mode: 'markers',
     marker: {
     color: 'black',
     size: 4
     },
     xaxis: 'x' + (i+1),
     yaxis: 'y' + (i*len+k),
     hoverinfo: 'none'
     };
     var hiLine = {
     x: xv,
     y: uyv,
     mode: 'lines',
     line: {
     color: 'gray',
     width: 1
     },
     xaxis: 'x' + (i+1),
     yaxis: 'y' + (i*len+k),
     opacity: 0.25,
     hoverinfo: 'none'
     };
     var lowLine = {
     x: xv,
     y: lyv,
     mode: 'lines',
     fill: 'tonexty',
     line: {
     color: 'gray',
     width: 1
     },
     xaxis: 'x' + (i+1),
     yaxis: 'y' + (i*len+k),
     opacity: 0.25,
     hoverinfo: 'none'
     };
     Traces.push(tracePoints);
     Traces.push(hiLine);
     Traces.push(lowLine);
     j++;
     k++;
     }
     i++;
     k=1;
     }
     Plotly.addTraces(el.id, Traces);
     
     var idRows = []
     for (a=0; a<data.dat.length; a++){
     idRows.push(data.dat[a]['ID'])
     }
     
     var nseltrace = 0;
     el.on('plotly_selected', function(e) {
     
     numSel = e.points.length
     cN = e.points[0].curveNumber;
     
     var pointNumbers = [];
     for (a=0; a<numSel; a++){
     pointNumbers.push(e.points[a].pointNumber)
     }
     
     // Determine which subplot was selected
     subPlot = (cN - Math.pow(len,2))/3+1
     
     var selData = []
     var selDots = []
     for (a=0; a<pointNumbers.length; a++){
     var selDot = SubPoints[subPlot-1][pointNumbers[a]]
     selData.push(myDat[selDot])
     selDots.push(selDot)
     }
     
     Shiny.onInputChange('selDots', selDots);
     
     if (nseltrace>0){
     Plotly.deleteTraces(el.id,range((-1*len*(len-1)/2),-1,1))
     }
     
     var Traces = [];
     var i=0;
     var k=1;
     while ((i*len+k)<=Math.pow((len-1),2)) {
     var xArr = [];
     for (a=0; a<selData.length; a++){
     xArr.push(selData[a][AxisNames[i]])
     }
     while ((i+k)<len){
     var yArr = [];
     for (a=0; a<selData.length; a++){
     yArr.push(selData[a][AxisNames[(len-k)]])
     }
     var trace = {
     x: xArr,
     y: yArr,
     mode: 'markers',
     marker: {
     color: 'red',
     size: 4
     },
     xaxis: 'x' + (i+1),
     yaxis: 'y' + (i*len+k),
     hoverinfo: 'none'
     };
     Traces.push(trace);
     k++;
     }
     i++;
     k=1;
     }
     
     nseltrace = nseltrace+1
     Plotly.addTraces(el.id, Traces);
     })
     })} ", data = list(dat=dat, b0=b0, b1=b1, sse=sse, st=st)))
    
    selID <- reactive(input$selDots)
    
    pcpDat <- reactive(dat[selID()+1, ])
    #output$selectedValues <- renderPrint({str(pcpDat())})
    output$selectedValues <- renderPrint({pcpDat()$ID})
    
    colNms <- colnames(dat[, c(2:(ncol(dat)))])
    
    boxDat <- dat %>% gather(key, val, -c(ID))
    BP <- ggplot(boxDat, aes(x = key, y = val)) + geom_boxplot()
    ggBP <- ggplotly(BP)
    
    output$boxPlot <- renderPlotly({
      ggBP %>% onRender("
        function(el, x, data) {
        
        var Traces = [];
        
        var dLength = data.pcpDat.length
        var vLength = data.nVar
        var cNames = data.colNms
        
        for (a=0; a<dLength; a++){
        xArr = [];
        yArr = [];
        for (b=0; b<vLength; b++){
        xArr.push(b+1)
        yArr.push(data.pcpDat[a][cNames[b]]);
        }
        
        var traceHiLine = {
        x: xArr,
        y: yArr,
        mode: 'lines',
        line: {
        color: 'red',
        width: 1
        },
        opacity: 0.9,
        }
        Traces.push(traceHiLine);
        }
        Plotly.addTraces(el.id, Traces);
        
        }", data = list(pcpDat = pcpDat(), nVar = p$nrow, colNms = colNms))})
  })
  
  shinyApp(ui, server)
    }
