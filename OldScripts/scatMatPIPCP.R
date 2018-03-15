#' Scatterplot prediction-interval matrix linked to parallel coordinate plot
#' 
#' Scatterplot matrix has contains genes above a given prediction interval. Selecting a subset of genes in any scatterplot leads to the parallel coordinate plot values of its contained genes to be overlaid onto the boxplots below.
#' @param data the data frame
#' @importFrom plotly plotlyOutput ggplotly renderPlotly config
#' @importFrom ggplot2 ggplot aes_string xlim ylim geom_boxplot
#' @importFrom shiny verbatimTextOutput fluidPage reactive renderPrint shinyUI sliderInput shinyServer shinyApp
#' @importFrom htmlwidgets onRender
#' @importFrom tidyr gather
#' @importFrom utils str
#' @importFrom stats qt lm coef
#' @export
#' @examples
#' data(soybean_cn)
#' soybean_cn <- soybean_cn[1:1000,c(1:7)]
#' soybean_cn <- soybean_cn[,c(1:7)]
#' scatMatPIPCP(data = soybean_cn)
#' 
#' data(soybean_ir)
#' soybean_ir[,2:ncol(soybean_ir)] <- log(soybean_ir[,2:ncol(soybean_ir)]+1)
#' scatMatPIPCP(data = soybean_ir)
#' 
#' set.seed(1)
#' f = function(){abs(rnorm(50000))}
#' data = data.frame(ID = paste0("ID", 1:50000), A=f(), B=f(), C=f(), D=f())
#' data$ID <- as.character(data$ID)
#' scatMatPIPCP(data = data)
scatMatPIPCP = function(data){
  
  ID <- key <- val <- NULL
  
  ui <- shinyUI(fluidPage(
    sliderInput("cio", "Prediction interval:", min = 0.8, max = 0.99, value=0.99, step=0.01),
    plotlyOutput("myPlot", height = 620),
    plotlyOutput("boxPlot"),
    verbatimTextOutput("selectedValues")
  ))
  
  server <- shinyServer(function(input, output, session) {
    
    dat <- data
    nCol = ncol(dat)
    conf=seq(0,0.99,.01)
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
    
    p <- ggpairs(dat[1,-1], lower = list(continuous = my_fn))
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
    
    output$myPlot <- renderPlotly(ggPS %>% config(displayModeBar = F) %>%
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
    output$selectedValues <- renderPrint({str(pcpDat())})
    
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