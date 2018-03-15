#' Scatterplot fold change matrix linked to parallel coordinate plot
#' 
#' Scatterplot matrix has contains genes above a given fold change. Selecting a subset of genes in any scatterplot leads to the parallel coordinate plot values of its contained genes to be overlaid onto the boxplots below.
#' @param data the data frame
#' @importFrom plotly plotlyOutput ggplotly renderPlotly
#' @importFrom ggplot2 ggplot aes_string xlim ylim geom_boxplot
#' @importFrom shiny verbatimTextOutput fluidPage reactive renderPrint shinyUI sliderInput shinyServer shinyApp
#' @importFrom htmlwidgets onRender
#' @importFrom utils str
#' @importFrom tidyr gather
#' @export
#' @examples
#' set.seed(1)
#' f = function(){abs(rnorm(50000))}
#' data = data.frame(ID = paste0("ID", 1:50000), A=f(), B=f(), C=f(), D=f())
#' data$ID <- as.character(data$ID)
#' scatMatFCPCP(data = data)
scatMatFCPCP = function(data){

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

  my_fn <- function(data, mapping, ...){
    x = data[,c(as.character(mapping$x))]
    y = data[,c(as.character(mapping$y))]
    p <- ggplot(data = data, aes(x=x, y=y)) + geom_point(alpha=0) + coord_cartesian(xlim = c(minVal, maxVal), ylim = c(minVal, maxVal))
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