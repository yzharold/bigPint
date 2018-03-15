#' Scatterplot orthogonal-distance matrix linked to parallel coordinate plot
#' 
#' Scatterplot matrix has contains genes above a given orthogonal distance. Selecting a subset of genes in any scatterplot leads to the parallel coordinate plot values of its contained genes to be overlaid onto the boxplots below.
#' @param data the data frame
#' @param maxOrthThresh maximum orthogonal threshold value (default 20)
#' @importFrom plotly plotlyOutput ggplotly renderPlotly config
#' @importFrom ggplot2 ggplot aes_string xlim ylim geom_boxplot
#' @importFrom shiny verbatimTextOutput fluidPage reactive renderPrint shinyUI sliderInput shinyServer shinyApp
#' @importFrom htmlwidgets onRender
#' @importFrom tidyr gather
#' @importFrom utils str
#' @export
#' @examples
#' set.seed(1)
#' f = function(){abs(rnorm(50000))}
#' data = data.frame(ID = paste0("ID", 1:50000), A=f(), B=f(), C=f(), D=f())
#' data$ID <- as.character(data$ID)
#' scatMatOrthPCP(data = data)
#' 
#' data(soybean_cn)
#' data = soybean_cn
#' scatMatOrthPCP(data = data[,c(1:7)])
scatMatOrthPCP = function(data, maxOrthThresh=20){

  key <- val <- ID <- NULL
  
ui <- shinyUI(fluidPage(
  sliderInput("threshold", "Threshold:", min = 0, max = maxOrthThresh, value=1, step=0.1),
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

  output$myPlot <- renderPlotly(ggPS %>% config(displayModeBar = F) %>%
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

//if (x.data.length > noPoint){
//Plotly.deleteTraces(el.id, x.data.length-1);
//}
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