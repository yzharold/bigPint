#' Scatterplot hexagon matrix linked to parallel coordinate plot
#' 
#' Scatterplot matrix has contains hexagons that bin the number of genes. Clicking on a hexagon of interest leads to the parallel coordinate plot values of its contained genes to be overlaid onto the boxplots below.
#' @param data the data frame
#' @param xbins the number of bins partitioning the range of the plot (default 16)
#' @importFrom plotly plotlyOutput ggplotly renderPlotly
#' @importFrom ggplot2 ggplot aes xlim ylim geom_boxplot
#' @importFrom shiny verbatimTextOutput fluidPage reactive renderPrint shinyUI HTML shinyServer shinyApp
#' @importFrom htmlwidgets onRender
#' @importFrom hexbin hexbin hcell2xy
#' @importFrom tidyr gather
#' @importFrom utils str
#' @export
#' @examples
#' data(soybean_cn)
#' scatMatHexPCP(data = soybean_cn[,c(1:4,8:10)])
scatMatHexPCP = function(data, xbins=16){

  counts <- hexID <- key <- val <- ID <- NULL
  
ui <- shinyUI(fluidPage(
  plotlyOutput("scatMatPlot", width = 700, height = 700),
  plotlyOutput("boxPlot", width = 700),
  HTML("<br><br><br>"),
  verbatimTextOutput("selectedValues")
))

server <- shinyServer(function(input, output) {
  
  ################################ Prepare scatterplot matrix
  ###########################################################
  
  maxVal = max(abs(data[,-1]))
  maxRange = c(-1*maxVal, maxVal)
  #maxRange = c(0, maxVal)
  
  my_fn <- function(data, mapping, ...){
    x = data[,c(as.character(mapping$x))]
    y = data[,c(as.character(mapping$y))]
    h <- hexbin(x=x, y=y, xbins=xbins, shape=1, IDs=TRUE, xbnds=maxRange, ybnds=maxRange)
    hexdf <- data.frame (hcell2xy (h),  hexID = h@cell, counts = h@count)
    attr(hexdf, "cID") <- h@cID
    p <- ggplot(hexdf, aes(x=x, y=y, fill = counts, hexID=hexID)) + geom_hex(stat="identity") + geom_abline(intercept = 0, color = "red", size = 0.25) + coord_cartesian(xlim = c(-0.5, maxRange[2]+0.5), ylim = c(-0.5, maxRange[2]+0.5))
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
  output$selectedValues <- renderPrint({str(pcpDat())})
  colNms <- colnames(data[, c(2:(p$nrow+1))])

  boxDat <- data[, c(1:(p$nrow+1))] %>% gather(key, val, -c(ID))
  BP <- ggplot(boxDat, aes(x = key, y = val)) + geom_boxplot()
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
