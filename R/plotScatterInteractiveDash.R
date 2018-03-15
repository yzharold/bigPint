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
#' @importFrom ggplot2 ggplot aes_string aes xlim ylim geom_boxplot theme
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
#' plotScatterInteractiveDash(soybean_cn)
plotScatterInteractiveDash = function(data=data, outDir=getwd(), threshOrth=3, xbins=10, option="hexagon"){
    scatMatHexPCP(data=data, xbins=xbins)
} 

##############################################################################
##############################################################################
##############################################################################
# HELPER FUNCTIONS


scatMatHexPCP = function(data=data, xbins=xbins){
  
  counts <- hexID <- key <- val <- ID <- NULL
  
  sidebar <- shinydashboard::dashboardSidebar(
    width = 180,
    shiny::hr(),
    shinydashboard::sidebarMenu(id="tabs",
    shinydashboard::menuItem("Application", tabName="scatMatPlot", selected=TRUE), #hexPlot
    shinydashboard::menuItem("About", tabName = "about") #boxPlot
    )
  )
  
  body <- shinydashboard::dashboardBody(
    shinydashboard::tabItems(
      
    shinydashboard::tabItem(tabName = "scatMatPlot",
    shiny::fluidRow(
      shiny::column(width = 12, shinydashboard::box(width = 660, height = 660, plotly::plotlyOutput("scatMatPlot"), collapsible = FALSE, background = "black", title = "Binned scatterplot", status = "primary", solidHeader = TRUE))),
    
    shiny::fluidRow(
    shiny::column(width = 12,
    shinydashboard::box(width = NULL, plotly::plotlyOutput("boxPlot"), collapsible = FALSE, background = "black", title = "Boxplot", status = "primary", solidHeader = TRUE))),
    
    shiny::fluidRow(
    shiny::column(width = 12,
    shinydashboard::box(width = NULL, shiny::verbatimTextOutput("selectedValues"), collapsible = TRUE, title = "Selected Gene IDs", status = "primary", solidHeader = TRUE)))),
    
    shinydashboard::tabItem(tabName = "about",
    shiny::column(width = 12,
    shiny::fluidRow("This application allows you to examine the relationship between all variables in your dataset with an interactive scatterplot matrix. Plotting points can obscure the number of genes in a given area due to overplotting. As a result, we use hexagon bins in the scatterplot matrix. If you hover over a given hexagon bin of interest, you can determine the number of genes in its area (Figure 1)."),
    br(),
    shiny::fluidRow("You can also click on a given hexagon bin of interest to overlay the genes it contains across all scatterplots as orange points (Figure 2). Doing so will automatically overlay these same genes as orange parallel coordinate lines across a side-by-side boxplot of your data immediately below (Figure 3). Moreover, beneath that, you will see an output of the IDs of theses selection genes (Figure 4)."),
    br(),
    shiny::fluidRow("The four figures below were created in simulated data drawn from the normal distribution for didactic purposes. We hope to improve upon this application by allowing users more customizing options, such as selecting hexagon bin size, changing color mappings, and providing a clear color legend."),
    #shiny::fluidRow(img(src='Figure1.png'))
    #tags$img(src='www/Figure1.png',height='60',width='200')
    shiny::imageOutput("fig1")
    ))
    )
  )
  
  ui <- shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(title = "Overlaying cases", titleWidth = 180),
    sidebar,
    body
  )
  
  server <- function(input, output, session) {
    
    ################################ Prepare scatterplot matrix
    ###########################################################
    
    output$fig1 <- renderImage({
      filename <- normalizePath(file.path("www/Figure1.png"))})
  
    
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
}
  
  shiny::shinyApp(ui = ui, server = server)
  
}

