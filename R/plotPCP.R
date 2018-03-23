#' Highlight parallel coordinate plot lines inside selection box
#' 
#' @param pcpDat the data frame that contains the parallel coordinate plot values
#' @param option the interactivity option ("deleteInteger", "delete", "highlight"); default ("deleteInteger")
#' @importFrom plotly plotlyOutput ggplotly renderPlotly layout
#' @importFrom ggplot2 ggplot aes_string geom_point xlim ylim scale_x_discrete
#' @importFrom shiny verbatimTextOutput fluidPage reactive renderPrint shinyApp bootstrapPage
#' @importFrom htmlwidgets onRender
#' @importFrom utils str
#' @export
#' @examples
#' set.seed(3)
#' f = function(){1.3*rnorm(500)}
#' pcpDat = data.frame(ID = paste0("ID", 1:500), A.1=f(), A.2=f(), A.3=f(), B.1=f(), B.2=f(), B.3=f())
#' pcpDat$ID = as.character(pcpDat$ID)
#' plotPCP(pcpDat = pcpDat)
plotPCP = function(pcpDat, option = "deleteInteger"){
  
  if (option=="delete"){
    selDelPCP(pcpDat)
  }
  else if (option=="highlight"){
    selPCP(pcpDat)
  }  
  else{
    selDelIntShadePCP(pcpDat)  
  }
} 

##############################################################################
##############################################################################
##############################################################################
# HELPER FUNCTIONS

#' Delete parallel coordinate plot lines inside integers of selection box
#' 
#' @param pcpDat the data frame that contains the parallel coordinate plot values
selDelIntPCP = function(pcpDat){
  
  ui <- basicPage(
    plotlyOutput("plot1")
  )
  
  server <- function(input, output) {
    
    colNms <- colnames(pcpDat[, c(2:(ncol(pcpDat)))])
    nVar <- length(colNms)
    browser()
    
    p <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point(alpha=0) + xlim(0,(nVar-1)) +ylim(min(pcpDat[,2:(nVar+1)]),max(pcpDat[,2:(nVar+1)])) + xlab("Sample") + ylab("Count")
    gp <- ggplotly(p)
    
    output$plot1 <- renderPlotly({
      gp %>% onRender("
      function(el, x, data) {
      
      var origPcpDat = data.pcpDat
      var pcpDat = data.pcpDat
      
      var Traces = [];
      var dLength = pcpDat.length
      var vLength = data.nVar
      var cNames = data.colNms

      xArr = [];
      for (b=0; b<vLength; b++){
        xArr.push(b)
      }

      for (a=0; a<dLength; a++){
        yArr = [];
        for (b=0; b<vLength; b++){
          yArr.push(pcpDat[a][cNames[b]]);
        }
        var pcpLine = {
          x: xArr,
          y: yArr,
          mode: 'lines',
          line: {
            color: 'red',
            width: 1
          },
          opacity: 0.9,
        }
        Traces.push(pcpLine);
      }
      Plotly.addTraces(el.id, Traces);
      
      var selectTime = 1
      
      el.on('plotly_selected', function(e) {
      if (pcpDat.length>0){
      
      var dLength = pcpDat.length
      var selectedPCP = []
      var xMin = e.range.x[0]
      var xMax = e.range.x[1]
      var xMinC = Math.abs(Math.ceil(xMin))
      var xMaxF = Math.floor(xMax)
      var yMin = e.range.y[0]
      var yMax = e.range.y[1]
      var integers = []
      
      if (!((xMax<0) || (xMin>(vLength-1)))){
      for (a=xMinC; a<(xMaxF+1); a++){
      integers.push(a)
      }
      }
      var iLength = integers.length
      //console.log(['iLength', iLength])
      
      var selectedPCP = []
      var notSelectedPCP = []
      
      if (selectTime==1){
      if (iLength > 0){
      for (a=0; a<dLength; a++){
      var dat = pcpDat[a]
      var isOut = 0;
      for (b=0; b<iLength; b++){
      var yVal = dat[cNames[integers[b]]]
      if (!(yMin < yVal && yVal < yMax)){
      isOut = 1;
      }
      }
      if (isOut==1){
      selectedPCP.push(a)
      }
      }
      //console.log(['selectedPCP', selectedPCP])
      
      var updateSPCP = []
      var selectedPCPL = selectedPCP.length
      for (a=0; a<selectedPCPL; a++){
      updateSPCP[a]=selectedPCP[a]+1
      }
      
      var update = {
      line: {
      color: 'blue',
      width: 1
      }
      }
      if (selectedPCPL!=0){
      Plotly.deleteTraces(el.id, updateSPCP);
      }
      
      var newDat = []
      var selectedPCPL = selectedPCP.length
      for (a=0; a<dLength; a++){
      var equal = 0;
      for (b=0; b<selectedPCPL; b++){
      if (a==selectedPCP[b]){
      equal=1
      }
      }
      if (equal==0){
      newDat.push(pcpDat[a])
      }
      }
      pcpDat = newDat
      }
      }
      
      if (selectTime>1){
      if (iLength > 0){
      for (a=0; a<dLength; a++){
      var dat = pcpDat[a]
      var isOut = 0;
      for (b=0; b<iLength; b++){
      var yVal = dat[cNames[integers[b]]]
      if (!(yMin < yVal && yVal < yMax)){
      isOut = 1;
      }
      }
      
      if (isOut==0){
      selectedPCP.push(a)
      }
      else{
      notSelectedPCP.push(a)
      }
      }
      //console.log(['notSelectedPCP', notSelectedPCP])
      
      var updateNSPCP = []
      var notSelectedPCPL = notSelectedPCP.length
      for (a=0; a<notSelectedPCPL; a++){
      updateNSPCP[a]=notSelectedPCP[a]+1
      }
      
      if (notSelectedPCPL!=0){
      //console.log(['deleting'], updateNSPCP)
      Plotly.deleteTraces(el.id, updateNSPCP);
      }
      
      var newDat = []
      var selectedPCPL = selectedPCP.length
      for (a=0; a<dLength; a++){
      var equal = 0;
      for (b=0; b<selectedPCPL; b++){
      if (a==selectedPCP[b]){
      equal=1
      }
      }
      if (equal==1){
      newDat.push(pcpDat[a])
      }
      }
      pcpDat = newDat
      }
      
      else{
      for (a=0; a<dLength; a++){
      notSelectedPCP.push(a)
      }
      
      var updateNSPCP = []
      var notSelectedPCPL = notSelectedPCP.length
      for (a=0; a<notSelectedPCPL; a++){
      updateNSPCP[a]=notSelectedPCP[a]+1
      }
      
      var update = {
      line: {
      color: 'red',
      width: 1
      }
      }
      if (notSelectedPCPL!=0){
      Plotly.deleteTraces(el.id, update, updateNSPCP);
      }
      pcpDat = []
      }
      }
      
      var sel1 = {
      x: [xMin, xMax],
      y: [yMin, yMin],
      mode: 'lines',
      //fill: 'tonexty',
      line: {
      color: 'black',
      width: 0.5,
      dash: 'dot'
      },
      hoverinfo: 'none',
      }
      var sel2 = {
      x: [xMax, xMax],
      y: [yMin, yMax],
      mode: 'lines',
      //fill: 'tonexty',
      line: {
      color: 'black',
      width: 0.5,
      dash: 'dot'
      },
      hoverinfo: 'none'
      }
      var sel3 = {
      x: [xMin, xMax],
      y: [yMax, yMax],
      mode: 'lines',
      //fill: 'tonexty',
      line: {
      color: 'black',
      dash: 'dot',
      width: 0.5
      },
      hoverinfo: 'none'
      }
      var sel4 = {
      x: [xMin, xMin],
      y: [yMin, yMax],
      mode: 'lines',
      //fill: 'tonexty',
      line: {
      color: 'black',
      dash: 'dot',
      width: 0.5
      },
      hoverinfo: 'none'
      }
      Plotly.addTraces(el.id, [sel1, sel2, sel3, sel4]);
      
      //console.log(['pcpDat', pcpDat])
      
      selectTime++
      }
      })
      }", data = list(pcpDat = pcpDat, nVar = nVar, colNms = colNms))})
    }
  shinyApp(ui, server)
  }

#' Delete parallel coordinate plot lines inside integers of shaded selection box
#' 
#' @param pcpDat the data frame that contains the parallel coordinate plot values
selDelIntShadePCP = function(pcpDat){
  
  ui <- basicPage(
    plotlyOutput("plot1"),
    verbatimTextOutput("rectdf")
  )
  
  server <- function(input, output) {
    
    colNms <- colnames(pcpDat[, c(2:(ncol(pcpDat)))])
    nVar <- length(colNms)

    p <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point(alpha=0) + xlim(0,(nVar-1)) +ylim(min(pcpDat[,2:(nVar+1)]),max(pcpDat[,2:(nVar+1)])) + xlab("Sample") + ylab("Count") + scale_x_discrete(limits=colnames(pcpDat[-1]))
    gp <- ggplotly(p)
    
    inputRectDf <- reactive({
      req(input$rects)
      # data comes back as a big character vector
      # so we reformat it as a dataframe here
      df <- data.frame(t(matrix(input$rects,nrow=8)))
      names(df) <- names(input$rects)[1:8]
      return(df)
    })
    output$rectdf <- renderPrint({print(inputRectDf())})
    
    output$plot1 <- renderPlotly({
      gp %>% onRender("
        function(el, x, data) {
        var rects = [];
        var origPcpDat = data.pcpDat
        var pcpDat = data.pcpDat
        
        var Traces = [];
        var dLength = pcpDat.length
        var vLength = data.nVar
        var cNames = data.colNms
        console.log(['cNames', cNames])
        for (a=0; a<dLength; a++){
        xArr = [];
        yArr = [];
        for (b=0; b<vLength; b++){
        // used to be b
        xArr.push(b+1)
        yArr.push(pcpDat[a][cNames[b]]);
        }
        var pcpLine = {
        x: xArr,
        y: yArr,
        mode: 'lines',
        line: {
        color: 'red',
        width: 1
        },
        opacity: 0.9,
        }
        Traces.push(pcpLine);
        }
        Plotly.addTraces(el.id, Traces);
        
        var selectTime = 1
        
        el.on('plotly_selected', function(e) {
        if (pcpDat.length>0){
        
        var dLength = pcpDat.length
        var selectedPCP = []
        var xMin = e.range.x[0]
        var xMax = e.range.x[1]
        var xMinC = Math.abs(Math.ceil(xMin))
        var xMaxF = Math.floor(xMax)
        var yMin = e.range.y[0]
        var yMax = e.range.y[1]
        var integers = []
        
        // used to be vLength-1
        if (!((xMax<0) || (xMin>(vLength)))){
        for (a=xMinC; a<(xMaxF+1); a++){
        // used to be a
        integers.push(a-1)
        }
        }
        var iLength = integers.length
        //console.log(['iLength', iLength])
        
        var selectedPCP = []
        var notSelectedPCP = []
        
        if (selectTime==1){
        if (iLength > 0){
        for (a=0; a<dLength; a++){
        var dat = pcpDat[a]
        var isOut = 0;
        for (b=0; b<iLength; b++){
        var yVal = dat[cNames[integers[b]]]
        if (!(yMin < yVal && yVal < yMax)){
        isOut = 1;
        }
        }
        if (isOut==1){
        selectedPCP.push(a)
        }
        }
        //console.log(['selectedPCP', selectedPCP])
        
        var updateSPCP = []
        var selectedPCPL = selectedPCP.length
        for (a=0; a<selectedPCPL; a++){
        updateSPCP[a]=selectedPCP[a]+1
        }
        
        var update = {
        line: {
        color: 'blue',
        width: 1
        }
        }
        if (selectedPCPL!=0){
        Plotly.deleteTraces(el.id, updateSPCP);
        }
        
        var newDat = []
        var selectedPCPL = selectedPCP.length
        for (a=0; a<dLength; a++){
        var equal = 0;
        for (b=0; b<selectedPCPL; b++){
        if (a==selectedPCP[b]){
        equal=1
        }
        }
        if (equal==0){
        newDat.push(pcpDat[a])
        }
        }
        pcpDat = newDat
        }
        }
        
        if (selectTime>1){
        if (iLength > 0){
        for (a=0; a<dLength; a++){
        var dat = pcpDat[a]
        var isOut = 0;
        for (b=0; b<iLength; b++){
        var yVal = dat[cNames[integers[b]]]
        if (!(yMin < yVal && yVal < yMax)){
        isOut = 1;
        }
        }
        
        if (isOut==0){
        selectedPCP.push(a)
        }
        else{
        notSelectedPCP.push(a)
        }
        }
        //console.log(['notSelectedPCP', notSelectedPCP])
        
        var updateNSPCP = []
        var notSelectedPCPL = notSelectedPCP.length
        for (a=0; a<notSelectedPCPL; a++){
        updateNSPCP[a]=notSelectedPCP[a]+1
        }
        //console.log(['updateNSPCP'], updateNSPCP)
        
        if (notSelectedPCPL!=0){
        //console.log(['deleting'], updateNSPCP)
        Plotly.deleteTraces(el.id, updateNSPCP);
        }
        
        var newDat = []
        var selectedPCPL = selectedPCP.length
        for (a=0; a<dLength; a++){
        var equal = 0;
        for (b=0; b<selectedPCPL; b++){
        if (a==selectedPCP[b]){
        equal=1
        }
        }
        if (equal==1){
        newDat.push(pcpDat[a])
        }
        }
        pcpDat = newDat
        }
        
        else{
        for (a=0; a<dLength; a++){
        notSelectedPCP.push(a)
        }
        
        var updateNSPCP = []
        var notSelectedPCPL = notSelectedPCP.length
        for (a=0; a<notSelectedPCPL; a++){
        updateNSPCP[a]=notSelectedPCP[a]+1
        }
        
        var update = {
        line: {
        color: 'red',
        width: 1
        }
        }
        if (notSelectedPCPL!=0){
        Plotly.deleteTraces(el.id, update, updateNSPCP);
        }
        pcpDat = []
        }
        }
        
        
        var drawRect = {
        type: 'rect',
        x0: xMin,
        y0: yMin,
        x1: xMax,
        y1: yMax,
        line: {
        color: 'gray',
        width: 1
        },
        fillcolor: 'gray',
        opacity: 0.25
        }
        rects.push(drawRect);
        var update = {
        shapes:rects
        }
        Plotly.relayout(el.id, update)
        Shiny.onInputChange('rects', rects); // make the rects available to shiny
        
        selectTime++
        }
        })
        }", data = list(pcpDat = pcpDat, nVar = nVar, colNms = colNms))})
    }
  shinyApp(ui, server)
  }

#' Delete parallel coordinate plot lines inside selection box
#' 
#' @param pcpDat the data frame that contains the parallel coordinate plot values
selDelPCP = function(pcpDat){
  
  ui <- fluidPage(
    plotlyOutput("plot1", height = 650)
  )
  
  server <- function(input, output) {
    colNms <- colnames(pcpDat[, c(2:(ncol(pcpDat)))])
    nVar <- length(colNms)
    
    p <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point(alpha=0) + xlim(0,(nVar-1)) +ylim(min(pcpDat[,2:(nVar+1)]),max(pcpDat[,2:(nVar+1)])) + xlab("Sample") + ylab("Count")
    gp <- ggplotly(p)
    
    output$plot1 <- renderPlotly({
      
      gp %>% onRender("
      function(el, x, data) {
      var pcpDat = data.pcpDat
      function range(start, stop, step){
      var a=[start], b=start;
      while(b<stop){b+=step;a.push(b)}
      return a;
      };
      var Traces = [];
      var dLength = pcpDat.length
      var vLength = data.nVar
      var cNames = data.colNms
      for (a=0; a<dLength; a++){
      xArr = [];
      yArr = [];
      for (b=0; b<vLength; b++){
      xArr.push(b)
      yArr.push(pcpDat[a][cNames[b]]);
      }
      var pcpLine = {
      x: xArr,
      y: yArr,
      mode: 'lines',
      line: {
      color: 'red',
      width: 1.5
      },
      opacity: 0.9,
      }
      Traces.push(pcpLine);
      }
      Plotly.addTraces(el.id, Traces);
      el.on('plotly_selected', function(e) {
      //console.log(pcpDat)
      var dLength = pcpDat.length
      var selectedPCP = []
      var xMin = e.range.x[0]
      var xMax = e.range.x[1]
      var xMinF = Math.floor(xMin)
      var xMinC = Math.ceil(xMin)
      var xMaxF = Math.floor(xMax)
      var xMaxC = Math.ceil(xMax)
      var yMin = e.range.y[0]
      var yMax = e.range.y[1]
      if (xMin<0){
      xMin = 0
      xMinF = 0
      xMinC = 1
      }
      if (xMax>(vLength-1)){
      xMax = vLength-1
      xMaxF = vLength -2
      xMaxC = vLength-1
      }
      if (!((xMax<0) || (xMin>(vLength-1)))){
      for (a=0; a<dLength; a++){
      var dat = pcpDat[a]
      //console.log(pcpDat[a])
      var yAtXminF = dat[cNames[xMinF]]
      var yAtXminC = dat[cNames[xMinF+1]]
      var yAtXmaxF = dat[cNames[xMaxF+1]]
      var yAtXmaxC = dat[cNames[xMaxF]]
      leftInt = xMinF
      while (leftInt < xMaxC){
      var rightInt = leftInt+1
      var yAtXmin = (xMin-xMinF)*(yAtXminC-yAtXminF)/(xMinC-xMinF) + yAtXminF
      var yAtXmax = (xMax-xMaxF)*(yAtXmaxF-yAtXmaxC)/(xMaxC-xMaxF) + yAtXmaxC
      if (leftInt == xMinF && yMin < yAtXmin && yAtXmin < yMax){
      selectedPCP.push(a)
      leftInt = xMaxC-1
      }
      else if (xMinF == xMaxF){
      if (Math.sign(yMin-yAtXmin)!=Math.sign(yMin-yAtXmax)){
      selectedPCP.push(a)
      leftInt = xMaxC-1
      }
      else if (Math.sign(yMax-yAtXmin)!=Math.sign(yMax-yAtXmax)){
      selectedPCP.push(a)
      leftInt = xMaxC-1
      }
      }
      else if (leftInt == xMinF && xMinF != xMaxF){
      var yLeftInt = dat[cNames[leftInt]]
      var yRightInt = dat[cNames[rightInt]]
      if ((Math.sign(yMin-yLeftInt)!=Math.sign(yMin-yRightInt) ||
      Math.sign(yMax-yLeftInt)!=Math.sign(yMax-yRightInt)) &&
      (Math.sign(yMin-yAtXmin) == Math.sign(yMin-yLeftInt) &&
      Math.sign(yMax-yAtXmin) == Math.sign(yMax-yLeftInt))){
      selectedPCP.push(a)
      leftInt = xMaxC-1
      }
      }
      else if (leftInt == xMaxF && xMinF != xMaxF){
      var yLeftInt = dat[cNames[leftInt]]
      var yRightInt = dat[cNames[rightInt]]
      if ((Math.sign(yMin-yLeftInt)!=Math.sign(yMin-yRightInt) ||
      Math.sign(yMax-yLeftInt)!=Math.sign(yMax-yRightInt)) &&
      (Math.sign(yMin-yAtXmax) == Math.sign(yMin-yRightInt) &&
      Math.sign(yMax-yAtXmax) == Math.sign(yMax-yRightInt))){
      selectedPCP.push(a)
      leftInt = xMaxC-1
      }
      }
      else if (leftInt >= xMin && rightInt <= xMax){
      var yLeftInt = dat[cNames[leftInt]]
      var yRightInt = dat[cNames[rightInt]]
      if (Math.sign(yMin-yLeftInt)!=Math.sign(yMin-yRightInt) ||
      Math.sign(yMax-yLeftInt)!=Math.sign(yMax-yRightInt)){
      selectedPCP.push(a)
      leftInt = xMaxC-1
      }
      }
      leftInt++;
      }
      }
      }
      var updateSPCP = []
      var selectedPCPL = selectedPCP.length
      for (a=0; a<selectedPCPL; a++){
      updateSPCP[a]=selectedPCP[a]+1
      }
      //console.log(selectedPCP)
      Plotly.deleteTraces(el.id, updateSPCP);
      var newDat = []
      var selectedPCPL = selectedPCP.length
      for (a=0; a<dLength; a++){
      var equal = 0;
      for (b=0; b<selectedPCPL; b++){
      if (a==selectedPCP[b]){
      equal=1
      }
      }
      if (equal==0){
      newDat.push(pcpDat[a])
      }
      }
      pcpDat = newDat
      })
      }", data = list(pcpDat = pcpDat, nVar = nVar, colNms = colNms))})
  }
  shinyApp(ui, server)
    }

#' Delete parallel coordinate plot lines inside shaded selection box
#' 
#' @param pcpDat the data frame that contains the parallel coordinate plot values
selDelShadePCP = function(pcpDat){
  
  ui <- basicPage(
    plotlyOutput("plot1", height = 650)
  )
  
  server <- function(input, output) {
    colNms <- colnames(pcpDat[, c(2:(ncol(pcpDat)))])
    nVar <- length(colNms)    
    p <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point(alpha=0) + xlim(0,(nVar-1)) +ylim(min(pcpDat[,2:(nVar+1)]),max(pcpDat[,2:(nVar+1)])) + xlab("Sample") + ylab("Count")
    gp <- ggplotly(p)
    
    output$plot1 <- renderPlotly({
      gp %>% layout(dragmode="select") %>% onRender("
      function(el, x, data) {
      
      var origPcpDat = data.pcpDat
      var pcpDat = data.pcpDat
      
      var Traces = [];
      var dLength = pcpDat.length
      var vLength = data.nVar
      var cNames = data.colNms
      for (a=0; a<dLength; a++){
      xArr = [];
      yArr = [];
      for (b=0; b<vLength; b++){
      xArr.push(b+1)
      yArr.push(pcpDat[a][cNames[b]]);
      }
      var pcpLine = {
      x: xArr,
      y: yArr,
      mode: 'lines',
      line: {
      color: 'red',
      width: 1.5
      },
      opacity: 0.9,
      }
      Traces.push(pcpLine);
      }
      Plotly.addTraces(el.id, Traces);
      
      var selectTime = 1
      
      el.on('plotly_selected', function(e) {
      if (pcpDat.length>0){
      
      var dLength = pcpDat.length
      var selectedPCP = []
      var xMin = e.range.x[0]
      var xMax = e.range.x[1]
      var xMinC = Math.abs(Math.ceil(xMin))
      var xMaxF = Math.floor(xMax)
      var yMin = e.range.y[0]
      var yMax = e.range.y[1]
      var integers = []
      
      if (!((xMax<0) || (xMin>(vLength-1)))){
      for (a=xMinC; a<(xMaxF+1); a++){
      integers.push(a)
      }
      }
      var iLength = integers.length
      //console.log(['iLength', iLength])
      
      var selectedPCP = []
      var notSelectedPCP = []
      
      if (selectTime==1){
      if (iLength > 0){
      for (a=0; a<dLength; a++){
      var dat = pcpDat[a]
      var isOut = 0;
      for (b=0; b<iLength; b++){
      var yVal = dat[cNames[integers[b]]]
      if (!(yMin < yVal && yVal < yMax)){
      isOut = 1;
      }
      }
      if (isOut==1){
      selectedPCP.push(a)
      }
      }
      //console.log(['selectedPCP', selectedPCP])
      
      var updateSPCP = []
      var selectedPCPL = selectedPCP.length
      for (a=0; a<selectedPCPL; a++){
      updateSPCP[a]=selectedPCP[a]+1
      }
      
      var update = {
      line: {
      color: 'blue',
      width: 1
      }
      }
      if (selectedPCPL!=0){
      Plotly.deleteTraces(el.id, updateSPCP);
      }
      
      var newDat = []
      var selectedPCPL = selectedPCP.length
      for (a=0; a<dLength; a++){
      var equal = 0;
      for (b=0; b<selectedPCPL; b++){
      if (a==selectedPCP[b]){
      equal=1
      }
      }
      if (equal==0){
      newDat.push(pcpDat[a])
      }
      }
      pcpDat = newDat
      }
      }
      
      if (selectTime>1){
      if (iLength > 0){
      for (a=0; a<dLength; a++){
      var dat = pcpDat[a]
      var isOut = 0;
      for (b=0; b<iLength; b++){
      var yVal = dat[cNames[integers[b]]]
      if (!(yMin < yVal && yVal < yMax)){
      isOut = 1;
      }
      }
      
      if (isOut==0){
      selectedPCP.push(a)
      }
      else{
      notSelectedPCP.push(a)
      }
      }
      //console.log(['notSelectedPCP', notSelectedPCP])
      
      var updateNSPCP = []
      var notSelectedPCPL = notSelectedPCP.length
      for (a=0; a<notSelectedPCPL; a++){
      updateNSPCP[a]=notSelectedPCP[a]+1
      }
      //console.log(['updateNSPCP'], updateNSPCP)
      
      if (notSelectedPCPL!=0){
      //console.log(['deleting'], updateNSPCP)
      Plotly.deleteTraces(el.id, updateNSPCP);
      }
      
      var newDat = []
      var selectedPCPL = selectedPCP.length
      for (a=0; a<dLength; a++){
      var equal = 0;
      for (b=0; b<selectedPCPL; b++){
      if (a==selectedPCP[b]){
      equal=1
      }
      }
      if (equal==1){
      newDat.push(pcpDat[a])
      }
      }
      pcpDat = newDat
      }
      
      else{
      for (a=0; a<dLength; a++){
      notSelectedPCP.push(a)
      }
      
      var updateNSPCP = []
      var notSelectedPCPL = notSelectedPCP.length
      for (a=0; a<notSelectedPCPL; a++){
      updateNSPCP[a]=notSelectedPCP[a]+1
      }
      
      var update = {
      line: {
      color: 'red',
      width: 1
      }
      }
      if (notSelectedPCPL!=0){
      Plotly.deleteTraces(el.id, update, updateNSPCP);
      }
      pcpDat = []
      }
      }
      
      var Traces = []
      var drawRect = {
      type: 'rect',
      x0: xMin,
      y0: yMin,
      x1: xMax,
      y1: yMax,
      line: {
      color: 'gray',
      width: 1
      },
      fillcolor: 'gray',
      opacity: 0.3
      }
      var update = {
      shapes:[drawRect],
      persistent: 'true'
      }
      Plotly.relayout(el.id, update)
      
      selectTime++
      }
      })
      }", data = list(pcpDat = pcpDat, nVar = nVar, colNms = colNms))})

    }
  shinyApp(ui, server)
  }

#' Highlight parallel coordinate plot lines inside selection box
#' 
#' @param pcpDat the data frame that contains the parallel coordinate plot values
selPCP = function(pcpDat){
  
  ui <- fluidPage(
    plotlyOutput("plot1", height = 500),
    verbatimTextOutput("selectedPCP")
  )
  
  server <- function(input, output) {
    colNms <- colnames(pcpDat[, c(2:(ncol(pcpDat)))])
    nVar <- length(colNms)
    
    p <- ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point(alpha=0) + xlim(0,(nVar-1)) +ylim(min(pcpDat[,2:(nVar+1)]),max(pcpDat[,2:(nVar+1)])) + xlab("Sample") + ylab("Count") + scale_x_discrete(limits=colnames(pcpDat[-1]))
    gp <- ggplotly(p)
    
    output$plot1 <- renderPlotly({
      gp %>% onRender("
      function(el, x, data) {
      var pcpDat = data.pcpDat
      function range(start, stop, step){
      var a=[start], b=start;
      while(b<stop){b+=step;a.push(b)}
      return a;
      };
      var Traces = [];
      var dLength = pcpDat.length
      var vLength = data.nVar
      var cNames = data.colNms
      for (a=0; a<dLength; a++){
      xArr = [];
      yArr = [];
      for (b=0; b<vLength; b++){
      xArr.push(b+1)
      yArr.push(pcpDat[a][cNames[b]]);
      }
      var pcpLine = {
      x: xArr,
      y: yArr,
      mode: 'lines',
      line: {
      color: 'red',
      width: 1.5
      },
      }
      Traces.push(pcpLine);
      }
      Plotly.addTraces(el.id, Traces);
      el.on('plotly_selected', function(e) {
      //console.log(pcpDat)
      var dLength = pcpDat.length
      var selectedPCP = []
      var xMin = e.range.x[0]-1
      var xMax = e.range.x[1]-1
      var xMinF = Math.floor(xMin)
      var xMinC = Math.ceil(xMin)
      var xMaxF = Math.floor(xMax)
      var xMaxC = Math.ceil(xMax)
      var yMin = e.range.y[0]
      var yMax = e.range.y[1]
      if (xMin<0){
      xMin = 0
      xMinF = 0
      xMinC = 1
      }
      if (xMax>(vLength-1)){
      xMax = vLength-1
      xMaxF = vLength -2
      xMaxC = vLength-1
      }
      if (!((xMax<0) || (xMin>(vLength)))){
      for (a=0; a<dLength; a++){
      var dat = pcpDat[a]
      //console.log(pcpDat[a])
      var yAtXminF = dat[cNames[xMinF]]
      var yAtXminC = dat[cNames[xMinF+1]]
      var yAtXmaxF = dat[cNames[xMaxF+1]]
      var yAtXmaxC = dat[cNames[xMaxF]]
      leftInt = xMinF
      while (leftInt < xMaxC){
      var rightInt = leftInt+1
      var yAtXmin = (xMin-xMinF)*(yAtXminC-yAtXminF)/(xMinC-xMinF) + yAtXminF
      var yAtXmax = (xMax-xMaxF)*(yAtXmaxF-yAtXmaxC)/(xMaxC-xMaxF) + yAtXmaxC
      if (leftInt == xMinF && yMin < yAtXmin && yAtXmin < yMax){
      selectedPCP.push(a)
      leftInt = xMaxC
      }
      else if (xMinF == xMaxF){
      if (Math.sign(yMin-yAtXmin)!=Math.sign(yMin-yAtXmax)){
      selectedPCP.push(a)
      leftInt = xMaxC-1
      }
      else if (Math.sign(yMax-yAtXmin)!=Math.sign(yMax-yAtXmax)){
      selectedPCP.push(a)
      leftInt = xMaxC-1
      }
      }
      else if (leftInt == xMinF && xMinF != xMaxF){
      var yLeftInt = dat[cNames[leftInt]]
      var yRightInt = dat[cNames[rightInt]]
      if ((Math.sign(yMin-yLeftInt)!=Math.sign(yMin-yRightInt) ||
      Math.sign(yMax-yLeftInt)!=Math.sign(yMax-yRightInt)) &&
      (Math.sign(yMin-yAtXmin) == Math.sign(yMin-yLeftInt) &&
      Math.sign(yMax-yAtXmin) == Math.sign(yMax-yLeftInt))){
      selectedPCP.push(a)
      leftInt = xMaxC-1
      }
      }
      else if (leftInt == xMaxF && xMinF != xMaxF){
      var yLeftInt = dat[cNames[leftInt]]
      var yRightInt = dat[cNames[rightInt]]
      if ((Math.sign(yMin-yLeftInt)!=Math.sign(yMin-yRightInt) ||
      Math.sign(yMax-yLeftInt)!=Math.sign(yMax-yRightInt)) &&
      (Math.sign(yMin-yAtXmax) == Math.sign(yMin-yRightInt) &&
      Math.sign(yMax-yAtXmax) == Math.sign(yMax-yRightInt))){
      selectedPCP.push(a)
      leftInt = xMaxC-1
      }
      }
      else if (leftInt >= xMin && rightInt <= xMax){
      var yLeftInt = dat[cNames[leftInt]]
      var yRightInt = dat[cNames[rightInt]]
      if (Math.sign(yMin-yLeftInt)!=Math.sign(yMin-yRightInt) ||
      Math.sign(yMax-yLeftInt)!=Math.sign(yMax-yRightInt)){
      selectedPCP.push(a)
      leftInt = xMaxC-1
      }
      }
      leftInt++;
      }
      }
      }
      var updateSPCP = []
      var selectedPCPL = selectedPCP.length
      for (a=0; a<selectedPCPL; a++){
      updateSPCP[a]=selectedPCP[a]+1
      }
      //console.log(selectedPCP)
      //Plotly.deleteTraces(el.id, updateSPCP);
      
      var update = {
      line: {
      color: 'blue',
      width: 1
      },
      opacity: 1,
      }
      Plotly.restyle(el.id, update, updateSPCP);
      Shiny.onInputChange('updateSPCP', updateSPCP);
      })
      }", data = list(pcpDat = pcpDat, nVar = nVar, colNms = colNms))})

    selID <- reactive(input$updateSPCP)
    output$selectedValues <- renderPrint({str(selID())})
    output$selectedPCP <- renderPrint({str(pcpDat[selID(),])})
    }
  
  shinyApp(ui, server)  
}
