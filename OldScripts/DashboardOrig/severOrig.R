library(shinydashboard)
library(mlxR)

shinyServer(function(input, output){
  r <- reactive({ 
    p <- c(ka=input$ka, Tk0=input$Tk0, alpha=input$alpha, 
           F0=input$F0, Vm=input$Vm,  Km=input$Km, V=1, k=input$k)
    
    t.value=seq(0,24,length.out=241)
    out <- list(name=c("C1", "C2", "C3", "C4", "C5", "C6"), time=t.value)
    
    t1=input$tfd
    t2=input$ii*(input$nd-1)+t1
    if (t2>=t1){
      t.dose=seq(t1,t2,by=input$ii)
      adm <- list(time=t.dose, amount=input$amt)
    }else{
      adm <- list(time=t1, amount=0)
    }
    #----------------------------------------------------  
    res <- simulx(model     = "absorptionModel.txt", 
                  parameter = p, 
                  output    = out, 
                  treatment = adm)
    #----------------------------------------------------    
    return(res)
  })
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]}
  vc=gg_color_hue(6)[c(1,4,2,3,6,5)]
  names(vc)=letters[1:6]
  lc <- c("first order", "zero order", "alpha order",
          "sequential 0-1", "simultaneous 0-1", "saturated")
  names(vc)=letters[1:6]
  
  output$plot <- renderPlot({
    r <- r()
    npl <- 0
    vdisp <- rep(FALSE,6)
    pl=ggplotmlx()
    if (input$first==TRUE){
      pl=pl + geom_line(data=r$C1, aes(x=time, y=C1, colour="a"), size=1)  
      vdisp[1] <- TRUE}
    if (input$zero==TRUE){
      pl=pl + geom_line(data=r$C2, aes(x=time, y=C2, colour="b"), size=1) 
      vdisp[2] <- TRUE}
    if (input$al==TRUE){
      pl=pl + geom_line(data=r$C3, aes(x=time, y=C3, colour="c"), size=1)  
      vdisp[3] <- TRUE}
    if (input$sequential==TRUE){
      pl=pl + geom_line(data=r$C4, aes(x=time, y=C4, colour="d"), size=1)  
      vdisp[4] <- TRUE}
    if (input$mixed==TRUE){
      pl=pl + geom_line(data=r$C5, aes(x=time, y=C5, colour="e"), size=1)  
      vdisp[5] <- TRUE}
    if (input$saturated==TRUE){
      pl=pl + geom_line(data=r$C6, aes(x=time, y=C6, colour="f"), size=1)  
      vdisp[6] <- TRUE}
    pl <- pl + ylab("Amount = Concentration (V=1)") 
    pl <- pl + scale_colour_manual(values=vc[vdisp], labels=lc[vdisp])
    if (input$legend==TRUE){
      pl <- pl + theme(legend.position=c(.65, 0.95), legend.justification=c(0,1), legend.title=element_blank())
    }else{
      pl <- pl + theme(legend.position="none")
    }   
    if (input$nd==1)
      pl <- pl +ylim(c(0,input$amt))
    print(pl)
  }) 
  
  rr <- reactive({
    r <- r()
    rr <- r[[1]]
    for (k in (2:6))
      rr <- merge(rr, r[[k]])
    return(rr)
  })
  
  output$table <- renderTable({ 
    rr()
  })
  
  output$downloadTable <- downloadHandler(
    filename = "table.csv",
    content = function(file) {
      write.csv(rr(), file)
    }
  )
  
})