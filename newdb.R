#

library(shiny)
library(tidyverse)
library(ggplot2)

ui <- fluidPage(
  

  titlePanel("Analyzing VDX breast cancer dataset"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("thegene","Gene to Analyse",
                  choices=c("ESR1","ERBB2","PTEN"),
                  selected  = "ESR1"),
      radioButtons("colour","Colour of histogram",choices=c("orange","green","blue"),selected = "orange"),
      downloadButton("plotPDF",label = "Download"),
      a("dataset source", href="10.18129/B9.bioc.breastCancerVDX") ,
      
    ),
    
 
    mainPanel(
      plotOutput("boxplot"),
      verbatimTextOutput("testResult")
    )
  )
)



server <- function(input, output) {
  
  library(breastCancerVDX)
  library(Biobase)
  library(tools)
  library(edgeR)
  
  data(vdx)
  expression.values <- exprs(vdx)
  features <- fData(vdx)
  er.status <- pData(vdx)$er
  
  output$plotPDF <- downloadHandler(
    filename = "boxplot.pdf",
    content = function(file){
      pdf(file)
      
   
      filterByExpr <- reactive({
        gene <- input$thegene
        probe.id <- as.character(features$probe[match(gene, features$Gene.symbol)])
        Sys.sleep(10)
        expression.values[probe.id,]
      }) 
      
      output$boxplot <- renderPlot({
        
        values <- filterByExpr()
        
        
        boxplot(values ~ er.status,col = input$colour)
        
      }) #output boxplot ends
      
      
      
      
      output$testResult <- renderPrint(
        { 
          values <- filterByExpr()
          t.test(values ~ er.status)
          wilcox.test(values ~ er.status, alternative = "two.sided")
        }
      ) #render print ends
      
      dev.off()
    }
    
  )  
  
}

# Run the application 
shinyApp(ui = ui, server = server)

