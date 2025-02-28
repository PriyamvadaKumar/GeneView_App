
#Load required libraries
library(shiny)
library(tidyverse)

#Define UI
ui <- fluidPage(
  

  titlePanel("Analyzing VDX breast cancer dataset"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("thegene","Gene to Analyse",
                  choices=c("BRCA1", "BRCA2", "ATM","ESR1","ERBB2","PTEN"),
                  selected  = "ESR1"),
      radioButtons("colour","Colour of histogram",choices=c("orange","green","blue","pink", "yellow","red"),selected = "blue"),
      downloadButton("plotPDF",label = "Download"),
      a("dataset source", href="10.18129/B9.bioc.breastCancerVDX") ,
      
    ),
    
 
    mainPanel(
      plotOutput("boxplot"),
      verbatimTextOutput("testResult")
    )
  )
)


#Define server
server <- function(input, output) {
  
  library(breastCancerVDX)
  library(Biobase)
  library(tools)
  library(edgeR)

  
  #Load Breast cancer dataset
  data(vdx)
  expression.values <- exprs(vdx)
  features <- fData(vdx)
  er.status <- pData(vdx)$er


#Define Downloadable Plot Functionality

  output$plotPDF <- downloadHandler(
    filename = "boxplot.pdf",
    content = function(file){
      pdf(file)
      
      # Define Reactive Function to Filter Gene Expression Data
      filterByExpr <- reactive({
        gene <- input$thegene
        probe.id <- as.character(features$probe[match(gene, features$Gene.symbol)])
        Sys.sleep(10)
        expression.values[probe.id,]
      }) 
      
      #Render the Boxplot
      output$boxplot <- renderPlot({
        
        values <- filterByExpr()
        
        
        boxplot(values ~ er.status,col = input$colour)
        
      }) #output boxplot ends
      
      
      
       #Render Statistical Test Results
      output$testResult <- renderPrint(
        { 
          values <- filterByExpr()
          t.test(values ~ er.status)
          wilcox.test(values ~ er.status, alternative = "two.sided")
        }
      ) 

      #Close the PDF Output
      dev.off()
    }
    
  )  
  
}

# Run the application 
shinyApp(ui = ui, server = server)

