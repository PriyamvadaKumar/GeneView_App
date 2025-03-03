library(shiny)
library(tidyverse)
library(ggplot2)
library(breastCancerVDX)

# Load dataset
data(vdx)
expression.values <- exprs(vdx)
features <- fData(vdx)
er.status <- pData(vdx)$er

# Define User Interface
ui <- fluidPage(
  titlePanel("GeneView: Gene Expression Analysis in Breast Cancer"),
  
  sidebarLayout(
    sidebarPanel(
      selectizeInput("genes", "Select Genes:", multiple = TRUE, choices = c("ESR1", "ERBB2", "PTEN", "BRCA1", "BRCA2", "ATM")),
      selectInput("gene_normality", "Select Gene for Normality Check:", choices = NULL),
      radioButtons("colour", "Colour of Boxplot:", choices = c("orange", "green", "blue"), selected = "orange"),
      downloadButton("downloadPlot", "Download Plot"),
      downloadButton("downloadResults", "Download Results")
    ),
    
    mainPanel(
      h3("Normality Test Results"),
      verbatimTextOutput("normalityResults"),
      
      h3("Normality Plots"),
      plotOutput("normalityPlot"),
      
      h3("Gene Expression Boxplot"),
      plotOutput("boxPlot"),
      
      h3("Per-Gene Statistical Test Results"),
      verbatimTextOutput("testResults"),
      
      h3("Group-Level Test for Selected Genes"),
      verbatimTextOutput("groupTestResults")
    )
  )
)

# Define Server Logic
server <- function(input, output, session) {
  
  # Function to filter expression values based on user-selected genes
  filterByExpr <- reactive({
    genes <- input$genes
    probe.ids <- features$probe[match(genes, features$Gene.symbol)]
    expression.values[probe.ids, , drop = FALSE]
  })
  
  observe({
    updateSelectInput(session, "gene_normality", choices = input$genes, selected = input$genes[1])
  })
  
  #check er status count
  output$erTable <- renderTable({
    # Count samples in each ER group
    er_counts <- table(er.status)
    
    # Convert table to data frame for display
    data.frame(ER_Status = names(er_counts), Count = as.vector(er_counts))
  })
  
  
  # Perform Shapiro-Wilk normality test
  output$normalityResults <- renderPrint({
    values <- filterByExpr()
    if (is.null(values) || nrow(values) == 0) return("No data available.")
    
    normality_results <- sapply(1:nrow(values), function(i) {
      shapiro.test(values[i, ])$p.value
    })
    gene_names <- features$Gene.symbol[match(rownames(values), features$probe)]
    names(normality_results) <- gene_names
    return(normality_results)
  })
  
  # Normality Plot
  output$normalityPlot <- renderPlot({
    values <- filterByExpr()
    if (is.null(values) || nrow(values) == 0) return(NULL)
    
    selected_gene <- input$gene_normality
    probe.id <- features$probe[match(selected_gene, features$Gene.symbol)]
    if (is.na(probe.id)) return(NULL)  # Avoid errors if gene is missing
    
    gene_values <- expression.values[probe.id, ]
    
    par(mfrow = c(1, 2))  
    hist(gene_values, probability = TRUE, main = paste("Histogram of", selected_gene), col = "lightblue", breaks = 20,  xlab = "Expression Level")
    curve(dnorm(x, mean = mean(gene_values), sd = sd(gene_values)), add = TRUE, col = "red", lwd = 2)
    lines(density(gene_values), col = "blue", lwd = 2)
    qqnorm(gene_values, main = paste("QQ Plot of", selected_gene))
    qqline(gene_values, col = "red")
  })
  
  
  #Box plot
  output$boxPlot <- renderPlot({
    values <- filterByExpr()
    if (is.null(values) || nrow(values) == 0) return(NULL)
    
    selected_gene <- input$gene_normality
    probe.id <- features$probe[match(selected_gene, features$Gene.symbol)]
    if (is.na(probe.id)) return(NULL)
    
    gene_values <- expression.values[probe.id, ]
    
    # Convert ER status to factor for proper grouping
    er.factor <- as.factor(er.status)
    
    # Create data frame for plotting
    plot_data <- data.frame(Expression = gene_values, ER_Status = er.factor)
    
    # Generate boxplot
    ggplot(plot_data, aes(x = ER_Status, y = Expression, fill = ER_Status)) +
      geom_boxplot() +
      scale_fill_manual(values = c("orange", "green", "blue")) + 
      labs(title = paste("Expression of", selected_gene, "by ER Status"),
           x = "ER Status", y = "Expression Level") +
      theme_minimal()
  })
  
  
  # Download Plots
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste0("normality_plot_", input$gene_normality, ".png")
    },
    content = function(file) {
      png(file)  
      
      values <- filterByExpr()
      if (is.null(values) || nrow(values) == 0) {
        dev.off()
        return()
      }
      
      selected_gene <- input$gene_normality
      probe.id <- features$probe[match(selected_gene, features$Gene.symbol)]
      if (is.na(probe.id)) {
        dev.off()
        return()
      }
      
      gene_values <- expression.values[probe.id, ]
      
      par(mfrow = c(1, 2))
      hist(gene_values, probability = TRUE, main = paste("Histogram of", selected_gene), col = "lightblue", breaks = 20)
      curve(dnorm(x, mean = mean(gene_values), sd = sd(gene_values)), add = TRUE, col = "red", lwd = 2)
      qqnorm(gene_values, main = paste("QQ Plot of", selected_gene))
      qqline(gene_values, col = "red")
      
      dev.off()  
    }
  )
  
  
  # Statistical Test Results
  output$testResults <- renderPrint({
    values <- filterByExpr()
    if (is.null(values) || nrow(values) == 0) return("No data available.")
    
    results <- list()
    probe_ids <- rownames(values)
    gene_names <- features$Gene.symbol[match(probe_ids, features$probe)]
    
    for (i in 1:nrow(values)) {
      gene_name <- gene_names[i]
      if (is.na(gene_name)) next 
      
      shapiro_result <- shapiro.test(values[i, ]) 
      test_used <- if (shapiro_result$p.value > 0.05) "Student's t-test" else "Wilcoxon rank-sum test"
      test_result <- if (shapiro_result$p.value > 0.05) t.test(values[i, ] ~ er.status) else wilcox.test(values[i, ] ~ er.status)
      
      results[[gene_name]] <- list("Test Used" = test_used, "Raw p-value" = test_result$p.value)
    }
    
    raw_p_values <- sapply(results, function(x) x$`Raw p-value`)
    adj_p_values <- p.adjust(raw_p_values, method = "BH")
    for (i in seq_along(results)) {
      results[[i]]$`Adjusted p-value` <- adj_p_values[i]
    }
    
    return(results)
  })
  
  # Group-Level Wilcoxon Test
  output$groupTestResults <- renderPrint({
    values <- filterByExpr()
    if (is.null(values) || nrow(values) == 0) return("No data available.")
    
    group_values <- colMeans(values)
    
    if (length(group_values) != length(er.status)) {
      return("Error: Mismatch in sample size between gene expression and ER status.")
    }
    
    test_result <- wilcox.test(group_values ~ er.status)
    list("Wilcoxon Test for Combined Gene Set" = test_result$p.value)
  })
#}

# Download Normality Test and Statistical Results as CSV

output$downloadResults <- downloadHandler(
  filename = function() {
    paste0("gene_expression_results.csv")
  },
  content = function(file) {
    values <- filterByExpr()
    if (is.null(values) || nrow(values) == 0) return()
    
    normality_results <- sapply(1:nrow(values), function(i) {
      shapiro.test(values[i, ])$p.value
    })
    
    gene_names <- features$Gene.symbol[match(rownames(values), features$probe)]
    
    # Construct results data frame
    results_df <- data.frame(Gene = gene_names, Normality_p_value = normality_results, stringsAsFactors = FALSE)
    write.csv(results_df, file, row.names = FALSE)  
  }
 )
}

# Run the Shiny App
shinyApp(ui = ui, server = server)

