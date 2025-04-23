library(shiny)
library(tidyverse)
library(ggplot2)
library(breastCancerVDX)
library(colourpicker)  
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(corrplot)
library(Rtsne)


# Load dataset
data(vdx)

# Clean column names - remove '' string 
colnames(vdx) <- str_remove_all(colnames(vdx), "'")

expression.values <- exprs(vdx)
features <- fData(vdx)
er.status <- pData(vdx)$er


# Define User Interface
ui <- fluidPage(
  titlePanel("GeneView: Gene Expression Analysis in Breast Cancer"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("geneFile", "Upload Custom Gene List (CSV or TXT)", accept = c(".csv", ".txt"), placeholder = "One gene symbol per line or in first column"),
      selectizeInput("genes", "Select Genes:", multiple = TRUE, choices = c("ESR1", "ERBB2", "PTEN", "BRCA1", "BRCA2", "ATM", "AAAS", "A4GNT", "AADAC", "AAK1")),
      colourInput("color_ER0", "Boxplot Color for ER- (0):", value = "green"),
      colourInput("color_ER1", "Boxplot Color for ER+ (1):", value = "blue"),
      colourInput("cluster_color", "Select cluster Color", value = "#2c7fb8"),
      checkboxInput("toggle_pdata", "Preview of Breast Cancer VDX Dataset", value = FALSE),
      radioButtons("dimred_method", "Select Dimensionality Reduction Method:",
                   choices = c("PCA", "t-SNE"), selected = "PCA"),
      actionButton("run_dimred", "Run Dimensionality Reduction"),
      downloadButton("downloadResults", "Download All Results")
    ),
    
    mainPanel(
      
      h3("Search for a gene in the Breast Cancer VDX dataset"),
      fluidRow(
        column(6,DT::dataTableOutput("allGenes")),
        column(6,tableOutput("erTable"))),
      
      conditionalPanel(condition = "input.toggle_pdata == true", verbatimTextOutput("head_pdata")),
      
      h3("Plot output for PCA/t-SNE"),
      plotOutput("dimred_plot"),
      
      h3("Correlation Plot of Selected Genes"),
      plotOutput("corPlot"),
      
      h3("Boxplots of Selected Gene Expression by ER Status"),
      plotOutput("boxPlot"),
      
      h3("Statistical Test (Auto-selected based on Normality)"),
      verbatimTextOutput("autoTestResult"),
      
      h3("Gene Ontology Enrichment Analysis Results"),
      verbatimTextOutput("goEnrichmentResults"),
      
      h3("Gene Ontology Enrichment Visualization"),
      plotOutput("goPlot")
      
    )
  )
)

# Define Server Logic
server <- function(input, output, session) {
  
  output$allGenes <- DT::renderDataTable({
    data.frame(Gene.Symbol = sort(unique(features$Gene.symbol)))
  })
  
  # Print the first few rows of pData(vdx) for reference 
  output$head_pdata <- renderPrint({ 
    head(pData(vdx))
  })
  
  
  # Read Uploaded Gene File with error handling 
  uploadedGenes <- reactive({
    req(input$geneFile)  # Ensure that the file is uploaded
    
    
    ext <- tools::file_ext(input$geneFile$name)  # Get the file extension
    tryCatch({
      if (ext == "csv") {  # If the file is a CSV
        # Read the CSV file and extract the first column of gene names
        genes <- read.csv(input$geneFile$datapath, header = TRUE)[, 1]
        if (length(genes) == 0) stop("CSV file is empty or does not contain genes.")
        return(genes)
      } else if (ext == "txt") {  # If the file is a TXT
        # Read the file as a list of gene names (one per line)
        genes <- readLines(input$geneFile$datapath)
        if (length(genes) == 0) stop("TXT file is empty or does not contain genes.")
        return(genes)
      } else {
        stop("Unsupported file format. Please upload a CSV or TXT file.")  # Handle unsupported formats
      }
    }, error = function(e) {  # Error handling for invalid file formats or read issues
      showModal(modalDialog(
        title = "File Error",
        paste("An error occurred while reading the file:", e$message),
        easyClose = TRUE,
        footer = NULL
      ))
      return(NULL)  # Return NULL if there is an error
    })
  })

  
  #Preserve Predefined Gene Choices as Fallback
  observe({
    genes <- uploadedGenes()
    if (length(genes) > 0) {
      updateSelectizeInput(session, "genes", choices = genes, server = TRUE)
    } else {
      updateSelectizeInput(session, "genes", choices = c("ESR1", "ERBB2", "PTEN", "BRCA1", "BRCA2", "ATM"))
    }
  })
  
  
  
  # Function to filter expression values based on user-selected genes with validation
  filterByExpr <- reactive({
    req(input$genes)# Ensures user has selected atleast 1 gene
    input_genes <- input$genes #Grabs the current gene list from user input (via selectize or File input)
    
    # Validate that input genes exist in the dataset-split into valid and invalid genes
    valid_genes <- input_genes[input_genes %in% features$Gene.symbol]
    invalid_genes <- setdiff(input_genes, valid_genes)
    
    if (length(valid_genes) == 0) {
      showNotification("None of the selected genes are present in the dataset.", type = "error")
      return(NULL)
    } else if (length(invalid_genes) > 0) {
      showNotification(paste("The following genes were not found and will be skipped:", 
                             paste(invalid_genes, collapse = ", ")), type = "warning")
    }
    
    # For valid genes, look up their corresponding probe IDs
    probe.ids <- features$probe[match(valid_genes, features$Gene.symbol)]
    # Extracts the expression matrix rows corresponding to the probe IDs.
    expression.values[probe.ids, , drop = FALSE]
  })
  
  # ER status count table
  output$erTable <- renderTable({
    er_counts <- table(er.status)
    data.frame(ER_Status = names(er_counts), Count = as.vector(er_counts))
  })
  
  ### run dimensionality reduction with a minimun of 10 genes
  
  # Define exprs_selected as a reactive expression
  exprs_selected <- reactive({
    req(filterByExpr())  # Ensure expression data is available
    return(filterByExpr())  # Return the filtered expression data
  })
  
  dimred_data <- eventReactive(input$run_dimred, {
    req(exprs_selected())
    
    # Require at least 10 genes
    validate(
      need(nrow(exprs_selected()) >= 10, "Please select at least 10 genes for clustering.")
    )
    
    #scale the expression data 
    exprs_scaled <- t(scale(t(exprs_selected())))
    n_samples <- ncol(exprs_scaled)
    
    if (input$dimred_method == "PCA") {
      pca <- prcomp(exprs_scaled, scale. = TRUE)
      df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                       Sample = rownames(pca$x))
      list(data = df, method = "PCA")
    } else {
      
      perplexity_val <- min(30, floor((n_samples - 1) / 3))
      
      if (perplexity_val < 1) {
        stop("Not enough samples to perform t-SNE. Please select more genes.")
      }
      tsne <- Rtsne(t(exprs_scaled), dims = 2, perplexity = perplexity_val, verbose = TRUE)
      df <- data.frame(Dim1 = tsne$Y[, 1], Dim2 = tsne$Y[, 2],
                       Sample = rownames(t(exprs_scaled)))
      list(data = df, method = "t-SNE")
    }
  })
  
  
  # Render clustering plot
  output$dimred_plot <- renderPlot({
    req(dimred_data())
    df <- dimred_data()$data
    method <- dimred_data()$method
    color3 <- input$cluster_color
    ggplot(df, aes_string(x = names(df)[1], y = names(df)[2])) +
      geom_point(color = color3, size = 3) +
      theme_minimal(base_size = 14) +
      labs(title = paste(method, "Sample Clustering"))
  })
  
  # Render heatmap
  output$heatmap_plot <- renderPlot({
    req(exprs_selected())
    pheatmap::pheatmap(exprs_selected(),
                       cluster_rows = TRUE, cluster_cols = TRUE,
                       show_colnames = FALSE, fontsize_row = 6)
  })
  

  # correlation plot 
  
  output$corPlot <- renderPlot({
    values <- filterByExpr()
    req(values)
    
    # Must have at least 2 genes to compute correlations
    if (nrow(values) < 2) {
      showNotification("Select at least 2 genes to view correlation plot.", type = "warning")
      return(NULL)
    }
    
    # Check gene symbols and their match in the features data
    selected_genes <- input$genes
    
    # Match the selected genes with the Gene symbols in features
    probe_ids <- match(selected_genes, features$Gene.symbol)# Get indices of matched genes

    # If any selected genes are not in the dataset, probe_ids will have NA
    if (any(is.na(probe_ids))) {
      showNotification("Some selected genes were not found in the dataset.", type = "error")
      return(NULL)
    }
    
    # Get the probe IDs corresponding to the valid genes
    probe_ids <- features$probe[probe_ids]

    
    # Subset the expression data by the matching probe IDs
    values <- values[probe_ids, , drop = FALSE]
    
    # Replace row names with gene symbols for readability
    rownames(values) <- selected_genes
    
    # Compute correlation between genes (rows are probes, so we transpose)
    corr_mat <- cor(t(values), method = "pearson")
    
    # Plot the correlation matrix
    ggcorrplot::ggcorrplot(corr_mat, 
                           lab = TRUE, 
                           type = "lower",
                           tl.cex = 10,
                           lab_size = 3)
  })
  

  
  
  # Box plot
  output$boxPlot <- renderPlot({
    expr_values <- filterByExpr()
    if (is.null(expr_values) || nrow(expr_values) == 0) return(NULL)
    
    # Convert ER status to factor for proper grouping
    er.factor <- as.factor(er.status)
    
    # Get matching gene symbols
    matched_symbols <- features$Gene.symbol[match(rownames(expr_values), features$probe)]
    
    # Create a long-format data frame for plotting
    plot_data <- data.frame(
      Expression = as.vector(t(expr_values)),
      ER_Status = rep(er.factor, each = nrow(expr_values)),
      Gene = rep(matched_symbols, times = ncol(expr_values))
    )
    
    # Get both selected colors
    color0 <- input$color_ER0
    color1 <- input$color_ER1
    fill_colors <- c("0" = color0, "1" = color1)
    
    ggplot(plot_data, aes(x = ER_Status, y = Expression, fill = ER_Status)) +
      geom_boxplot() +
      facet_wrap(~ Gene, scales = "free_y") +
      scale_fill_manual(values = fill_colors) +
      labs(x = "Estrogen Receptor Status (0 = negative, 1 = positive)", y = "Expression Level") +
      theme_minimal()
  })
  
  #Decide whether to run a t-test or Wilcoxon test based on the Shapiro-Wilk normality results
  output$autoTestResult <- renderPrint({
    values <- filterByExpr()
    req(values)
    
    # Check if there are selected genes
    selected_genes <- input$genes
    if (length(selected_genes) == 0) {
      return("No genes selected.")
    }
    
    # Loop through each selected gene and perform normality testing
    for (selected_gene in selected_genes) {
      cat("Shapiro-Wilk Test Results for Gene:", selected_gene, "\n\n")
      
      probe.id <- features$probe[match(selected_gene, features$Gene.symbol)]
      if (is.na(probe.id)) return("Gene not found in the dataset.\n")
      
      gene_values <- expression.values[probe.id, ]
      group1 <- gene_values[er.status == 0]
      group2 <- gene_values[er.status == 1]
      
      # Shapiro-Wilk test
      shapiro1 <- shapiro.test(group1)
      shapiro2 <- shapiro.test(group2)
      
      cat("ER- group p-value:", round(shapiro1$p.value, 4), "\n")
      cat("ER+ group p-value:", round(shapiro2$p.value, 4), "\n\n")
      
      if (shapiro1$p.value > 0.05 && shapiro2$p.value > 0.05) {
        # Use t-test
        test_result <- t.test(group1, group2)
        cat ("Both groups passed the normality test (p > 0.05).\n")
        cat ("Using t-test \n\n")
      } else {
        # Use Wilcoxon test
        test_result <- wilcox.test(group1, group2)
        cat("At least one group did not pass the normality test (p ≤ 0.05).\n")
        cat ("Using **Wilcoxon rank-sum test** \n\n")
      }
      print(test_result)
      
      
      # Interpret p-value in test_result
      cat("\nInterpretation:\n")
      cat ("-----------------------------------------------------------------\n")
      if (test_result$p.value < 0.001) {
        cat("Highly significant difference (p < 0.001)\n\n")
      } else if (test_result$p.value < 0.01) {
        cat("Very significant difference (p < 0.01)\n\n")
      } else if (test_result$p.value < 0.05) {
        cat("Statistically significant difference (p < 0.05)\n\n")
      } else {
        cat("No statistically significant difference (p ≥ 0.05)\n\n")
      }
    }
  })
  
  # GO Enrichment Analysis - After Statistical Testing
  output$goEnrichmentResults <- renderPrint({
    selected_genes <- input$genes
    gene_ids <- bitr(selected_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    
    if (nrow(gene_ids) == 0) {
      cat("No valid genes for GO enrichment analysis.\n")
      return(NULL)
    }
    
    ego <- enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
    
    if (is.null(ego) || nrow(ego) == 0) {
      cat("No enriched GO terms found.\n")
    } else {
      simplified <- ego@result[, c("ID", "Description", "p.adjust", "GeneRatio", "FoldEnrichment", "Count", "geneID")]
      print(head(simplified, 10))  # Show top 10 entries for readability
    }
  })
  
  # visualize enriched GO terms,
  output$goPlot <- renderPlot({
    selected_genes <- input$genes
    gene_ids <- bitr(selected_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    
    if (nrow(gene_ids) == 0) return(NULL)
    
    ego <- enrichGO(
      gene = gene_ids$ENTREZID,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05
    )
    
    if (is.null(ego) || nrow(ego) == 0) return(NULL)
    
    barplot(ego, showCategory = 10, title = "Top 10 GO Terms")
  })
  
  
  ## Download Results 
  ##### Download Handler
  output$downloadResults <- downloadHandler(
    filename = function() {
      paste0("gene_expression_results_", Sys.Date(), ".zip")  # Zip file with date in the filename
    },
    content = function(file) {
      # Create a temporary directory to store files
      temp_dir <- tempdir()
      setwd(temp_dir)
      
      
      # Extract necessary components from ExpressionSet
      expr_data <- exprs(vdx)
      pheno_data <- pData(vdx)
      feature_data <- fData(vdx)
      
      
      # Step 1: Prepare the CSV for Statistical Testing Results
      selected_genes <- input$genes
      print("Selected genes:")
      print(selected_genes)
      
      df <- vdx
      
      results <- lapply(selected_genes, function(gene) {
        cat("Processing gene:", gene, "\n")
        
        # Find probe(s) for the gene
        probes <- rownames(feature_data)[feature_data$Gene.symbol == gene]
        
        cat("Matched probes for", gene, ":", paste(probes, collapse = ", "), "\n")
        
        
        if (length(probes) == 0) {
          warning(paste("Gene", gene, "not found. Skipping."))
          return(data.frame(
            Gene = gene,
            Test_Used = NA,
            P.Value = NA,
            Mean.ERneg = NA,
            Mean.ERpos = NA
          ))
        }
        
        # Remove NA values from the matched probes
        probes <- probes[!is.na(probes)]
        cat("Filtered matched probes:", paste(probes, collapse = ", "), "\n")
        
        if (length(probes) == 0) {
          warning(paste("No valid probes found for gene", gene, ". Skipping."))
          return(data.frame(
            Gene = gene,
            Test_Used = NA,
            P.Value = NA,
            Mean.ERneg = NA,
            Mean.ERpos = NA
          ))
        }
        
        
        # Use first matching probe (or average across probes if preferred)
        #        expr_values <- expr_data[probes[1], ]
        probe <- probes[1]
        cat("Using probe:", probe, "\n")
        
        
        expr_values <- expr_data[probe, ]
        cat("Expression values (first 5):", paste(round(expr_values[1:5], 2), collapse = ", "), "\n")
        
        # Extract ER status
        group <- pheno_data$er
        cat("ER group summary:\n")
        print(table(group, useNA = "ifany"))
        
        if (is.null(group) || length(unique(group)) < 2) {
          warning("Insufficient ER group data. Skipping gene.")
          return(data.frame(
            Gene = gene,
            Test_Used = NA,
            P.Value = NA,
            Mean.ERneg = NA,
            Mean.ERpos = NA
          ))
        }
        
        
        group1 <- expr_values[group == 1]
        group2 <- expr_values[group == 0]
        cat("Group sizes: ER+ =", length(group1), ", ER− =", length(group2), "\n")
        
        
        # Shapiro test to check normality
        p1 <- shapiro.test(group1)$p.value
        p2 <- shapiro.test(group2)$p.value
        normal <- p1 > 0.05 && p2 > 0.05  # Check if both groups are normally distributed
        
        test_result <- if (normal) {
          t.test(group1, group2)  # t-test if normal
        } else {
          wilcox.test(group1, group2)  # Wilcoxon test if not normal
        }
        
        data.frame(
          Gene = gene,
          Test_Used = if (normal) "t-test" else "Wilcoxon",
          P.Value = test_result$p.value,
          Mean.ERpos = mean(group1, na.rm = TRUE),
          Mean.ERneg = mean(group2, na.rm = TRUE)
        )
      })
      
      result_df <- do.call(rbind, results)
      result_csv <- file.path(temp_dir, "statistical_test_results.csv")
      write.csv(result_df, result_csv, row.names = FALSE)
      print(paste("Saving results to:", result_csv))
      
      # Step 2: Prepare Gene Ontology (GO) Enrichment Analysis Results
      print(selected_genes)
      
      gene_ids <- bitr(selected_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
      
      # Check if gene IDs were successfully mapped
      if (nrow(gene_ids) == 0) {
        print("No gene IDs could be mapped. Please check the gene symbols or the OrgDb.")
        return(NULL)  # Exit if no valid gene IDs are found
      } else {
        print(paste("Mapped", nrow(gene_ids), "gene(s) to Entrez IDs."))
      }
      
      print(gene_ids)
      valid_gene_ids <- gene_ids$ENTREZID[!is.na(gene_ids$ENTREZID)]
      
      ego <- enrichGO(gene = valid_gene_ids, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
      
      # Check if GO enrichment results are empty
      if (length(ego) == 0 || nrow(ego) == 0) {
        print("No significant GO terms found. Check the enrichment parameters.")
        return(NULL)  # Exit if no enrichment results are found
      }
      
      if (nrow(ego) > 0) {
        print("GO enrichment results found.")
      } else {
        print("No GO enrichment results.")
      }
      
      
      go_results_df <- data.frame(
        ID = ego$ID,
        Description = ego$Description,
        p.adjust = ego$p.adjust,
        GeneRatio = ego$GeneRatio,
        FoldEnrichment = ego$FoldEnrichment,
        Count = ego$Count,
        geneID = ego$geneID
      )
      go_csv <- file.path(temp_dir, "go_enrichment_results.csv")
      write.csv(go_results_df, go_csv, row.names = FALSE)
      print(paste("Saving GO results to:", go_csv))
      
      # Step 3: Plot the GO Enrichment Barplot (Check data before plotting)
      if (nrow(go_results_df) > 0) {
        go_plot_file <- file.path(temp_dir, "go_enrichment_plot.png")
        png(go_plot_file)
        show_categories <- min(10, nrow(go_results_df))
        barplot(ego, showCategory = show_categories, title = "Top 10 GO Terms")
        dev.off()
        print(paste("GO enrichment barplot saved to:", go_plot_file))
      } else {
        print("No GO terms to plot. Skipping barplot.")
      }
      
      # Step 4 Get clsuter plots and heatmaps
 
       # Sample clustering plot (PCA/t-SNE)
       clustering_plot_path <- file.path(temp_dir, paste0("Sample_Clustering", ".png"))
       png(clustering_plot_path, width = 7, height = 5)
       print(ggplot(dimred_data()$data, aes_string(x = names(dimred_data()$data)[1],
                                                   y = names(dimred_data()$data)[2])) +
               geom_point(color = input$color, size = 3) +
               theme_minimal(base_size = 14) +
               labs(title = paste(dimred_data()$method, "Sample Clustering")))
       dev.off()
       
       # Heatmap plot
       heatmap_path <- file.path(temp_dir, "Gene_Expression_Heatmap.pdf")
       pdf(heatmap_path, width = 8, height = 10)
       pheatmap::pheatmap(exprs_selected(),
                          cluster_rows = TRUE, cluster_cols = TRUE,
                          show_colnames = FALSE, fontsize_row = 6)
       dev.off()
 
      
      # Step 5: Create Correlation Plot for Selected Genes
      selected_genes <- trimws(selected_genes)  # Remove leading/trailing spaces
      colnames(vdx) <- trimws(colnames(vdx))  # Remove spaces from column names in vdx
      selected_genes <- toupper(selected_genes)  # Convert to uppercase
      colnames(vdx) <- toupper(colnames(vdx))  # Convert column names to uppercase

      # Get gene symbols from features
      
      gene_symbols <- features$Gene.symbol  # Might also be 'GeneSymbol' or similar
      rownames(expression.values) <- gene_symbols
      
      # Remove double quotes from the selected genes
      selected_genes_clean <- gsub("\"", "", selected_genes)
      selected_genes_clean <- trimws(selected_genes)
      
      # Print cleaned list
      print("Selected genes (cleaned):")
      print(selected_genes_clean)
      
      # Check for available and missing genes
      available_genes <- selected_genes_clean[selected_genes_clean %in% rownames(expression.values)]
      missing_genes <- selected_genes_clean[!(selected_genes_clean %in% rownames(expression.values))]
      
      if(length(missing_genes) > 0) {
        writeLines(paste("Missing genes:", paste(missing_genes, collapse = ", ")), readme_file)
      }
      

      print("Checking which genes are available in rownames(vdx):")
      print(available_genes)
      
      # Then subsetting
      plot_data <- expression.values[available_genes, ]
      
      
      
      if (length(selected_genes_clean) >=3){
        correlation_plot_file <- file.path(temp_dir, "correlation_plot.png")
        png(correlation_plot_file)
        corr_matrix <- cor(t(plot_data), use = "complete.obs")
        corrplot(corr_matrix, method = "circle", title = "Correlation Plot of Selected Genes")
        dev.off()
        print(paste("Correlation plot saved to:", correlation_plot_file))
      } else {
        print("Not enough genes selected for correlation plot (need at least 3). Skipping correlation plot.")
      }
      
      
      
      # Step 6: Create Boxplot for Selected Genes by ER Status
      
      # Print cleaned list
      print("Selected genes (cleaned):")
      print(selected_genes_clean)
      
      # Check for available and missing genes
      available_genes <- selected_genes_clean[selected_genes_clean %in% rownames(expression.values)]
      missing_genes <- selected_genes_clean[!(selected_genes_clean %in% rownames(expression.values))]
      
      print("Missing genes:")
      print(missing_genes)
      
      print("Available genes:")
      print(available_genes)
      
      results <- lapply(available_genes, function(gene) {
        cat("Processing gene:", gene, "\n")
        
        # Find matching probes
        probes <- rownames(feature_data)[feature_data$Gene.symbol == gene]
        probes <- probes[!is.na(probes)]
        
        if (length(probes) == 0) {
          warning(paste("No probes found for gene:", gene))
          return(data.frame(
            Gene = gene,
            Test_Used = NA,
            P.Value = NA,
            Mean.ERneg = NA,
            Mean.ERpos = NA
          ))
        }
        
        probe <- probes[1]  # Just take the first matching probe
        expr_values <- as.numeric(expr_data[probe, ])
        group <- as.numeric(as.character(pheno_data$er))
        
        
        cat("ER group summary:\n")
        print(table(group, useNA = "ifany"))
        cat("Unique values in group:", unique(group), "\n")
        cat("Length of expr_values:", length(expr_values), "\n")
        cat("Length of group:", length(group), "\n")
        
        # Split expression values into two groups based on ER status
        group1 <- expr_values[group == 1]  # ER+ samples (group == 1)
        group2 <- expr_values[group == 0]  # ER- samples (group == 0)
        
        # Check dimensions of group1 and group2
        cat("Dimensions of group1 (ER+):", length(group1), "\n")
        cat("Dimensions of group2 (ER-):", length(group2), "\n")
        
        # Run appropriate test
        test_used <- NULL
        p_val <- NA
        
        if (length(group1) > 2 && length(group2) > 2) {
          if (shapiro.test(group1)$p.value > 0.05 && shapiro.test(group2)$p.value > 0.05) {
            test_used <- "t-test"
            p_val <- t.test(group1, group2)$p.value
          } else {
            test_used <- "wilcox"
            p_val <- wilcox.test(group1, group2)$p.value
          }
        }
        
        mean_neg <- mean(group2, na.rm = TRUE)
        mean_pos <- mean(group1, na.rm = TRUE)
        
        return(data.frame(
          Gene = gene,
          Test_Used = test_used,
          P.Value = p_val,
          Mean.ERneg = mean_neg,
          Mean.ERpos = mean_pos
        ))
      })
      # Combine results
      statistical_results <- do.call(rbind, results)
      write.csv(statistical_results, file.path(temp_dir, "statistical_test_results.csv"), row.names = FALSE)
      
      # Step 7: Create boxplot if >= 2 valid genes
      if (length(available_genes) >= 2) {
        expression_subset <- expression.values[available_genes, , drop = FALSE]
        
        plot_data <- data.frame(
          Expression = as.vector(expression_subset),
          Gene = rep(available_genes, each = ncol(expression.values)),
          ER_Status = rep(pheno_data$er, times = length(available_genes))
        )
        
        
        plot_data$ER_Status <- factor(plot_data$ER_Status, levels = c(0, 1), labels = c("ER-", "ER+"))
        
        cat("Preview of boxplot data:\n")
        print(head(plot_data))
        cat("Dimensions of plot_data:", dim(plot_data), "\n")
        
        boxplot_file <- file.path(temp_dir, "boxplot.png")
        png(boxplot_file, width = 900, height = 600)
        
        
        # Retrieve colors from user input
        er_minus_color <- input$color_ER0  # Color for ER-
        er_plus_color <- input$color_ER1   # Color for ER+
        
        print(
          ggplot(plot_data, aes(x = ER_Status, y = Expression, fill = ER_Status)) +
            geom_boxplot() +
            facet_wrap(~Gene, scales = "free_y") +
            scale_fill_manual(values = c("ER-" = er_minus_color, "ER+" = er_plus_color)) + # Apply the colors from input
            labs(title = "Expression by ER Status", x = "ER Status", y = "Expression Level") +
            theme_minimal()
        )
        dev.off()
        cat("Boxplot saved to:", boxplot_file, "\n")
        
      } else {
        cat("Not enough valid genes for boxplot (need at least 2). Skipping plot.\n")
      }
      
      # Step 8: Prepare the README File
      
      readme_text <- "
## GeneView: Gene Expression Analysis in Breast Cancer

This ZIP archive contains the results from your GeneView session, including statistical tests, Gene Ontology (GO) enrichment analysis, and various plots.

#### How to Download Results
To download the statistical test results and related files, open the app in **browser mode** (not in the Shiny app viewer) and click the **Download Results** button. This will generate a ZIP file containing the requested data.

#### Tips on Gene Selection 

For optimal results when using this app, we recommend selecting: >= 10 genes. Once selected, wait for a few seconds for the results to render.

#### Why this range?
Biological relevance: A focused subset of genes (e.g., ESR1, BRCA1, PGR, GATA3) helps interpret key pathways in estrogen receptor signaling.
Clean visualizations: Limits clutter in faceted boxplots, making trends easier to interpret.
Statistical reliability: Testing a moderate number of genes balances statistical power with reduced multiple testing burden.
App performance: Keeps rendering fast and responsive in interactive sessions.

- Choose genes with known roles in ER-positive or ER-negative breast cancer subtypes.
- Use literature, pathway databases, or your differential expression results to guide gene choices.

**For larger selection of genes (>20)**: 
- Skip faceted boxplots and use other visualizations like heatmaps or violin plots.
- Use dimensionality reduction (PCA, t-SNE) for overview.
- Apply multiple testing corrections (e.g., Benjamini-Hochberg FDR).

#### Files Included:
1. **statistical_test_results.csv**:
   - A table of results for each selected gene.
   - Columns:
     * Gene: Gene symbol
     * Test_Used: Statistical test applied (t-test if normal, Wilcoxon if not)
     * P.Value: p-value from the test
     * Mean.ERpos: Mean expression level in ER+ group
     * Mean.ERneg: Mean expression level in ER- group

2. **go_enrichment_results.csv**:
   - Results from Gene Ontology enrichment analysis, including the top enriched GO terms.
   - Columns:
     * ID: GO term ID
     * Description: Description of GO term
     * p.adjust: p-value adjusted for multiple testing
     * GeneRatio: Ratio of genes in the GO term to the total number of genes analyzed
     * FoldEnrichment: Fold enrichment for the GO term
     * Count: Number of genes in the GO term
     * geneID: List of genes in the GO term

3. **go_enrichment_plot.png**:
   - A barplot of the top 10 GO terms enriched in your gene set.

4. **correlation_plot.png**:
   - Correlation plot showing the relationships between the selected genes.

5. **boxplot.png**:
   - Boxplot comparing gene expression levels between ER+ and ER- groups for selected genes.
   
6. **Gene_Expression_Heatmap.pdf**:
   - A heatmap of gene expression across samples for the selected genes.
7. **Sample_Clustering.png**:
   - A plot showing hierarchical clustering of samples based on gene expression.
   
#### Statistical Test Notes:
- **Test Selection**: A Shapiro-Wilk test is performed on the gene expression values. If both groups (ER+ and ER-) pass normality (p > 0.05), a t-test is used. If either group fails normality, a Wilcoxon rank-sum test is applied.
- **P-Value Selection**: A p-value < 0.05 indicates a statistically significant difference between the two groups.
- **Null Hypothesis**: The null hypothesis for the t-test and Wilcoxon test is that there is no difference in expression between the two groups (ER+ vs ER-).
- **Alternative Hypothesis**: The alternative hypothesis is that there is a significant difference in gene expression between the ER+ and ER- groups.

Thank you for using GeneView! For questions, contact: pk563@snu.edu.in
 "
      
      readme_file <- file.path(temp_dir, "README.txt")
      writeLines(readme_text, readme_file)
      print("README file saved.")
      
      # List of files to include in the zip
      files_to_zip <- c(
        "statistical_test_results.csv",
        "go_enrichment_results.csv",
        "go_enrichment_plot.png",
        "correlation_plot.png",
        "boxplot.png",
        "Gene_Expression_Heatmap.pdf",
        "Sample_Clustering.png",
        "README.txt"
      )
      
      # Full paths in temp_dir
      full_paths <- file.path(temp_dir, files_to_zip)
      
      # Only include files that actually exist
      existing_files <- full_paths[file.exists(full_paths)]
      
      if (file.exists(boxplot_file)) {
        cat("Boxplot.png exists and will be added to the zip.\n")
      } else {
        cat("Boxplot.png does not exist. Skipping it in the zip.\n")
      }
      
      # Create the zip archive (flat structure using -j flag)
      utils::zip(zipfile = file, files = existing_files, flags = "-j")
      
      message("ZIP file created successfully.")
    }  
  )
}  

# Run the Shiny App
shinyApp(ui, server)   
