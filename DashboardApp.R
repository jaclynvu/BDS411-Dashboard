library(shiny)
library(reshape2)
library(corrplot)
library(ggplot2)
library(reshape2)
library(corrplot)

# Load the CSV files
pca_data <- read.csv("./AT_root_scRNAseq_PCA.csv")
cluster_data <- read.csv("./AT_root_scRNAseq_clustermap.csv")

# Merge datasets based on CellBarcode
merged_data <- merge(pca_data, cluster_data, by = "CellBarcode")

# Calculate cluster average PCA values
cluster_avg <- aggregate(. ~ seurat_clusters, data = merged_data[, -1], FUN = mean)
t_cluster <- t(cluster_avg)

# Calculate Pearson correlation
cor_matrix <- cor(t_cluster[, -1], method = "pearson")

# Convert correlation matrix to long format
cor_long <- melt(cor_matrix)

# Reshape correlation matrix to square form
cor_square <- acast(cor_long, Var1 ~ Var2, value.var = "value")

predicted <- read.csv("AT_root_scRNAseq_anticipatedmarkers.csv", header=T)
predicted <- as.data.frame(t(predicted))
predicted <- cbind(CellBarcode = rownames(predicted), predicted)
rownames(predicted) <- NULL

known <- read.csv("AT_root_scRNAseq_putativemarkers.csv")

umap <- read.csv("./AT_root_scRNAseq_UMAP.csv")


merged <- merge(umap, predicted, by="CellBarcode")
merged <- merge(merged, known, by="CellBarcode") #merge dfs for marker expression plot

# Define UI

ui <- fluidPage(
  titlePanel("Our Super Cool Baller BDS 411 Shiny App that is Fully Functional?"),
  tabsetPanel(
    tabPanel("UMAP vs PCA",
             fluidRow(
               column(6,
                      h2("UMAP"),
                      checkboxGroupInput("choice","Color Cell Cluster",0:20, selected = 0:20,
                                         inline=TRUE), #inline=TRUE makes it so it's not all one column
                      checkboxGroupInput("display", "Display Unselected Cell Clusters?",
                                         "Yes!", selected="Yes!")
                      
               ),
               column(6,
                      h2("PCA"),
                      h4('Please select the values for X and Y: PC_# vs PC_# (Numbers between 1 and 50)'),
                      numericInput(inputId = "pca_x", label = "For X values:", value = 1, min = 1, max = 50),
                      numericInput(inputId = "pca_y", label = "For Y values:", value = 1, min = 1, max = 50)
                      
               )
             ),
             fluidRow(
               column(6,
                      plotOutput("umap")
               ),
               column(6,
                      plotOutput("pca")
               )
             )
    ),
    tabPanel("Bar Chart & Heatmap",
             fluidRow(
               column(5,
                      h2("Bar Chart"),
                      h4("Counts the number of cells in each cell cluster"),
                      plotOutput("barchart")
                      
               ),
               column(5,
                      h2("Heatmap"),
                      h4("Pearson correlation in an all vs all comparison of cell clusters (based on cluster average PCA values)"),
                      plotOutput("heatmap")
               )
             )
    ),
    tabPanel("Predicted or Known Markers",
             column(2,
                    selectInput("category", "Select a category of markers:",
                                choices = c("Predicted", "Known")),
                    uiOutput("sub_category")
             ),
             column(8,
                    plotOutput("umap_expr")
             )
    )
  )
)

# server() defines the UMAP and PCA to pass to ui
server <- function(input, output, session) {
  
  #import ggplot2
  library (ggplot2)
  library(dplyr)
  
  #import data as "cell_pca"
  umap_df <- read.csv("./AT_root_scRNAseq_UMAP.csv")
  pca_df <- read.csv("./AT_root_scRNAseq_PCA.csv")
  
  cell_clusters <- read.csv("./AT_root_scRNAseq_clustermap.csv")
  
  #join clusters onto UMAP
  umap_clusters<-left_join(umap_df, cell_clusters, by='CellBarcode')
  pca_clusters <- left_join(pca_df, cell_clusters, by='CellBarcode')
  
  colors <- scales::hue_pal()(length(unique(umap_clusters$seurat_clusters))) #This gets you the default ggplot2 colormap
  names(colors) <- as.character(0:20) #Give each color a cluster ID it can be mapped to (like Python Dictionary)
  
  
  output$umap <- renderPlot({
    plot_mapping <- colors #Copy base color scale
    
    
    if(length(input$display) > 0){ #If the box is selected
      
      other_color <- "black"
      
    } else{
      
      other_color <- "transparent" #transparent makes the dots disappear
    }
    plot_mapping[!(names(colors) %in% input$choice)] <- other_color
    
    ggplot(data = umap_clusters, aes(UMAP_1, UMAP_2, color=as.factor(seurat_clusters))) +
      geom_point(alpha=0.6) +
      scale_color_manual(values = plot_mapping) +
      labs(color="Cell Cluster", x="UMAP Dimension 1", y="UMAP Dimension 2",
           title="UMAP Projection of Single Cell Expression Profiles") + 
      theme(text = element_text(size = 16), plot.title = element_text(hjust=0.5))
  })
  # Create a reactive expression for the PCA plot
  pca_plot <- reactive({
    
    plot_mapping <- colors #Copy base color scale
    
    
    if(length(input$display) > 0){ #If the box is selected
      
      other_color <- "black"
      
    } else{
      
      other_color <- "transparent" #transparent makes the dots disappear
    }
    plot_mapping[!(names(colors) %in% input$choice)] <- other_color
    x_val <- paste0("PC_", as.character(input$pca_x))
    y_val <- paste0("PC_", as.character(input$pca_y))
    pca_data <- data.frame (PC_x  = pca_clusters[x_val],
                            PC_y = pca_clusters[y_val],
                            seurat_clusters = pca_clusters['seurat_clusters']
    )
    pca_data$seurat_clusters <- as.factor(pca_data$seurat_clusters)
    
    ggplot(data = pca_data, aes_string(x = x_val, y = y_val, color='seurat_clusters')) +
      geom_point(alpha=0.6) + 
      labs(title = paste0("PCA plot of ", x_val, " vs ", y_val)) +
      scale_color_manual(values = plot_mapping) + 
      theme(text = element_text(size = 16), plot.title = element_text(hjust=0.5), legend.position="none")
  })
  
  # Render the PCA plot
  output$pca <- renderPlot({
    pca_plot()
  })
  
  # Update the plot when both dropdown menus are changed
  observeEvent((input$pca_x), {
    output$pca <- renderPlot({
      pca_plot()
    })
  })
  observeEvent(c(input$pca_y), {
    output$pca <- renderPlot({
      pca_plot()
    })
  })
  output$barchart <- renderPlot({
    
    clustermap <- read.csv("./AT_root_scRNAseq_clustermap.csv")
    grouping_df <- clustermap %>% group_by(seurat_clusters) %>% 
      summarise(total_count=n(), .groups = 'drop')
    # Barplot
    ggplot(grouping_df, aes(seurat_clusters, total_count, color = as.factor(seurat_clusters), fill = as.factor(seurat_clusters))) + 
      geom_bar(stat = "identity") + labs(title = "Number of Cells vs. Seurate Cell Cluster", x = "Seurate Cluster ID", y = "Number of Cells")
  })
  #new_pca <- PCA_clusters[-1]
  #new_pca %>%
  #group_by(seurat_clusters) %>%
  #summarise(across(everything(), mean), .groups = 'drop') %>%
  #head()
  output$heatmap <- renderPlot({
    # Generate heatmap using corrplot library
    corrplot(cor_square, method = "color", order = "hclust", tl.col = "black")
  })
  
  
  sub_categories <- reactiveValues(predicted = colnames(predicted)[-1],
                                   known = colnames(known)[-1])
  
  output$sub_category <- renderUI({
    if (input$category == "Predicted") {
      selectInput("predicted",
                  "Select a marker:", choices = sub_categories$predicted)
    } 
    else if (input$category == "Known") {
      selectInput("known", "Select a marker:", 
                  choices = sub_categories$known)
    }
  })
  
  
  output$umap_expr <- renderPlot({
    if (input$category == "Predicted") {
      ggplot(merged, aes_string(x = 'UMAP_1', y = 'UMAP_2', color=input$predicted, alpha=0.6)) +
        geom_point()
    } else if (input$category == "Known") {
      ggplot(merged, aes_string(x = 'UMAP_1', y = 'UMAP_2', color=input$known, alpha=0.6)) +
        geom_point()
    }
  })
  
}
# combine UI and server into Shiny
shinyApp(ui, server)
