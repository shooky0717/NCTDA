#' @title Visualize the expression levels of a specified gene across spots within the tissue
#' @description
#' visualize the expression levels of a specified gene across spatial locations within the tissue, providing an intuitive representation of spatial gene expression patterns.
#' @param object NCTDA object
#' @param gene_to_display character string claiming the gene to display.
#' @param scaled_gene_expression (default TURE) if TRUE, the gene expression levels
#' will be scaled to be the values between 0 and 1.
#' @param point_size (default 0.05) the point size on the pattern.
#' @param filling_colors (default c("white", "blue")) a 2-dimensional character string
#' vector to specify the "low" and "high" inputs respectively for "scale_fill_gradient()"
#' of ggplot2, which is to create a diverging colour gradient.
#' @param title_size (default 20) font size of the figure title.
#' @param legend_title_size (default 18) font size of the legend title.
#' @param legend_text_size (default 16) font size of the figure texts.
#' @return return NCTDA object.
#'
#' @import ggplot2
#'
#' @export

visualizeGenePattern <- function(object = NCTDA,
                                 gene_to_display = gene_name,
                                 celltype = 1,
                                 scaled_gene_expression = TRUE,
                                 point_size = 4,
                                 filling_colors = c("white", "darkblue"),
                                 title_size = 20, legend_title_size = 18, legend_text_size = 16
){
  
  if (gene_to_display %in% rownames(object@original_gene_expression) & length(gene_to_display) == 1) {
    gene_id <- which(rownames(object@original_gene_expression) %in% gene_to_display)
    
    gene_ex <- object@original_gene_expression[gene_id,]
    gene_ex <- gene_ex[intersect(rownames(object@location),
                                 colnames(object@original_gene_expression))]
    
    if (scaled_gene_expression) {
      gene_ex <- (gene_ex - min(gene_ex)) /(max(gene_ex) - min(gene_ex))
    }
    
    display_info <- data.frame(x = object@location[,1],
                               y = object@location[,2],
                               gene_expression = gene_ex)
    
    plot <- ggplot(display_info, aes(x, y)) +
      geom_point(aes(fill = gene_expression),
                 size = point_size, 
                 shape = 21,
                 stroke = 0.5 ) +
      theme_minimal() +
      scale_fill_gradient(low = filling_colors[1], high = filling_colors[2]) +
      theme(legend.title = element_text(size = legend_title_size),
            legend.key.height = unit(1, 'cm'),
            legend.key.width = unit(1, 'cm'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            plot.title = element_text(size = title_size),
            legend.text = element_text(size = legend_text_size),
            legend.position = "right",
            legend.title.align = 0.5)+
      labs(fill = " Expression level",
           title = paste0("Top gene in ", celltype,": ", gene_to_display)) +
      coord_fixed(xlim = c(0,1), ylim = c(0,1))
    
    print(plot)
  } else {
    stop('The gene you selected is not included in the dataset. Please select another!')
  }
}
