#' genes.size.plot
#' 
#' Function that plots the row gene number (reads greather then 0) for each cell and sample (if there 
#' are more samples).
#'
#' @param obj COTAN object
#' @param split.pattern pattern used to extract, from the column names, the 
#' sample field (default " ")
#' @param n.col ones the column names are splitted by split.pattern, the column 
#' number with the sample name (default 2)
#' 
#'@importFrom Matrix colSums
#'@import ggplot2
#'@importFrom stringr str_split
#' 
#' @return the violin-boxplot plot
#' @export
#'
#' @examples
genes.size.plot <- function(obj, split.pattern = " ", n.col=2){
  sizes <- Matrix::colSums(obj@raw > 0)
  sizes <- sort(sizes)
  sizes <- as.data.frame(sizes)
  sizes$n <- c(1:dim(sizes)[1])
  sizes$sample <- str_split(rownames(sizes),pattern = split.pattern,simplify = T)[,n.col]
  
  plot <- 
    sizes %>%
    ggplot(aes(x = sample, y = sizes, fill = sample)) +
    #geom_point(size = 0.5) +
    geom_flat_violin(position = position_nudge(x = .15, y = 0),adjust =2, alpha = 0.5) +
    geom_point(position = position_jitter(width = .1), size = .4, color="black", alpha= 0.5) +
    geom_boxplot(aes(x = sample, y = sizes, fill = sample),
                 outlier.shape = NA, alpha = .8, width = .15, colour = "gray65", size = 0.6)+
    labs(title = "Detected gene number", 
         y = "Size (number of genes)",
         x = "") +
    scale_y_continuous(expand = c(0, 0)) +
    ylim(0,max(sizes$sizes))+
    plotTheme("size-plot")

  return(plot)
  
}