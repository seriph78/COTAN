
#' scatter.plot
#' 
#' Plot that check the realtion between the library size and the number gene 
#' detected. 
#'
#' @param obj COTAN object
#' @param split.pattern pattern used to extract, from the column names, the 
#' sample field (default " ")
#' @param n.col ones the column names are splitted by split.pattern, the column 
#' number with the sample name (default 2)
#' @param split.samples T/F if we want to plot each sample in a different
#' panel (default F)
#'
#'@importFrom scales trans_breaks
#'@importFrom scales math_format
#'@importFrom scales trans_format
#'@importFrom Matrix colSums
#'@importFrom ggthemes theme_tufte
#'@import ggplot2
#'@importFrom stringr str_split
#' @return
#' @export
#'
#' @examples
scatter.plot <- function(obj,split.pattern = " ", n.col=2, split.samples = F){
  
  lib.size <- Matrix::colSums(obj@raw)
  gene.size <- Matrix::colSums(obj@raw > 0)
  to.plot <- cbind(lib.size, gene.size)
  to.plot <- as.data.frame(to.plot)
<<<<<<< HEAD
  to.plot$sample <- stringr::str_split(rownames(to.plot),pattern = split.pattern,simplify = T)[,n.col]
  
  plot <- ggplot(to.plot, aes(x=lib.size,y=gene.size, color = sample)) + geom_point(size = 0.5, alpha= 0.8)+
    ggthemes::theme_tufte()+
=======
  to.plot$sample <- str_split(rownames(to.plot),pattern = split.pattern,simplify = T)[,n.col]
  
  plot <- ggplot(to.plot, aes(x=lib.size,y=gene.size, color = sample)) + geom_point(size = 0.5, alpha= 0.8)+
    theme_tufte()+
>>>>>>> In cotan-class changed a few fields. Added documentation. In houseKeepingGenes() removed a useless line and added a check for the nCells field in the object.
    labs(title = "Scatter plot of library size VS gene detected for each cell", 
         y = "Gene number",
         x = "Library size (UMI)") + 
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)))
  
  if(split.samples == T){
    plot <- plot + facet_grid(cols = vars(sample))
  }
  
  return(plot)
<<<<<<< HEAD
}
=======
}
>>>>>>> In cotan-class changed a few fields. Added documentation. In houseKeepingGenes() removed a useless line and added a check for the nCells field in the object.
