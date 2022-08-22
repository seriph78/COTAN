
#' mitochondrial.percentage.plot
#' 
#' Function that plots the row library size for each cell and sample (if there 
#' are more samples).
#'
#' @param obj COTAN object
#' @param split.pattern pattern used to extract, from the column names, the 
#' sample field (default " ")
#' @param n.col ones the column names are splitted by split.pattern, the column 
#' number with the sample name (default 2)
#' @param gene.prefix Prefix for the mitochondrial genes (default "^MT-" for 
#' Human, mouse "^mt-")
#' 
#'@importFrom Matrix colSums
#'@import ggplot2
#'@importFrom stringr str_split
#'@importFrom ggthemes theme_tufte
#' @return the violin-boxplot plot
#' @export
#'
#' @examples
mitochondrial.percentage.plot <- function(obj, split.pattern = " ", n.col=2, gene.prefix = "^MT-"){
  sizes <- Matrix::colSums(obj@raw)
  sizes <- as.data.frame(sizes)
  sizes$n <- c(1:dim(sizes)[1])
  sizes$sample <- stringr::str_split(rownames(sizes),
                                     pattern = split.pattern,simplify = T)[,n.col]
  mit.genes <- as.data.frame(obj@raw[rownames(obj@raw) %in% 
<<<<<<< HEAD
                         rownames(obj@raw)[stringr::str_detect(rownames(obj@raw),pattern = "^MT-")],])
  #colnames(mit.genes) <- rownames(obj@raw)[stringr::str_detect(rownames(obj@raw),pattern = "^MT-")]
  if(!identical(colnames(mit.genes),rownames(sizes))){
    print("Problem cell oreder!")
  }
  mit.genes <- t(mit.genes)
=======
                         rownames(obj@raw)[str_detect(rownames(obj@raw),pattern = "^MT-")],])
  colnames(mit.genes) <- rownames(obj@raw)[str_detect(rownames(obj@raw),pattern = "^MT-")]
  if(!identical(rownames(mit.genes),rownames(sizes))){
    print("Problem cell oreder!")
  }
>>>>>>> In cotan-class changed a few fields. Added documentation. In houseKeepingGenes() removed a useless line and added a check for the nCells field in the object.
  sizes$sum.mit <- Matrix::rowSums(mit.genes)
  sizes$mit.percentage <- round(sizes$sum.mit/sizes$sizes*100,digits = 2)
  
  plot <- 
    sizes %>%
    ggplot(aes(x = sample, y = mit.percentage, fill = sample)) +
    #geom_point(size = 0.5) +
    geom_flat_violin(position = position_nudge(x = .15, y = 0),adjust =2, alpha = 0.5) +
    geom_point(position = position_jitter(width = .1), size = .4, color="black", alpha= 0.5) +
    geom_boxplot(aes(x = sample, y = mit.percentage, fill = sample),
                 outlier.shape = NA, alpha = .8, width = .15, colour = "gray65", size = 0.6)+
    labs(title = "Mitochondrial percentage of reads", 
         y = "% (mit. reads / tot reads * 100)",
         x = "") +
    scale_y_continuous(expand = c(0, 0)) +
    #ylim(0,max(sizes$sizes))+
<<<<<<< HEAD
    ggthemes::theme_tufte()+
=======
    theme_tufte()+
>>>>>>> In cotan-class changed a few fields. Added documentation. In houseKeepingGenes() removed a useless line and added a check for the nCells field in the object.
    theme(legend.position = "none")#,
  #axis.text.x=element_blank(),
  #axis.ticks.x=element_blank())
  
  return(list("plot" = plot,"sizes"= sizes))
   
}


  
