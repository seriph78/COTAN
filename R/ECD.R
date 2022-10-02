

#' ECD.plot
#' 
#' This function plots the empirical distribution function of library sizes  
#' (UMI number). It helps to define where to drop "cells" that are simple 
#' background signal.
#'
#' @param obj COTAN object
#' @param y_cut y threshold of library size to drop
#'
#' @return
#' @export
#'
#' @examples
ECD.plot <- function(obj, y_cut){
  lib.size <- colSums(obj@raw)
  lib.size <- sort(lib.size,decreasing = T)
  lib.size <- as.data.frame(lib.size)
  lib.size$n <- c(1:length(lib.size$lib.size))
  ggplot(lib.size, aes(y=log(lib.size), x=log(n)))+geom_point()+
    geom_hline(yintercept=log(y_cut), linetype="dashed",color = "red")

}
