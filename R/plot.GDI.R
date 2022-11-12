#' plot_GDI
#'
#' This function directly evaluate and plot the GDI for a sample.
#'
#' @param object A COTAN object
#' @param cond A string corresponding to the condition/sample (it is used only for the title). Default is empty.
#' @param type Type of statistic to be used. Default is "S":
#' Pearson's chi-squared test statistics. "G" is G-test statistics
#' @param genes a named list of genes to label. Each array will have different color. Default is empty.
#' @param GDI.df when the GDI data frame was already calculated, it can be put here to speed up the process. Default is NULL.
#' @return A ggplot2 object
#' @export
#' @import ggplot2
#' @import RColorBrewer
#' @import ggrepel
#' @importFrom  stats quantile
#' @importFrom Matrix forceSymmetric
#' @rdname plot_GDI
#' @examples
#' data("ERCC.cotan")
#' plot_GDI(ERCC.cotan, cond = "ERCC")
setGeneric("plot_GDI", function(object, cond = NULL,genes = NULL,type="S", GDI.df = NULL) standardGeneric("plot_GDI"))
#' @rdname plot_GDI
setMethod("plot_GDI","scCOTAN",
          function(object, cond, genes, type="S",GDI.df) {
            ET <- sum.raw.norm <- NULL
            
            if(is(class(object@coex)[1], "dtCMatrix")){
              print("COTAN object in the old format! Converting...")
              object <- get.coex(object)
              print(paste0("Saving as new file as ", dir, ET, "new.cotan.RDS"))
              saveRDS(object,paste0(dir, ET, "new.cotan.RDS"))
              
            }
            
            print("GDI plot ")
            if (is.null(GDI.df)) {
              if (type=="S") {
                GDI <- get.GDI(object,type="S")
              }else if(type=="G"){
                print("Using G")
                GDI <- get.GDI(object,type="G")
              }
              
            }else{
              GDI <- GDI.df
            }
            
            text.size <- 10
            GDI$colors <- "none"
            for (n in names(genes)) {
              GDI[rownames(GDI) %in% genes[[n]],]$colors <- n
              
            }
            
            qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
            col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
            
            mycolours <- col_vector[seq_along(names(genes))]
            names(mycolours) <- names(genes)
            
            base.plot <- ggplot(subset(GDI,colors == "none" ), aes(x=sum.raw.norm, y=GDI)) +  geom_point(alpha = 0.3, color = "#8491B4B2", size=2.5)
            
            textdf <- GDI[!GDI$colors == "none",]
            
            #############
            
            GDI_plot = base.plot +  geom_point(data = subset(GDI,colors != "none"  ), aes(x=sum.raw.norm, y=GDI, colour=colors),size=2.5,alpha = 0.8) +
              geom_hline(yintercept=quantile(GDI$GDI)[4], linetype="dashed", color = "darkblue") +
              geom_hline(yintercept=quantile(GDI$GDI)[3], linetype="dashed", color = "darkblue") +
              geom_hline(yintercept=1.5, linetype="dotted", color = "red", size= 0.5) +
              scale_color_manual("Status", values = mycolours)  +
              scale_fill_manual("Status", values = mycolours)  +
              xlab("log normalized counts")+ylab("GDI")+
              geom_label_repel(data =textdf , aes(x=sum.raw.norm, y=GDI, label = rownames(textdf),fill=colors),
                               label.size = NA,max.overlaps = 40,
                               alpha = 0.8,
                               direction ="both",
                               na.rm=TRUE,
                               seed = 1234) +
              ggtitle(paste("GDI plot ", cond))+
              theme(axis.text.x = element_text(size = text.size, angle = 0, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF" ),
                    axis.text.y = element_text( size = text.size, angle = 0, hjust = 0, vjust = .5, face = "plain", colour ="#3C5488FF"),
                    axis.title.x = element_text( size = text.size, angle = 0, hjust = .5, vjust = 0, face = "plain", colour ="#3C5488FF"),
                    axis.title.y = element_text( size = text.size, angle = 90, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF"),
                    legend.title = element_blank(),
                    legend.text = element_text(color = "#3C5488FF",face ="italic" ),
                    legend.position = "right")
            
            return(GDI_plot)
            
          }
)

