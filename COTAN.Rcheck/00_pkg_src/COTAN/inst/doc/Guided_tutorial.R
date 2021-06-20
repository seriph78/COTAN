## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 7
)
#options(rmarkdown.html_vignette.check_title = F)

## ----setup--------------------------------------------------------------------
library(COTAN)
library(data.table)
library(Matrix)
library(ggrepel)
library(factoextra)
library(Rtsne)
library(utils)
library(plotly)
library(tidyverse)
library(htmlwidgets)
library(MASS)
library(dendextend)

## -----------------------------------------------------------------------------
mycolours <- c("A" = "#8491B4B2","B"="#E64B35FF")
my_theme = theme(axis.text.x = element_text(size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF" ),
                 axis.text.y = element_text( size = 14, angle = 0, hjust = 0, vjust = .5, face = "plain", colour ="#3C5488FF"),
                 axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain", colour ="#3C5488FF"),
                 axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF"))
data_dir = tempdir()

## ----eval=FALSE, include=FALSE------------------------------------------------
#  Download the dataset for mouse cortex E17.5.
#  if (!file.exists(paste("inst/extdata/","E175_only_cortical_cells.txt.gz", sep = "/"))) {
#    download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/
#                  samples/GSM2861nnn/GSM2861514/suppl/GSM2861514_E175_Only_Cortical_Cells_DGE.txt.gz",
#                paste(data_dir,"E175_only_cortical_cells.txt.gz", sep = "/"),method = "wget", quiet = FALSE)
#  
#  }

## -----------------------------------------------------------------------------
#load(system.file("data", "ERCCraw.rda", package = "COTAN", mustWork = TRUE))
data("ERCCraw", package = "COTAN")

ERCCraw = as.data.frame(ERCCraw)
rownames(ERCCraw) = ERCCraw$V1
ERCCraw = ERCCraw[,2:ncol(ERCCraw)]
ERCCraw[1:5,1:5]

## -----------------------------------------------------------------------------
#out_dir = out_dir
out_dir = tempdir()

## -----------------------------------------------------------------------------
obj = new("scCOTAN",raw = ERCCraw)
#obj = initRaw(obj,GEO="GSM2861514" ,sc.method="Drop_seq",cond = "mouse cortex E17.5")
obj = initRaw(obj,GEO="ERCC" ,sc.method="10X",cond = "negative ERCC dataset")

## -----------------------------------------------------------------------------
genes_to_rem = rownames(obj@raw[grep('^mt', rownames(obj@raw)),]) 
obj@raw = obj@raw[!rownames(obj@raw) %in% genes_to_rem,]
cells_to_rem = colnames(obj@raw[which(colSums(obj@raw) == 0)])
obj@raw = obj@raw[,!colnames(obj@raw) %in% cells_to_rem]

## -----------------------------------------------------------------------------
#t = "E17.5_cortex"
t = "ERCC"

print(paste("Condition ",t,sep = ""))
#--------------------------------------
n_cells = length(colnames(obj@raw))
print(paste("n cells", n_cells, sep = " "))

n_it = 1

## -----------------------------------------------------------------------------
if(!file.exists(out_dir)){
  dir.create(file.path(out_dir))
}

if(!file.exists(paste(out_dir,"cleaning", sep = ""))){   
  dir.create(file.path(out_dir, "cleaning"))
}

## -----------------------------------------------------------------------------
ttm = clean(obj)

obj = ttm$object

ttm$pca.cell.2

## ----eval=F, include=T--------------------------------------------------------
#  pdf(file.path(out_dir,"cleaning",paste(t,"_",n_it,"_plots_before_cells_exlusion.pdf", sep = "")))
#  ttm$pca.cell.2
#  ggplot(ttm$D, aes(x=n,y=means)) + geom_point() +
#    geom_text_repel(data=subset(ttm$D, n > (max(ttm$D$n)- 15) ), aes(n,means,label=rownames(ttm$D[ttm$D$n > (max(ttm$D$n)- 15),])),
#                    nudge_y      = 0.05,
#                    nudge_x      = 0.05,
#                    direction    = "x",
#                    angle        = 90,
#                    vjust        = 0,
#                    segment.size = 0.2)+
#    ggtitle("B cell group genes mean expression")+my_theme +
#    theme(plot.title = element_text(color = "#3C5488FF", size = 20, face = "italic",vjust = - 5,hjust = 0.02 ))
#  dev.off()
#  
#  if (length(ttm$cl1) < length(ttm$cl2)) {
#    to_rem = ttm$cl1
#  }else{
#    to_rem = ttm$cl2
#  }
#  n_it = n_it+1
#  obj@raw = obj@raw[,!colnames(obj@raw) %in% to_rem]
#  #obj@raw = obj@raw[rownames(obj@raw) %in% rownames(cells),colnames(obj@raw) %in% colnames(cells)]
#  gc()
#  
#  ttm = clean(obj)
#  #ttm = clean.sqrt(obj, cells)
#  obj = ttm$object
#  
#  ttm$pca.cell.2
#  

## -----------------------------------------------------------------------------
pdf(file.path(out_dir,"cleaning",paste(t,"_",n_it,"_plots_before_cells_exlusion.pdf", sep = "")))
ttm$pca.cell.2
ggplot(ttm$D, aes(x=n,y=means)) + geom_point() +
  geom_text_repel(data=subset(ttm$D, n > (max(ttm$D$n)- 15) ), aes(n,means,label=rownames(ttm$D[ttm$D$n > (max(ttm$D$n)- 15),])),
                  nudge_y      = 0.05,
                  nudge_x      = 0.05,
                  direction    = "x",
                  angle        = 90,
                  vjust        = 0,
                  segment.size = 0.2)+
  ggtitle(label = "B cell group genes mean expression", subtitle = " - B group NOT removed -")+my_theme +
  theme(plot.title = element_text(color = "#3C5488FF", size = 20, face = "italic",vjust = - 10,hjust = 0.02 ),
        plot.subtitle = element_text(color = "darkred",vjust = - 15,hjust = 0.01 ))

dev.off()

## -----------------------------------------------------------------------------
nu_est = round(obj@nu, digits = 7)

plot.nu <-ggplot(ttm$pca_cells,aes(x=PC1,y=PC2, colour = log(nu_est)))

plot.nu = plot.nu + geom_point(size = 1,alpha= 0.8)+
  scale_color_gradient2(low = "#E64B35B2",mid =  "#4DBBD5B2", high =  "#3C5488B2" ,
                        midpoint = log(mean(nu_est)),name = "ln(nu)")+
  ggtitle("Cells PCA coloured by cells efficiency") +
  my_theme +  theme(plot.title = element_text(color = "#3C5488FF", size = 20),
                    legend.title=element_text(color = "#3C5488FF", size = 14,face = "italic"),
                    legend.text = element_text(color = "#3C5488FF", size = 11),
                    legend.key.width = unit(2, "mm"),
                    legend.position="right")

pdf(file.path(out_dir,"cleaning",paste(t,"_plots_PCA_efficiency_colored.pdf", sep = "")))
plot.nu
dev.off()

plot.nu

## -----------------------------------------------------------------------------
nu_df = data.frame("nu"= sort(obj@nu), "n"=c(1:length(obj@nu)))

ggplot(nu_df, aes(x = n, y=nu)) + 
    geom_point(colour = "#8491B4B2", size=1)+
    my_theme #+ ylim(0,1) + xlim(0,70)


## -----------------------------------------------------------------------------
yset = 0.4#threshold to remove low UDE cells
plot.ude <- ggplot(nu_df, aes(x = n, y=nu)) + 
    geom_point(colour = "#8491B4B2", size=1) + 
    my_theme + ylim(0,1) + xlim(0,400) +
    geom_hline(yintercept=yset, linetype="dashed", color = "darkred") +
    annotate(geom="text", x=200, y=0.25, 
             label=paste("to remove cells with nu < ",yset,sep = " "), 
             color="darkred", size=4.5)


pdf(file.path(out_dir,"cleaning",paste(t,"_plots_efficiency.pdf", sep = "")))
plot.ude
dev.off()

plot.ude

## -----------------------------------------------------------------------------
obj@meta[(nrow(obj@meta)+1),1:2] = c("Threshold low UDE cells:",yset)

to_rem = rownames(nu_df[which(nu_df$nu < yset),])
obj@raw = obj@raw[, !colnames(obj@raw) %in% to_rem]

## -----------------------------------------------------------------------------
ttm = clean(obj)
obj = ttm$object
ttm$pca.cell.2

## -----------------------------------------------------------------------------
pdf(file.path(out_dir,"cleaning",paste(t,"_plots_efficiency.pdf", sep = "")))
ggplot(nu_df, aes(x = n, y=nu)) + geom_point(colour = "#8491B4B2", size=1) +my_theme + #xlim(0,100)+
  annotate(geom="text", x=50, y=0.25, label="nothing to remove ", color="darkred")
dev.off()

## -----------------------------------------------------------------------------
nu_est = round(obj@nu, digits = 7)
plot.nu <-ggplot(ttm$pca_cells,aes(x=PC1,y=PC2, colour = log(nu_est)))
plot.nu = plot.nu + geom_point(size = 2,alpha= 0.8)+
  scale_color_gradient2(low = "#E64B35B2",mid =  "#4DBBD5B2", high =  "#3C5488B2" ,
                        midpoint = log(mean(nu_est)),name = "ln(nu)")+
  ggtitle("Cells PCA coloured by cells efficiency: last") +
  my_theme +  theme(plot.title = element_text(color = "#3C5488FF", size = 20),
                    legend.title=element_text(color = "#3C5488FF", size = 14,face = "italic"),
                    legend.text = element_text(color = "#3C5488FF", size = 11),
                    legend.key.width = unit(2, "mm"),
                    legend.position="right")

pdf(file.path(out_dir,"cleaning",paste(t,"_plots_PCA_efficiency_colored_FINAL.pdf", sep = "")))
plot.nu
dev.off()

plot.nu

## -----------------------------------------------------------------------------
# obj@yes_yes = c()
obj = cotan_analysis(obj,cores = 2)

## -----------------------------------------------------------------------------
obj = get.coex(obj)

# saving the structure 
saveRDS(obj,file = paste(out_dir,t,".cotan.RDS", sep = ""))

## -----------------------------------------------------------------------------
plot_GDI(obj, cond = "ERCC")

## ----echo=T-------------------------------------------------------------------
quant.p = get.GDI(obj)

head(quant.p)

## ----echo=T-------------------------------------------------------------------
AA=c("ERCC-00012","ERCC-00013","ERCC-00014") 
BB=c("ERCC-00016","ERCC-00017","ERCC-00019")
CC=c("ERCC-00022","ERCC-00024","ERCC-00028")

text.size = 12

quant.p$highlight = with(quant.p, ifelse(rownames(quant.p) %in% AA, "AA",
                                                 ifelse(rownames(quant.p) %in% CC,"Constitutive" ,
                                                               ifelse(rownames(quant.p) %in% BB,"BB" , "normal"))))

textdf <- quant.p[rownames(quant.p) %in% c(AA,CC,BB), ]
mycolours <- c("Con" = "#00A087FF","AA"="#E64B35FF","BB"="#F39B7FFF","normal" = "#8491B4B2")
f1 = ggplot(subset(quant.p,highlight == "normal" ), aes(x=sum.raw.norm, y=GDI)) +  geom_point(alpha = 0.1, color = "#8491B4B2", size=2.5)
GDI_plot = f1 +  geom_point(data = subset(quant.p,highlight != "normal"  ), aes(x=sum.raw.norm, y=GDI, colour=highlight),size=2.5,alpha = 0.8) +
  geom_hline(yintercept=quantile(quant.p$GDI)[4], linetype="dashed", color = "darkblue") +
  geom_hline(yintercept=quantile(quant.p$GDI)[3], linetype="dashed", color = "darkblue") +
  geom_hline(yintercept=1.5, linetype="dotted", color = "red", size= 0.5) +
  scale_color_manual("Status", values = mycolours)  +
  scale_fill_manual("Status", values = mycolours)  +
  xlab("log normalized counts")+ylab("GDI")+
  geom_label_repel(data =textdf , aes(x=sum.raw.norm, y=GDI, label = rownames(textdf),fill=highlight),
                   label.size = NA, 
                   alpha = 0.5, 
                   direction ="both",
                   na.rm=TRUE,
                   seed = 1234) +
  geom_label_repel(data =textdf , aes(x=sum.raw.norm, y=GDI, label = rownames(textdf)),
                   label.size = NA, 
                   segment.color = 'black',
                   segment.size = 0.5,
                   direction = "both",
                   alpha = 0.8, 
                   na.rm=TRUE,
                   fill = NA,
                   seed = 1234) +
  theme(axis.text.x = element_text(size = text.size, angle = 0, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF" ),
        axis.text.y = element_text( size = text.size, angle = 0, hjust = 0, vjust = .5, face = "plain", colour ="#3C5488FF"),  
        axis.title.x = element_text( size = text.size, angle = 0, hjust = .5, vjust = 0, face = "plain", colour ="#3C5488FF"),
        axis.title.y = element_text( size = text.size, angle = 90, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF"),
        legend.title = element_blank(),
        legend.text = element_text(color = "#3C5488FF",face ="italic" ),
        legend.position = "bottom")  # titl)
legend <- cowplot::get_legend(GDI_plot)
GDI_plot =GDI_plot + theme(
        legend.position = "none") 
GDI_plot


## -----------------------------------------------------------------------------
list.genes = list("Ref.col"= BB, "AA"=AA, "Const."=CC )

plot_heatmap(df_genes = list.genes,sets = c(1:3),conditions = "ERCC",dir = out_dir)

## -----------------------------------------------------------------------------
plot_general.heatmap(prim.markers = c("ERCC-00014","ERCC-00019"), condition = "ERCC",dir = out_dir, p_value = 0.05,)

## -----------------------------------------------------------------------------
get.observed.ct(object = obj, g1 = "ERCC-00014",g2 = "ERCC-00019")

## -----------------------------------------------------------------------------
get.expected.ct(object = obj, g1 = "ERCC-00014",g2 = "ERCC-00019")

## -----------------------------------------------------------------------------
if (file.exists(paste(out_dir,t,".cotan.RDS", sep = ""))) {
  #Delete file if it exists
  file.remove(paste(out_dir,t,".cotan.RDS", sep = ""))
}
unlink(paste(out_dir,"cleaning", sep = ""),recursive = TRUE)

print(sessionInfo())

