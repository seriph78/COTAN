pkgname <- "COTAN"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('COTAN')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("automatic.COTAN.object.creation")
### * automatic.COTAN.object.creation

flush(stderr()); flush(stdout())

### Name: automatic.COTAN.object.creation
### Title: automatic.COTAN.object.creation
### Aliases: automatic.COTAN.object.creation
###   automatic.COTAN.object.creation,data.frame-method

### ** Examples


data("raw.dataset")
obj = automatic.COTAN.object.creation(df= raw.dataset,
out_dir =  tempdir(),
GEO = "test_GEO",
sc.method = "test_method",
cond = "test")




cleanEx()
nameEx("clean")
### * clean

flush(stderr()); flush(stdout())

### Name: clean
### Title: clean
### Aliases: clean clean,scCOTAN-method

### ** Examples

data("ERCC.cotan")
ttm = clean(ERCC.cotan)




cleanEx()
nameEx("cotan_analysis")
### * cotan_analysis

flush(stderr()); flush(stdout())

### Name: cotan_analysis
### Title: cotan_analysis
### Aliases: cotan_analysis cotan_analysis,scCOTAN-method

### ** Examples

data("ERCC.cotan")
ERCC.cotan = cotan_analysis(ERCC.cotan)




cleanEx()
nameEx("get.GDI")
### * get.GDI

flush(stderr()); flush(stdout())

### Name: get.GDI
### Title: get.GDI
### Aliases: get.GDI get.GDI,scCOTAN-method

### ** Examples

data("ERCC.cotan")
quant.p = get.GDI(ERCC.cotan)



cleanEx()
nameEx("get.coex")
### * get.coex

flush(stderr()); flush(stdout())

### Name: get.coex
### Title: get.coex
### Aliases: get.coex get.coex,scCOTAN-method

### ** Examples


data("ERCC.cotan")
obj = get.coex(ERCC.cotan)




cleanEx()
nameEx("get.expected.ct")
### * get.expected.ct

flush(stderr()); flush(stdout())

### Name: get.expected.ct
### Title: get.expected.ct
### Aliases: get.expected.ct get.expected.ct,scCOTAN-method

### ** Examples

data("ERCC.cotan")
g1 = rownames(ERCC.cotan@raw)[sample(nrow(ERCC.cotan@raw), 1)]
g2 = rownames(ERCC.cotan@raw)[sample(nrow(ERCC.cotan@raw), 1)]
while (g1 %in% ERCC.cotan@hk) {
g1 = rownames(ERCC.cotan@raw)[sample(nrow(ERCC.cotan@raw), 1)]
}

while (g2 %in% ERCC.cotan@hk) {
   g2 = rownames(ERCC.cotan@raw)[sample(nrow(ERCC.cotan@raw), 1)]
}
get.expected.ct(object = ERCC.cotan, g1 = g1, g2 = g2)



cleanEx()
nameEx("get.gene.coexpression.space")
### * get.gene.coexpression.space

flush(stderr()); flush(stdout())

### Name: get.gene.coexpression.space
### Title: get.gene.coexpression.space
### Aliases: get.gene.coexpression.space
###   get.gene.coexpression.space,scCOTAN-method

### ** Examples

data("ERCC.cotan")
df = get.gene.coexpression.space(ERCC.cotan, n.genes.for.marker = 10,
primary.markers=rownames(ERCC.cotan@raw[sample(nrow(ERCC.cotan@raw),5),]))



cleanEx()
nameEx("get.observed.ct")
### * get.observed.ct

flush(stderr()); flush(stdout())

### Name: get.observed.ct
### Title: get.observed.ct
### Aliases: get.observed.ct get.observed.ct,scCOTAN-method

### ** Examples

data("ERCC.cotan")
g1 = rownames(ERCC.cotan@raw)[sample(nrow(ERCC.cotan@raw), 1)]
g2 = rownames(ERCC.cotan@raw)[sample(nrow(ERCC.cotan@raw), 1)]
get.observed.ct(object = ERCC.cotan, g1 = g1, g2 = g2)



cleanEx()
nameEx("get.pval")
### * get.pval

flush(stderr()); flush(stdout())

### Name: get.pval
### Title: get.pval
### Aliases: get.pval get.pval,scCOTAN-method

### ** Examples


data("ERCC.cotan")
ERCC.cotan = get.pval(ERCC.cotan,type_stat="S")




cleanEx()
nameEx("initRaw")
### * initRaw

flush(stderr()); flush(stdout())

### Name: initRaw
### Title: initRaw
### Aliases: initRaw initRaw,scCOTAN-method

### ** Examples


data("raw.dataset")
obj = new("scCOTAN", raw = raw.dataset)
obj = initRaw(obj, GEO="code" , sc.method="10X",cond = "mouse dataset")





cleanEx()
nameEx("plot_GDI")
### * plot_GDI

flush(stderr()); flush(stdout())

### Name: plot_GDI
### Title: plot_GDI
### Aliases: plot_GDI plot_GDI,scCOTAN-method

### ** Examples

data("ERCC.cotan")
plot_GDI(ERCC.cotan, cond = "ERCC")



cleanEx()
nameEx("plot_general.heatmap")
### * plot_general.heatmap

flush(stderr()); flush(stdout())

### Name: plot_general.heatmap
### Title: plot_general.heatmap
### Aliases: plot_general.heatmap plot_general.heatmap,ANY-method

### ** Examples

## Not run: 
##D plot_general.heatmap(dir=input_dir,
##D condition = "E17.5",
##D prim.markers  = c("Mef2c","Mef2a","Mef2d"),
##D symmetric = FALSE,
##D markers.list = c("Reln","Satb2","Cux1","Bcl11b","Tbr1","Sox5","Foxp2","Slc17a6","Slc17a7"),
##D p_value = 0.05)
## End(Not run)



cleanEx()
nameEx("plot_heatmap")
### * plot_heatmap

flush(stderr()); flush(stdout())

### Name: plot_heatmap
### Title: plot_heatmap
### Aliases: plot_heatmap plot_heatmap,ANY-method

### ** Examples

## Not run: 
##D # some genes
##D primary.markers = c("Tbr1","Tubb3","Neurod1", "Stmn1","Notch1","Vim","Sox2","Pax6","Hes5")
##D a example of named list of different gene set
##D gene.sets.list = list("primary.markers"=primary.markers,
##D                    "2.Radial Glia" = c("Vim","Sox2","Pax6","Hes5","Hes1","Fabp7"),
##D                    "PNGs"=c("Map2","Tubb3","Neurod1","Nefm","Nefl","Dcx","Tbr1"),
##D                    "constitutive" = c("Calm1","Cox6b1","Ppia","Rpl18","Cox7c",
##D                                       "Erh","H3f3a","Taf1b","Taf2",
##D                                       "Gapdh","Actb", "Golph3", "Mtmr12",
##D                                       "Zfr", "Sub1", "Tars", "Amacr"),
##D                    "4.Mat.neu."= c("Map2","Rbfox3","Nefl","Nefh","Nefm","Mapt"))
##D plot_heatmap(p_v = 0.05, df_genes =gene.sets.list ,
##D sets =c(2,3,4,6) ,conditions =c("E11.5","E13.5","E14.5") ,dir = input_dir)
## End(Not run)



cleanEx()
nameEx("scCOTAN-class")
### * scCOTAN-class

flush(stderr()); flush(stdout())

### Name: scCOTAN-class
### Title: scCOTAN-class
### Aliases: scCOTAN-class scCOTAN

### ** Examples


data("ERCCraw")
obj = new("scCOTAN",raw = data)





### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
