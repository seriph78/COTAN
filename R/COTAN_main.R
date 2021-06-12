## Main functions

#' scCOTAN-class
#'
#' Define my COTAN structure
#' @slot raw ANY. To store the raw data matrix
#' @slot raw.norm ANY. To store the raw data matrix divided for the cell
#' efficiency estimated (nu)
#' @slot coex ANY. The coex matrix (sparce)
#' @slot nu vector.
#' @slot lambda vector.
#' @slot a vector.
#' @slot hk vector.
#' @slot n_cells numeric.
#' @slot meta data.frame.
#' @slot yes_yes ANY.
#' @slot clusters vector.
#' @slot cluster_data data.frame.
#'
#' @return the object class
#' @export
#' @examples
#'
#' data("ERCCraw")
#' obj = new("scCOTAN",raw = data)
#'
#'
setClass("scCOTAN", slots =
                        c(raw="ANY",raw.norm="ANY", coex="ANY",
                            nu="vector",lambda="vector",a="vector",
                            hk = "vector", n_cells = "numeric",
                            meta="data.frame",yes_yes="ANY", clusters="vector",
                            cluster_data="data.frame")) -> scCOTAN

#' initRaw
#'
#' It starts to fill some fields of cotan object.
#' @param object the dataframe containing the raw data: it should be a data.frame, not a matrix.
#' @param GEO a code reporting the GEO identification or other specific dataset
#' code
#' @param sc.method a string reporting the method used for the sequencing
#' @param cond a string reporting the specific sample condition or time point
#'
#' @return A COTAN object to store all information
#' @import methods
#' @export
#' @rdname initRaw
#' @examples
#'
#' data("raw.dataset")
#' obj = new("scCOTAN", raw = raw)
#' obj = initRaw(obj, GEO="code" , sc.method="10X",cond = "mouse dataset")
#'
#'
setGeneric("initRaw", function(object,GEO,sc.method="10X", cond)
    standardGeneric("initRaw"))
#' @rdname initRaw
setMethod("initRaw","scCOTAN",
            function(object,GEO,sc.method,cond) {
                print("Initializing S4 object")
                if( (!attr(object@raw, "class")[1] == "dgCMatrix") |
                   is.null(attr(object@raw, "class")) ){
                    object@raw = methods::as(as.matrix(object@raw),
                                             "sparseMatrix")
                }
                if(! all((object@raw - round(object@raw)) == 0)){
                    print("WARNING! Input data contains not integer numbers!")

                }
                object@meta[1,1:2] = c("GEO:",GEO)
                object@meta[2,1:2] = c("scRNAseq method:",sc.method)
                object@meta[3,1] = "starting n. of cells:"
                object@meta[3,2] = ncol(object@raw)
                object@meta[4,1:2] = c("Condition sample:",cond)

              #object@clusters = rep(NA,ncol(object@raw))
              #names(object@clusters)=colnames(object@raw)
                return(object)
          }
)


#' clean
#'
#' Main function that can be used to check and clean the dataset. It also
#' produce
#' (using the function fun_linear) and store
#' the estimators for nu and lambda. It also fill the raw.norm (raw / nu) and
#' n_cell
#' (the initial number of cells in the dataset)
#' @param object COTAN object
#' @return a list of objects containing: "cl1" is the first cell cluster, "cl2"
#' is the
#' second cell cluster,
#' "pca.cell.2" is a ggplot2 cell pca plot, "object" is the COTAN object with
#' saved the
#' estimated lambda and mu, "mu_estimator", "D"
#' "pca_cells" pca numeric data.
#' @export
#'
#' @import ggplot2
#' dplyr
#'
#' @importFrom tibble rownames_to_column
#' @importFrom stats hclust
#' @importFrom stats cutree
#'
#' @importFrom Matrix rowMeans
#' @importFrom utils head
#' @importFrom Matrix colMeans
#' @rdname clean
#' @examples
#' data("ERCC.cotan")
#' ttm = clean(ERCC.cotan)
#'
setGeneric("clean", function(object) standardGeneric("clean"))
#' @rdname clean
setMethod("clean","scCOTAN",

              function(object) {

                 mycolours <- c("A" = "#8491B4B2","B"="#E64B35FF")

                 my_theme = theme(axis.text.x = element_text(size = 14,
                                                             angle = 0, hjust = .5,
                                                             vjust = .5,
                                                              face = "plain",
                                                             colour ="#3C5488FF" ),
                                   axis.text.y = element_text( size = 14, angle = 0, hjust = 0,
                                                               vjust = .5,
                                                               face = "plain",
                                                               colour ="#3C5488FF"),
                                   axis.title.x = element_text( size = 14,
                                                                angle = 0, hjust = .5,
                                                                vjust = 0,
                                                                face = "plain",
                                                                colour ="#3C5488FF"),
                                   axis.title.y = element_text( size = 14, angle = 90, hjust = .5,
                                                                vjust = .5,
                                                                face = "plain",
                                                                colour ="#3C5488FF"))


                cells = get.zero_one.cells(object)

                object@raw = object@raw[rownames(cells), colnames(cells)]

                list1 = fun_linear(object)

                gc()

                dist_cells = list1$dist_cells
                pca_cells = list1$pca_cells
                t_to_clust = as.matrix(list1$t_to_clust)
                mu_estimator = list1$mu_estimator
                object = list1$object
                to_clust = list1$to_clust


                raw_norm <- t(t(as.matrix(object@raw)) *
                                  (1/(as.vector(object@nu))))
                #raw_norm <- as.matrix(object@raw) %b/% matrix(as.vector(object@nu), nrow = 1)
                raw_norm <- as(as.matrix(raw_norm), "sparseMatrix")

                object@raw.norm = raw_norm


                hc_cells = hclust(dist_cells, method = "complete")

                groups <- cutree(hc_cells, k=2)

                if(length(groups[groups == 1]) < length(groups[groups == 2])  ){
                  groups[groups == 1] <- "B"
                  groups[groups == 2] <- "A"

                }else{
                  groups[groups == 1] <- "A"
                  groups[groups == 2] <- "B"

                }


                cl2 = names(which(cutree(hc_cells, k = 2) == 2))
                cl1 = names(which(cutree(hc_cells, k = 2) == 1))

                if (length(cl2) > length(cl1) ) {
                  cl2 = names(which(cutree(hc_cells, k = 2) == 1))
                  cl1 = names(which(cutree(hc_cells, k = 2) == 2))
                }


                t_to_clust = cbind(as.data.frame(t_to_clust),groups)

              # ---- next: to check which genes are specific for the B group of cells
                to_clust = as.matrix(to_clust)
                B = as.data.frame(to_clust[,colnames(to_clust) %in% cl2])
                colnames(B)=cl2
                B = rownames_to_column(B)
                if (dim(B)[2]>2) {
                  B = B[order(rowMeans(B[,2:length(colnames(B))]),
                              decreasing = TRUE), ]
                }else{
                  B = B[order(B[,2],decreasing = TRUE), ]
                }

                print(utils::head(B, 15))

                C = arrange(B,rowMeans(B[2:length(colnames(B))]))
                rownames(C) = C$rowname
                D = data.frame("means"=rowMeans(C[2:length(colnames(C))]),
                               "n"=NA )
                D = D[D$means>0,]
                D$n = c(1:length(D$means))

                # check if the pca plot is clean enought and from the printed genes,
                #if the smalest group of cells are caratterised by particular genes

                pca_cells = cbind(pca_cells,"groups"=t_to_clust$groups)

                pca.cell.1 = ggplot(subset(pca_cells,groups == "A" ),
                                    aes(x=PC1, y=PC2,colour =groups)) +
                  geom_point(alpha = 0.5, size=3)

                pca.cell.2=  pca.cell.1 + geom_point(data = subset(pca_cells,
                            groups != "A" ),aes(x=PC1, y=PC2,colour =groups),
                                                   alpha = 0.8, size=3)+
                  scale_color_manual("groups", values = mycolours)  +
                  my_theme + theme(legend.title = element_blank(),
                                   legend.text = element_text( size = 12,
                                                            color = "#3C5488FF",
                                                            face ="italic" ),
                                   legend.position="bottom")

                object@n_cells = length(colnames(object@raw))

                output = list("cl1"=cl1,"cl2"=cl2,"pca.cell.2"=pca.cell.2,
                              "object"=object,
                              "mu_estimator"=mu_estimator, "D"=D,
                              "pca_cells"=pca_cells)
                return(output)
          }
)


#' cotan_analysis
#'
#' This is the main function that estimates the a vector to store all the
#' negative binomial
#' dispersion factors. It need to be run after \code{\link{clean}}
#' @param object A COTAN object
#' @param cores number of cores to use. Default is 11.
#' @importFrom parallel mclapply
#' @return a COTAN object
#' @export
#' @rdname cotan_analysis
#' @examples
#' data("ERCC.cotan")
#' ERCC.cotan = cotan_analysis(ERCC.cotan)
#'
setGeneric("cotan_analysis", function(object, cores= 1)
    standardGeneric("cotan_analysis"))
#' @rdname cotan_analysis
setMethod("cotan_analysis","scCOTAN",
          function(object, cores= 1) {

              print("cotan analysis")

              cells=as.matrix(object@raw)
              #---------------------------------------------------
              # Cells matrix : formed by row data matrix changed to 0-1 matrix
              cells[cells > 0] <- 1
              cells[cells <= 0] <- 0

              # exlude the effective ubiqutarius genes and saved in a separate file
              mu_estimator = mu_est(object)

              object = hk_genes(object)

              hk =object@hk

              mu_estimator = mu_estimator[!rownames(mu_estimator) %in% hk,]
              cells = cells[!rownames(cells) %in% hk, ]

              mu_estimator = as.matrix(mu_estimator)
              print("start a minimization")
              gc()
              p=1
              tot = list()


              while(p <= length(rownames(mu_estimator))) {
                  #bb=Sys.time()
                  if((p+200) <= length(rownames(mu_estimator))){
                      tot1 =  parallel::mclapply(rownames(mu_estimator)[p:(p+200)],
                                       fun_my_opt, ce.mat= cells, mu_est= mu_estimator,
                                       mc.cores = cores)

                  }else{
                      print("Final trance!")
                      tot1 = parallel::mclapply(rownames(mu_estimator)[p:length(rownames(mu_estimator))],
                                      fun_my_opt, ce.mat=cells, mu_est= mu_estimator, mc.cores = cores)
                  }
                  tot = append(tot, tot1)
                  p=p+200+1
                  if((p %% 10)==0){
                      #print(Sys.time()-bb)
                      print(paste("Next gene:",rownames(mu_estimator)[p],"number",p, sep = " "))
                  }
              }
              gc()

              tot2 = tot[[1]]
              for (tt in 2:length(tot)) {
                  tot2= rbind(tot2,tot[[tt]])

              }

              object@a = tot2$a
              names(object@a) = rownames(tot2)
              #save(tot2 , file = paste(out_dir,"a_minimization_",t, sep = ""))
              print(paste("a min:",min(tot2$a) ,"| a max",max(tot2$a) , "| negative a %:",
                          sum(tot2$a <0)/nrow(tot2)*100 ,sep=" "))

              gc()

              return(object)
          }
)


#' get.coex
#'
#' This function estimates and stores the coex*sqrt(n_cells) matrix in the coex field.
#' It need to be run after \code{\link{cotan_analysis}}
#' @param object A COTAN object
#'
#' @return It returns a COTAN object
#' @export
#' @rdname get.coex
#' @examples
#'
#' data("ERCC.cotan")
#' obj = get.coex(ERCC.cotan)
#'
setGeneric("get.coex", function(object) standardGeneric("get.coex"))
#' @rdname get.coex
setMethod("get.coex","scCOTAN",
          function(object) {
              print("coex dataframe creation")
              hk = object@hk
              ll = obs_ct(object)

              object = ll$object

              ll$no_yes= ll$no_yes[!rownames(ll$no_yes) %in% hk,!colnames(ll$no_yes) %in% hk]
              ll$no_no = ll$no_no[!rownames(ll$no_no) %in% hk,!colnames(ll$no_no) %in% hk]
              si_si = object@yes_yes[!rownames(object@yes_yes) %in% hk,!colnames(object@yes_yes)
                                     %in% hk]
              ll$yes_no = ll$yes_no[!rownames(ll$yes_no) %in% hk,!colnames(ll$yes_no) %in% hk]

              est = expected_ct(object)

              new_estimator_si_si = as.matrix(est$estimator_yes_yes)
              new_estimator_si_si[new_estimator_si_si < 1] <- 1
              new_estimator_si_no = as.matrix(est$estimator_yes_no)
              new_estimator_si_no[new_estimator_si_no < 1] <- 1
              new_estimator_no_no = as.matrix(est$estimator_no_no)
              new_estimator_no_no[new_estimator_no_no < 1] <- 1
              new_estimator_no_si = as.matrix(est$estimator_no_yes)
              new_estimator_no_si[new_estimator_no_si < 1] <- 1
              print("coex estimation")
              coex = ((as.matrix(si_si) - as.matrix(est$estimator_yes_yes))/new_estimator_si_si) +
                  ((as.matrix(ll$no_no) - as.matrix(est$estimator_no_no))/new_estimator_no_no) -
                  ((as.matrix(ll$yes_no) - as.matrix(est$estimator_yes_no))/new_estimator_si_no) -
                  ((as.matrix(ll$no_yes) - as.matrix(est$estimator_no_yes))/new_estimator_no_si)

              coex = coex / sqrt(1/new_estimator_si_si + 1/new_estimator_no_no + 1/
                                     new_estimator_si_no + 1/new_estimator_no_si)

              #print("coex low values substituted with 0")
              #coex[abs(coex) <= 1]=0
              print("Diagonal coex values substituted with 0")
              diag(coex) = 0
              coex = coex / sqrt(object@n_cells)
              object@coex = spMat(coex)

              return(object)
          }
)


#' get.pval
#'
#' This function computes the p-values for genes in the COTAN object. It can be used genome-wide
#' or setting some specific genes of interest. By default it computes the p-values using the S
#' statistics (\eqn{\chi^{2}})
#' @param object a COTAN object
#' @param gene.set.col an array of genes. It will be put in columns.
#' If left empty the function will do it genome-wide.
#' @param gene.set.row an array of genes. It will be put in rows.
#' If left empty the function will do it genome-wide.
#' @param type_stat By default it computes the S (\eqn{\chi^{2}})
#'
#' @return a p-value matrix
#' @export
#' @rdname get.pval
#' @importFrom Matrix forceSymmetric
#' @examples
#'
#' data("ERCC.cotan")
#' ERCC.cotan = get.pval(ERCC.cotan,type_stat="S")
#'
setGeneric("get.pval", function(object, gene.set.col=c(),gene.set.row=c(), type_stat="S" )
    standardGeneric("get.pval"))
#' @rdname get.pval
setMethod("get.pval","scCOTAN",
          function(object, gene.set.col=c(),gene.set.row=c(), type_stat="S") {
              object@coex = Matrix::forceSymmetric(object@coex, uplo="L" )
              print(gene.set.col)

              if (!is.null(gene.set.row)) {

                  # a set for rows, not Genome Wide
                  cond.row = "on a set of genes on rows"
                  if (is.null(gene.set.col)) {
                      #print("Error: can't have genome wide on columns and not rows! Use a
                      #subset on gene.set.col, not on rows.")
                      #break()
                      stop("Error: can't have genome wide on columns and not rows! Use a
                           subset on gene.set.col, not on rows.")
                  }
                  cond.col = "on a set of genes on columns"

              }else{
                  cond.row = "genome wide on rows"
                  if (is.null(gene.set.col)) {
                      cond.col = "genome wide on columns"
                  }else{
                      cond.col = "on a set of genes on columns"
                  }

              }

            print(paste("Get p-values", cond.col, cond.row, sep = " "))

            if (type_stat == "S") {
                print("Using function S")
                S = get.S(object)
            }else if(type_stat == "G"){
                print("Using function G")
                S = get.G(object)
            }else if(type_stat == "S2"){
                print("Using function G")
                S = get.S2(object)
            }


              if(cond.col == "on a set of genes on columns"){
                  S = S[,colnames(S) %in% gene.set.col]
              }
              if(cond.row == "on a set of genes on rows"){
                  S = S[rownames(S) %in% gene.set.row,]
              }


              p_value = pchisq(as.matrix(S), df=1, lower.tail=FALSE)


              return(p_value)
          }
)





#' plot_heatmap
#'
#' This is the function that create the heatmap of one or more COTAN object.
#'
#' @param p_val.tr p-value threshold. Default is 0.05
#' @param df_genes this is a list of gene array. The first array will define genes in the columns.
#' @param sets This is a numeric array indicating from which fields of the previous list
#' will be considered
#' @param conditions An array of prefixes indicating the different files.
#' @param dir The directory in which are all COTAN files (corresponding to the previous prefixes)
#'
#' @return a ggplot object
#' @importFrom Matrix forceSymmetric
#' @export
#'
#' @import ggplot2
#' @import tidyr
#' @import scales
#' @rdname plot_heatmap
#' @examples
#' \dontrun{
#' # some genes
#' primary.markers = c("Tbr1","Tubb3","Neurod1", "Stmn1","Notch1","Vim","Sox2","Pax6","Hes5")
#' a example of named list of different gene set
#' gene.sets.list = list("primary.markers"=primary.markers,
#'                    "2.Radial Glia" = c("Vim","Sox2","Pax6","Hes5","Hes1","Fabp7"),
#'                    "PNGs"=c("Map2","Tubb3","Neurod1","Nefm","Nefl","Dcx","Tbr1"),
#'                    "constitutive" = c("Calm1","Cox6b1","Ppia","Rpl18","Cox7c",
#'                                       "Erh","H3f3a","Taf1b","Taf2",
#'                                       "Gapdh","Actb", "Golph3", "Mtmr12",
#'                                       "Zfr", "Sub1", "Tars", "Amacr"),
#'                    "4.Mat.neu."= c("Map2","Rbfox3","Nefl","Nefh","Nefm","Mapt"))
#' plot_heatmap(p_v = 0.05, df_genes =gene.sets.list ,
#' sets =c(2,3,4,6) ,conditions =c("E11.5","E13.5","E14.5") ,dir = input_dir)
#' }
setGeneric("plot_heatmap", function(p_val.tr = 0.05, df_genes , sets, conditions, dir)
    standardGeneric("plot_heatmap"))
#' @rdname plot_heatmap
setMethod("plot_heatmap","ANY",
          function(p_val.tr = 0.05, df_genes , sets, conditions, dir) {
              print("plot heatmap")

              gr = df_genes[[1]]

              ge = unique(array(sort(unlist(df_genes[sets]))))

              df.to.print = data.frame()

              for(ET in conditions){
                  print(paste("Loading condition",ET,sep=" "))
                  obj = readRDS(paste(dir,ET,".cotan.RDS", sep = ""))
                  obj@coex = Matrix::forceSymmetric(obj@coex, uplo="L" )
                  if(any(gr %in% rownames(obj@coex)) == FALSE){
                      #print(paste("primary markers all absent in ", ET, sep = " "))
                      #break()
                      paste0("primary markers all absent in ", ET)
                      stop()

                  }
                  p_val = get.pval(obj,gene.set.col = gr, gene.set.row = ge)
                  p_val = as.data.frame(p_val)

                  #this to add some eventually effective housekeeping genes

                  if (any(ge %in% obj@hk)) {
                      genes.to.add = ge[ge %in% obj@hk]
                      temp.hk.rows = as.data.frame(matrix(ncol = ncol(p_val), nrow =
                                                              length(genes.to.add)))
                      rownames(temp.hk.rows) =genes.to.add
                      colnames(temp.hk.rows) = colnames(p_val)
                      temp.hk.rows = 1
                      p_val = rbind(p_val,temp.hk.rows)
                  }

                  if (any(gr %in% obj@hk)) {
                      genes.to.add = gr[gr %in% obj@hk]
                      temp.hk.cols = as.data.frame(matrix(ncol = length(genes.to.add), nrow =
                                                              nrow(p_val) ))
                      colnames(temp.hk.cols) =genes.to.add
                      rownames(temp.hk.cols) = rownames(p_val)
                      temp.hk.cols = 1
                      p_val = cbind(p_val,temp.hk.cols)
                  }

                  #------------------------

                  p_val$g2 = as.vector(rownames(p_val))
                  #df.temp.pval <- pivot_longer(p_val, cols=1:(ncol(p_val)-1), names_to = "g1",
                  #                             values_to = "p_val")
                  df.temp.pval <- pivot_longer(p_val, cols=seq_along(colnames(p_val))-1, names_to = "g1",
                                               values_to = "p_val")

                  coex = obj@coex[rownames(obj@coex) %in% ge,colnames(obj@coex) %in% gr]
                  #coex = coex / sqrt(obj@n_cells)
                  coex = as.data.frame(as.matrix(coex))
                  #this to add some eventually effective housekeeping genes
                  if (any(ge %in% obj@hk)) {
                      temp.hk.rows = 0
                      coex = rbind(coex,temp.hk.rows)
                  }

                  if (any(gr %in% obj@hk)) {
                      temp.hk.cols = 0
                      coex = cbind(coex,temp.hk.cols)
                  }
                  #---------------------------------------------------------

                  coex$g2 = as.vector(rownames(coex))

                  df.temp.coex <- pivot_longer(coex, cols=seq_along(colnames(p_val))-1, names_to = "g1",
                                               values_to = "coex")


                  #df.temp.coex <- pivot_longer(coex, cols=1:(ncol(p_val)-1), names_to = "g1",
                   #                            values_to = "coex")

                  df.temp = merge(df.temp.coex, df.temp.pval)
                  df.temp$time = ET
                  df.temp$type = NA
                  df.temp$absent = NA
                  df.temp2 = data.frame()
                  for (type in names(df_genes)[sets]) {

                      for (g1 in gr) {
                          #for (g in df_genes[[type]]) {
                          #if(all(is.na(df.temp[df.temp$g2 ==g, "type"]))){
                          # df.temp[df.temp$g2 ==g, "type"] = type
                          #}else{
                          #tt = df.temp[df.temp$g2 ==g, ]
                          tt = df.temp[df.temp$g2 %in% df_genes[[type]] & df.temp$g1 == g1, ]
                          #control if the subset is smaller than the number of wanted genes
                          if(dim(tt)[1] < length(df_genes[[type]])){
                              n.row = length(df_genes[[type]]) - dim(tt)[1]
                              t.rows = as.data.frame(matrix(nrow = n.row, ncol=7))
                              colnames(t.rows)=colnames(tt)
                              t.rows[,"g1"]=g1
                              t.rows[,"time"]= ET
                              t.rows[,"absent"]= "yes"
                              #t.rows[,"coex"]= 0
                              t.rows[,"p_val"]= 1
                              t.rows[,"g2"]=df_genes[[type]][!df_genes[[type]] %in% tt$g2]
                              tt = rbind(tt,t.rows)
                          }
                          tt$type = type
                          #df.temp = rbind(df.temp,tt)
                          df.temp2 =rbind(df.temp2,tt)
                          #}
                      }
                      print(type)

                  }
                  df.temp = df.temp2
                  df.temp$t_hk = ifelse((df.temp$g2  %in% obj@hk) | (df.temp$g1  %in% obj@hk),
                                        "hk", "n")

                  df.temp[df.temp$p_val > p_val.tr,]$coex = 0

                  df.to.print =  rbind(df.to.print,df.temp)

              }

              print(paste("min coex:",min(df.to.print$coex, na.rm = TRUE), "max coex",
                          max(df.to.print$coex, na.rm = TRUE),sep = " "))

              heatmap = ggplot(data = subset(df.to.print,type %in%  names(df_genes)[sets] ),
                               aes(time, factor(g2, levels = rev(levels(factor(g2)))))) +

                  geom_tile(aes(fill = coex),colour = "black", show.legend = TRUE) +
                  facet_grid( type ~ g1  ,scales = "free", space = "free") +
                  scale_fill_gradient2(low = "#E64B35FF", mid = "gray93",   high = "#3C5488FF",
                                       midpoint = 0,
                                       na.value = "grey80", space = "Lab", guide = "colourbar",
                                       aesthetics = "fill", oob=scales::squish)+ #, limits = lim_coex
                  theme(axis.title.x = element_blank(),
                        panel.spacing = unit(0, "lines"),
                        strip.background = element_rect(fill="#8491B44C"),
                        strip.text.y  = element_text(size = 9, colour = "#3C5488FF"),
                        strip.text.x  = element_text(size = 9,angle= 90 ,colour = "#3C5488FF"),
                        axis.title.y = element_blank(),
                        # axis.text.x = element_blank(),
                        axis.text.y = element_text( size = 9, angle = 0, hjust = 0, vjust = .5,
                                                    face = "plain", colour ="#3C5488FF"),
                        #legend.title = element_blank(),
                        legend.text = element_text(color = "#3C5488FF",face ="italic" ),
                        legend.position = "bottom",
                        legend.title=element_blank(),
                        legend.key.height = unit(2, "mm"))

              heatmap
              return(heatmap)
          }
)



#' plot_general.heatmap
#'
#' This function is used to plot an heatmap made using only some genes, as markers,
#' and collecting all other genes correlated
#' with these markers with a p-value smaller than the set threshold. Than all relations are plotted.
#' Primary markers will be plotted as groups of rows. Markers list will be plotted as columns.
#' @param prim.markers A set of genes plotted as rows.
#' @param markers.list A set of genes plotted as columns.
#' @param dir The directory where the COTAN object is stored.
#' @param condition The prefix for the COTAN object file.
#' @param p_value The p-value threshold
#' @param symmetric A boolean: default F. If T the union of prim.markers and marker.list
#' is sets as both rows and column genes
#'
#' @return A ggplot2 object
#' @import ComplexHeatmap
#' @importFrom grid gpar
#' @importFrom  stats quantile
#' @importFrom circlize colorRamp2
#' @importFrom Matrix forceSymmetric
#' @importFrom rlang is_empty
#' @export
#' @rdname plot_general.heatmap
#' @examples
#' \dontrun{
#' plot_general.heatmap(dir=input_dir,
#' condition = "E17.5",
#' prim.markers  = c("Mef2c","Mef2a","Mef2d"),
#' symmetric = FALSE,
#' markers.list = c("Reln","Satb2","Cux1","Bcl11b","Tbr1","Sox5","Foxp2","Slc17a6","Slc17a7"),
#' p_value = 0.05)
#' }
setGeneric("plot_general.heatmap", function(prim.markers =c("Satb2","Bcl11b","Cux1","Fezf2","Tbr1"),
                                            markers.list= c(),
                                            dir, condition,
                                            p_value = 0.001,
                                            symmetric = TRUE) standardGeneric("plot_general.heatmap"))
#' @rdname plot_general.heatmap
setMethod("plot_general.heatmap","ANY",
          function(prim.markers =c("Satb2","Bcl11b","Cux1","Fezf2","Tbr1"),markers.list=c(), dir,
                   condition,p_value = 0.001, symmetric = TRUE) {
              print("ploting a general heatmap")

              if(symmetric == TRUE){
                  markers.list = as.list(c( unlist(prim.markers),unlist(markers.list)))
              }

              if (is.null(markers.list)) {
                  markers.list = as.list(prim.markers)
              }else{
                  markers.list = as.list(markers.list)
              }

              obj = readRDS(paste(dir,condition,".cotan.RDS", sep = ""))
              obj@coex = Matrix::forceSymmetric(obj@coex, uplo="L" )
              coex = as.matrix(obj@coex)
              no_genes = unique(c(unlist(markers.list),prim.markers))[!unique(c(unlist(markers.list),
                                                                                prim.markers))
                                                                      %in% colnames(coex)]

              if(! rlang::is_empty(no_genes)){
                  print(paste(no_genes,"not present!",sep = " "))
              }

              coex = coex[,colnames(coex) %in% unique(c(unlist(markers.list),prim.markers)) ]
              coex =as.data.frame(as.matrix(coex))

              #prim.markers = c("Satb2","Bcl11b","Cux1","Fezf2","Tbr1")
              pval = get.pval(object = obj,gene.set.col = unique(c(unlist(markers.list),prim.markers)))
              pval.red = apply(pval, 1, FUN=min)
              genes.row = names(pval.red[pval.red < p_value])
              genes.row = unique(c(colnames(coex),genes.row))

              if(symmetric == TRUE){
                  coex = obj@coex

                  coex = as.data.frame(as.matrix(coex))

                  coex = coex[rownames(coex) %in% genes.row, colnames(coex) %in%  genes.row]

              }else{
                  coex = obj@coex
                  coex = as.data.frame(as.matrix(coex))

                  coex = coex[rownames(coex) %in% genes.row,]
              }

              #coex = coex#/sqrt(obj@n_cells)
              #coex = as.data.frame(as.matrix(coex))
              list.rows = c()
              for (m in unlist(markers.list)) {
                  genes = rownames(pval[pval[,m] < p_value,])
                  genes = genes[genes %in% rownames(coex[coex[,m] > 0,]) ]
                  list.rows[[m]]=genes

              }

              list.cols = c()
              for (m in prim.markers) {
                  genes = rownames(pval[pval[,m] < p_value,])
                  genes = genes[genes %in% rownames(coex[coex[,m] > 0,]) ]
                  list.cols[[m]]=genes

              }

              #coex_17.2 = coex_17.2[,colnames(coex_17.2) %in% rownames(coex_17.2)]
              #cl.genes.rows = cl.genes.cols
              cl.genes.rows = c()
              for (ll in names(list.rows)) {
                  tmp = data.frame("genes"=list.rows[[ll]],"cl"=rep(ll,length(list.rows[[ll]])))
                  cl.genes.rows = rbind(cl.genes.rows,tmp)
              }

              cl.genes.rows = cl.genes.rows[cl.genes.rows$genes %in% rownames(coex),]

              reorder_idx_row <- match(cl.genes.rows$gene,rownames(coex))


              if (symmetric == TRUE) {
                  cl.genes.cols = data.frame()
                  for (ll in names(list.rows)) {
                      tmp = data.frame("genes"=list.rows[[ll]],"cl"=rep(ll,length(list.rows[[ll]])))
                      cl.genes.cols = rbind(cl.genes.cols,tmp)
                  }
              }else{
                  cl.genes.cols = data.frame()
                  for (ll in names(list.cols)) {
                      tmp = data.frame("genes"=list.cols[[ll]],"cl"=rep(ll,length(list.cols[[ll]])))
                      cl.genes.cols = rbind(cl.genes.cols,tmp)
                  }

              }
              #cl.genes = cl.genes[cl.genes$gene %in% colnames(coex_17.2),]
              cl.genes.cols = cl.genes.cols[cl.genes.cols$genes %in% colnames(coex),]

              reorder_idx_col <- match(cl.genes.cols$gene,colnames(coex))


              to.plot <- coex[reorder_idx_row,reorder_idx_col]

              col_fun = circlize::colorRamp2(c(round(stats::quantile(as.matrix(to.plot),probs =0.001),
                                                     digits = 3),0,
                                     round(stats::quantile(as.matrix(to.plot),probs =0.999),
                                           digits = 3)),
                                   c("#E64B35FF", "gray93", "#3C5488FF"))

              #The next line is to set the columns and raws order
              # need to be implemented
              #cl.genes.rows$cl =factor(cl.genes.rows$cl,c("Reln","Satb2","Sox5","Bcl11b"))

              part1 = ComplexHeatmap::Heatmap(as.matrix(to.plot),
                              cluster_rows = FALSE,
                              cluster_columns = FALSE ,
                              row_split = cl.genes.rows$cl,
                              column_split = cl.genes.cols$cl ,
                              col = col_fun,
                              show_row_names = FALSE,
                              show_column_names = FALSE,
                              column_title_gp = grid::gpar(fill = "#8491B44C", font = 3,
                                                           col= "#3C5488FF"),
                              row_title_gp = grid::gpar(fill = "#8491B44C",font = 3, col= "#3C5488FF"))
              lgd = ComplexHeatmap::Legend(col_fun = col_fun, title = "coex",grid_width =
                                               unit(0.3, "cm"),
                           direction = "horizontal", title_position = "topcenter",
                           title_gp = grid::gpar(fontsize = 10, fontface = "bold",col="#3C5488FF"),
                           labels_gp = grid::gpar(col = "#3C5488FF", font = 3) )

              #part1 =
                  ComplexHeatmap::draw(part1,show_heatmap_legend = FALSE,
                           annotation_legend_list = lgd,annotation_legend_side = "bottom")


              #return(part1)
          }
)



#' get.gene.coexpression.space
#'
#'To make the GDI more specific, it may be desirable to restrict the set of genes against which
#'GDI is computed to a selected subset V with the recommendation to include a
#'consistent fraction of cell-identity genes, and possibly focusing on markers specific
#'for the biological question of interest (for instance neural cortex layering markers).
#'In this case we denote it as local differentiation index (LDI) relative to V.
#'
#' @param object The COTAN object.
#' @param n.genes.for.marker The number of genes correlated with the primary markers that
#' we want to consider.
#' By default this is 25.
#' @param primary.markers A vector of primary marker names.
#'
#' @return A dataframe
#' @export
#' @importFrom stats quantile
#' @importFrom Matrix rowSums
#' @importFrom Matrix colMeans
#' @importFrom Matrix forceSymmetric
#' @rdname get.gene.coexpression.space
#' @examples
#' data("ERCC.cotan")
#' df = get.gene.coexpression.space(ERCC.cotan, n.genes.for.marker = 10,
#' primary.markers=rownames(ERCC.cotan@raw[sample(nrow(ERCC.cotan@raw),5),]))
setGeneric("get.gene.coexpression.space", function(object, n.genes.for.marker = 25, primary.markers)
    standardGeneric("get.gene.coexpression.space"))
#' @rdname get.gene.coexpression.space
setMethod("get.gene.coexpression.space","scCOTAN",
          #da sistemare
          function(object, n.genes.for.marker = 25, primary.markers) {
              print("calculating gene coexpression space: output tanh of reduced coex matrix" )

              object@coex = Matrix::forceSymmetric(object@coex, uplo="L" )

              p.val.matrix = get.pval(object,gene.set.col = primary.markers)
              if(!length(primary.markers) == ncol(p.val.matrix)){
                  print(paste("Gene", primary.markers[!primary.markers %in% colnames(p.val.matrix)],
                              "not present!", sep = " "))
                  primary.markers = primary.markers[primary.markers %in% colnames(p.val.matrix)]
              }


              all.genes.to.an = vector()
              for (m in primary.markers) {
                  #print(m)
                  #tm =rownames(p.val.matrix[order(p.val.matrix[,m]),])[1:n.genes.for.marker]
                  tm =rownames(p.val.matrix[order(p.val.matrix[,m]),])[seq_len(n.genes.for.marker)]
                  all.genes.to.an = c(all.genes.to.an,tm)
                  all.genes.to.an =unique(all.genes.to.an)
              }
              all.genes.to.an =unique(c(all.genes.to.an,primary.markers))
              print(paste0("Secondary markers:", length(all.genes.to.an)))
              tmp = p.val.matrix[all.genes.to.an,]
              for (m in primary.markers) {
                  tmp = as.data.frame(tmp[order(tmp[,m]),])
                  #tmp$rank = c(1:nrow(tmp))
                  tmp$rank = c(seq_len(nrow(tmp)))
                  colnames(tmp)[ncol(tmp)] = paste("rank",m,sep = ".")
              }
              rank.genes = tmp[,(length(primary.markers)+1):ncol(tmp)]
            #for (c in c(1:length(colnames(rank.genes)))) {
            for (c in seq_along(colnames(rank.genes))) {
                colnames(rank.genes)[c] =strsplit(colnames(rank.genes)[c], split='.',
                                                    fixed = TRUE)[[1]][2]
                }

              S = get.S(object)

              S = S[,colnames(S) %in% all.genes.to.an]
              S = as.matrix(S)


              # This is the LDI 5%
              CD.sorted <- t(apply(t(S),2,sort,decreasing=TRUE))
              #CD.sorted = CD.sorted[,1:round(ncol(CD.sorted)/5, digits = 0)] #20
              CD.sorted = CD.sorted[,1:round(ncol(CD.sorted)/10, digits = 0)] #20
              CD.sorted = pchisq(as.matrix(CD.sorted), df=1, lower.tail=FALSE)

              quant.p.val2 = rowMeans(CD.sorted)
              quant.p.val2 =as.data.frame(quant.p.val2)
              colnames(quant.p.val2) = "loc.GDI"

              quant.p.val2$names = rownames(quant.p.val2)

              quant.p.val2 = quant.p.val2[quant.p.val2$loc.GDI <=
                                              stats::quantile(quant.p.val2$loc.GDI,probs = 0.1),] #0.9
              genes.raw = quant.p.val2$names #0.9
              #genes.raw = unique(c(genes.raw, all.genes.to.an ))
              # just to color

              to.plot_cl.genes = object@coex[rownames(object@coex) %in%
                                                 genes.raw,colnames(object@coex) %in% all.genes.to.an]

              to.plot_cl.genes = to.plot_cl.genes * sqrt(object@n_cells)

              to.plot_cl.genes = tanh(to.plot_cl.genes)
              print(paste0("Columns (V set) number: ",dim(to.plot_cl.genes)[2],
                           " Rows (U set) number: ",dim(to.plot_cl.genes)[1]))

              return(to.plot_cl.genes)
          }
)


#' get.GDI
#'
#' This function produce a dataframe with the GDI for each genes.
#'
#' @param object A COTAN object
#' @param type Type of statistic to be used. Default is "S":
#' Pearson's chi-squared test statistics. "G" is G-test statistics
#'
#' @return A dataframe
#' @export
#' @importFrom stats pchisq
#' @importFrom Matrix rowSums
#' @importFrom Matrix colMeans
#' @importFrom Matrix forceSymmetric
#' @rdname get.GDI
#' @examples
#' data("ERCC.cotan")
#' quant.p = get.GDI(ERCC.cotan)
setGeneric("get.GDI", function(object,type="S") standardGeneric("get.GDI"))
#' @rdname get.GDI
setMethod("get.GDI","scCOTAN",
          function(object,type="S") {

              object@coex = Matrix::forceSymmetric(object@coex, uplo="L" )

              print("function to generate GDI dataframe")
              if (type=="S") {
                  print("Using S")
                  S = get.S(object)
              }else if(type=="G"){
                  print("Using G")
                  S = get.G(object)
              }


              S = as.data.frame(as.matrix(S))
              #G = as.data.frame(as.matrix(G))
              CD.sorted <- apply(S,2,sort,decreasing=TRUE)
              CD.sorted = CD.sorted[1:round(nrow(CD.sorted)/20, digits = 0),]
              CD.sorted = pchisq(as.matrix(CD.sorted), df=1, lower.tail=FALSE)



              GDI = colMeans(CD.sorted)
              #GDI = colMeans(CD.sorted)
              GDI =as.data.frame(GDI)
              colnames(GDI) = "mean.pval"



              sum.raw.norm = log(rowSums(as.matrix(object@raw.norm)))

              cells=as.matrix(object@raw)
              cells[cells > 0] <- 1
              cells[cells <= 0] <- 0

              exp.cells = (rowSums(cells)/object@n_cells)*100

              GDI =  merge(GDI, as.data.frame(sum.raw.norm), by="row.names",all.x=TRUE)
              rownames(GDI)= GDI$Row.names
              GDI =  GDI[,2:ncol(GDI)]
              GDI =  merge(GDI, as.data.frame(exp.cells), by="row.names",all.x=TRUE)
              rownames(GDI)= GDI$Row.names
              GDI$log.mean.p = -log(GDI$mean.pval)
              GDI$GDI = log(GDI$log.mean.p)
              GDI = GDI[,c("sum.raw.norm","GDI","exp.cells")]

              return(GDI)

          }
)



#' plot_GDI
#'
#' This function directly evaluate and plot the GDI for a sample.
#'
#' @param object A COTAN object
#' @param cond A string corresponding to the condition/sample (it is used only for the title)
#' @param type Type of statistic to be used. Default is "S":
#' Pearson's chi-squared test statistics. "G" is G-test statistics
#'
#' @return A ggplot2 object
#' @export
#' @import ggplot2
#' @importFrom  stats quantile
#' @importFrom Matrix forceSymmetric
#' @rdname plot_GDI
#' @examples
#' data("ERCC.cotan")
#' plot_GDI(ERCC.cotan, cond = "ERCC")
setGeneric("plot_GDI", function(object, cond,type="S") standardGeneric("plot_GDI"))
#' @rdname plot_GDI
setMethod("plot_GDI","scCOTAN",
          function(object, cond,type="S") {

              object@coex = Matrix::forceSymmetric(object@coex, uplo="L" )

              print("GDI plot ")
              if (type=="S") {
                  GDI = get.GDI(object,type="S")
              }else if(type=="G"){
                  print("Using G")
                  GDI = get.GDI(object,type="G")
              }


              si = 12
              plot =  ggplot(GDI, aes(x = sum.raw.norm, y = GDI)) +
                  geom_point(size = 2, alpha=0.5, color= "#8491B4B2") +
                  geom_hline(yintercept=1.5, linetype="dotted", color = "darkred", size =1) +
                  geom_hline(yintercept=stats::quantile(GDI$GDI)[4], linetype="dashed",
                             color = "darkblue") +
                  geom_hline(yintercept=stats::quantile(GDI$GDI)[3], linetype="dashed",
                             color = "darkblue") +
                  xlab("log normalized reads sum")+
                  ylab("global p val index (GDI)")+
                  ggtitle(paste("GDI ",cond, sep = " "))+
                  theme(axis.text.x = element_text(size = si, angle = 0, hjust = .5, vjust = .5,
                                                   face = "plain", colour ="#3C5488FF" ),
                        axis.text.y = element_text( size = si, angle = 0, hjust = 0, vjust = .5,
                                                    face = "plain", colour ="#3C5488FF"),
                        axis.title.x = element_text( size = si, angle = 0, hjust = .5, vjust = 0,
                                                     face = "plain", colour ="#3C5488FF"),
                        axis.title.y = element_text( size = si, angle = 90, hjust = .5, vjust = .5,
                                                     face = "plain", colour ="#3C5488FF"),
                        legend.title = element_blank(),
                        plot.title = element_text(color="#3C5488FF", size=14, face="bold.italic"),
                        legend.text = element_text(color = "#3C5488FF",face ="italic" ),
                        legend.position = "bottom")  # titl)

              return(plot)

          }
)


#' automatic.COTAN.object.creation
#'
#' @param df dataframe with the row counts
#' @param out_dir directory for the output
#' @param GEO GEO or other code that identify the dataset
#' @param sc.method Type of single cell RNA-seq method used
#' @param cond A string that will identify the sample or condition. It will be part of the
#' final file name.
#' @param cores number of cores to be used
#' @param mt A boolean (default F). If T mitochondrial  genes will be kept in the analysis,
#' otherwise they will be removed.
#' @param mt_prefix is the prefix that identify the mitochondrial  genes (default is the mouse
#' prefix: "^mt")
#' @return It return the COTAN object. It will also store it directly in the output directory
#' @export
#'
#' @importFrom Matrix rowSums
#' @importFrom Matrix colSums
#' @importFrom utils write.csv
#' @importFrom methods new
#' @import grDevices
#' @import ggplot2
#' @import ggrepel
#' @rdname automatic.COTAN.object.creation
#' @examples
#'
#' data("raw.dataset")
#' obj = automatic.COTAN.object.creation(df= raw,
#' out_dir =  tempdir(),
#' GEO = "test_GEO",
#' sc.method = "test_method",
#' cond = "test")
#'
setGeneric("automatic.COTAN.object.creation", function(df, out_dir, GEO, sc.method,
                                                       cond, mt = FALSE, mt_prefix="^mt", cores = 1)
    standardGeneric("automatic.COTAN.object.creation"))
#' @rdname automatic.COTAN.object.creation
setMethod("automatic.COTAN.object.creation","data.frame",
          function(df, out_dir, GEO, sc.method, cond, mt = FALSE, mt_prefix="^mt", cores = 1) {
              start_time_all <- Sys.time()

              mycolours <- c("A" = "#8491B4B2","B"="#E64B35FF")
              my_theme = theme(axis.text.x = element_text(size = 14, angle = 0, hjust = .5,
                                                          vjust = .5,
                                                          face = "plain", colour ="#3C5488FF" ),
                               axis.text.y = element_text( size = 14, angle = 0, hjust = 0,
                                                           vjust = .5,
                                                           face = "plain", colour ="#3C5488FF"),
                               axis.title.x = element_text( size = 14, angle = 0, hjust = .5,
                                                            vjust = 0,
                                                            face = "plain", colour ="#3C5488FF"),
                               axis.title.y = element_text( size = 14, angle = 90, hjust = .5,
                                                            vjust = .5,
                                                            face = "plain", colour ="#3C5488FF"))
              obj = methods::new("scCOTAN",raw = df)
              obj = initRaw(obj,GEO = GEO ,sc.method = sc.method,cond = cond)
              if (mt == FALSE) {
                  genes_to_rem = rownames(obj@raw[grep(mt_prefix, rownames(obj@raw)),])
                  obj@raw = obj@raw[!rownames(obj@raw) %in% genes_to_rem,]
                  cells_to_rem = colnames(obj@raw[which(colSums(obj@raw) == 0)])
                  obj@raw = obj@raw[,!colnames(obj@raw) %in% cells_to_rem]
              }
              t = cond

              print(paste("Condition ",t,sep = ""))
              #--------------------------------------
              n_cells = length(colnames(obj@raw))
              print(paste("n cells", n_cells, sep = " "))

              n_it = 1

              if(!file.exists(out_dir)){
                  dir.create(file.path(out_dir))
              }

              if(!file.exists(paste(out_dir,"cleaning", sep = ""))){
                  dir.create(file.path(out_dir, "cleaning"))
              }

              ttm = clean(obj)

              obj = ttm$object

              #pdf(paste(out_dir,"cleaning/",t,"_",n_it,"_plots_without_cleaning.pdf", sep = ""))
              pdf(file.path(out_dir,"cleaning",paste(t,"_",n_it,"_plots_without_cleaning.pdf",
                                                     sep = "")))
              ttm$pca.cell.2
              ggplot(ttm$D, aes(x=n,y= means)) + geom_point() +
                  geom_text_repel(data=subset(ttm$D, n > (max(ttm$D$n)- 15) ),
                                  aes(n, means, label=rownames(ttm$D[ttm$D$n > (max(ttm$D$n)- 15),])),
                                  nudge_y      = 0.05,
                                  nudge_x      = 0.05,
                                  direction    = "x",
                                  angle        = 90,
                                  vjust        = 0,
                                  segment.size = 0.2)+
                  ggtitle(label = "B cell group genes mean expression", subtitle =
                              " - B group NOT removed -")+
                  my_theme + theme(plot.title = element_text(color = "#3C5488FF",
                                                  size = 20, face = "italic",
                                                  vjust = - 10,hjust = 0.02 ),
                        plot_subtitle = element_text(color = "darkred",vjust = - 15,hjust = 0.01 ))

              dev.off()

              nu_est = round(obj@nu, digits = 7)

              plot_nu <-ggplot(ttm$pca_cells,aes(x=PC1,y=PC2, colour = log(nu_est)))

              plot_nu = plot_nu + geom_point(size = 1,alpha= 0.8)+
                  scale_color_gradient2(low = "#E64B35B2",mid =  "#4DBBD5B2", high =  "#3C5488B2" ,
                                        midpoint = log(mean(nu_est)),name = "ln(nu)")+
                  ggtitle("Cells PCA coloured by cells efficiency") +
                  my_theme +  theme(plot.title = element_text(color = "#3C5488FF", size = 20),
                                    legend.title=element_text(color = "#3C5488FF", size = 14,
                                                              face = "italic"),
                                    legend.text = element_text(color = "#3C5488FF", size = 11),
                                    legend.key.width = unit(2, "mm"),
                                    legend.position="right")

              #pdf(paste(out_dir,"cleaning/",t,"_plots_PCA_efficiency_colored.pdf", sep = ""))
              pdf(file.path(out_dir,"cleaning",paste(t,"_plots_PCA_efficiency_colored.pdf", sep = "")))
              plot_nu
              dev.off()

              nu_df = data.frame("nu"= sort(obj@nu), "n"=c(1:length(obj@nu)))

              #pdf(paste(out_dir,"cleaning/",t,"_plots_efficiency.pdf", sep = ""))
              pdf(file.path(out_dir,"cleaning",paste(t,"_plots_efficiency.pdf", sep = "")))
              ggplot(nu_df, aes(x = n, y=nu)) + geom_point(colour = "#8491B4B2", size=1) +my_theme +
                  annotate(geom="text", x=50, y=0.25, label="nothing to remove ", color="darkred")
              dev.off()

              analysis_time <- Sys.time()

              print("Cotan analysis function started")
              obj = cotan_analysis(obj,cores = cores)

              coex_time <- Sys.time()
              analysis_time = difftime(Sys.time(), analysis_time, units = "mins")
              #analysis_time <- coex_time - analysis_time

              print(paste0("Only analysis time ", analysis_time))

              print("Cotan coex estimation started")
              obj = get.coex(obj)

              end_time <- Sys.time()

              #all.time = end_time - start_time_all
              all.time =  difftime(end_time, start_time_all, units = "mins")
              print(paste0("Total time ", all.time))

              #coex_time = end_time - coex_time
              coex_time = difftime(end_time, coex_time, units = "mins")

              print(paste0("Only coex time ",coex_time))
              # saving the structure
              utils::write.csv(data.frame("type" = c("tot_time","analysis_time","coex_time"),
                                   "times"=
                                       c(as.numeric(all.time),
                                         as.numeric(analysis_time),
                                         as.numeric(coex_time) ),
                                   "n.cells"=n_cells,"n.genes"=dim(obj@raw)[1]),
                        file = file.path(out_dir, paste(t,"_times.csv", sep = "")))

              print(paste0("Saving elaborated data locally at ", out_dir,t,".cotan.RDS"))
              saveRDS(obj,file = file.path(out_dir,paste(t,".cotan.RDS", sep = "")))

              return(obj)

          }
)


#' get.observed.ct
#'
#' This function is used to get the observed contingency table for a given pair of genes.
#'
#' @param object The cotan object
#' @param g1 A gene
#' @param g2 The other gene
#'
#' @return A contingency table as dataframe
#' @export
#' @rdname get.observed.ct
#' @examples
#' data("ERCC.cotan")
#' g1 = rownames(ERCC.cotan@raw)[sample(nrow(ERCC.cotan@raw), 1)]
#' g2 = rownames(ERCC.cotan@raw)[sample(nrow(ERCC.cotan@raw), 1)]
#' get.observed.ct(object = ERCC.cotan, g1 = g1, g2 = g2)
setGeneric("get.observed.ct", function(object,g1,g2) standardGeneric("get.observed.ct"))
#' @rdname get.observed.ct
setMethod("get.observed.ct","scCOTAN",
          function(object,g1,g2) {
              #cells=obj@raw
              #---------------------------------------------------
              # Cells matrix : formed by row data matrix changed to 0-1 matrix
              #cells[cells > 0] <- 1
              #cells[cells <= 0] <- 0

              si_si = object@yes_yes[g1,g2]
              n_cells = object@n_cells

              si_any = max(object@yes_yes[g1,])

              any_si = max(object@yes_yes[g2,])
              si_no = si_any - si_si
              no_si = any_si - si_si
              no_no = n_cells - (si_si + si_no + no_si)

              ct = as.data.frame(matrix(ncol = 2,nrow = 2))
              colnames(ct) = c(paste(g1,"yes",sep = "."),paste(g1,"no",sep = "."))
              rownames(ct) = c(paste(g2,"yes",sep = "."),paste(g2,"no",sep = "."))
              ct[1,1] = si_si
              ct[1,2] = no_si
              ct[2,2] = no_no
              ct[2,1] = si_no

              return(ct)

          }
)


#' get.expected.ct
#'
#' This function is used to get the expected contingency table for a given pair of genes.
#' @param object The cotan object
#' @param g1 A gene
#' @param g2 The other gene
#'
#' @return A contingency table as dataframe
#' @export
#'
#' @importFrom Matrix rowSums
#' @importFrom Matrix colSums
#' @rdname get.expected.ct
#' @examples
#' data("ERCC.cotan")
#' g1 = rownames(ERCC.cotan@raw)[sample(nrow(ERCC.cotan@raw), 1)]
#' g2 = rownames(ERCC.cotan@raw)[sample(nrow(ERCC.cotan@raw), 1)]
#' while (g1 %in% ERCC.cotan@hk) {
#'g1 = rownames(ERCC.cotan@raw)[sample(nrow(ERCC.cotan@raw), 1)]
#'}
#'
#'while (g2 %in% ERCC.cotan@hk) {
#'    g2 = rownames(ERCC.cotan@raw)[sample(nrow(ERCC.cotan@raw), 1)]
#'}
#' get.expected.ct(object = ERCC.cotan, g1 = g1, g2 = g2)
setGeneric("get.expected.ct", function(object,g1,g2) standardGeneric("get.expected.ct"))
#' @rdname get.expected.ct
setMethod("get.expected.ct","scCOTAN",
          function(object,g1,g2) {
              fun_pzero <- function(a,mu){
                  #a= as.numeric(a[,2])
                  #print(a)
                  (a <= 0)*(exp(-(1+abs(a))*mu)) + (a > 0)*(1+abs(a)*mu)^(-1/abs(a))
              }

              if(g1 %in% object@hk | g2 %in% object@hk ){
                  #print("Error. A gene is constitutive!")
                  #break()
                  stop("Error. A gene is constitutive!")
              }
              mu_estimator = object@lambda[c(g1,g2)] %*% t(object@nu)

              cells=as.matrix(object@raw)[c(g1,g2),]
              #---------------------------------------------------
              # Cells matrix : formed by row data matrix changed to 0-1 matrix
              cells[cells > 0] <- 1
              cells[cells <= 0] <- 0

              M = fun_pzero(object@a[c(g1,g2)],mu_estimator)
              rownames(M) = c(g1,g2)
              N = 1-M #fun_pzero(as.numeric(tot2[,2]),mu_estimator[,colnames(cells)])

              n_zero_esti = rowSums(M) # estimated number of zeros for each genes
              n_zero_obs = rowSums(cells[!rownames(cells) %in% object@hk,] == 0)
              # observed number of zeros for each genes
              dist_zeros = sqrt(sum((n_zero_esti - n_zero_obs)^2))



              if(any(is.na(M))){
                  #print(paste("Errore: some Na in matrix M", which(is.na(M),arr.ind = TRUE),sep = " "))
                  #break()
                  stop(paste("Errore: some Na in matrix M", which(is.na(M),arr.ind = TRUE),sep = " "))
              }

              gc()
              estimator_no_no = M %*% t(M)
              estimator_no_si = M %*% t(N)
              estimator_si_no = t(estimator_no_si)
              estimator_si_si = N %*% t(N)


              ct = as.data.frame(matrix(ncol = 2,nrow = 2))
              colnames(ct) = c(paste(g1,"yes",sep = "."),paste(g1,"no",sep = "."))
              rownames(ct) = c(paste(g2,"yes",sep = "."),paste(g2,"no",sep = "."))
              ct[1,1] = estimator_si_si[g1,g2]
              ct[1,2] = estimator_no_si[g1,g2]
              ct[2,2] = estimator_no_no[g1,g2]
              ct[2,1] = estimator_si_no[g1,g2]

              return(ct)

          }
)

