## Main functions

#' scCOTAN-class
#'
#' Define my COTAN structure
#' @slot raw ANY. To store the raw data matrix
#' @slot raw.norm ANY. To store the raw data matrix divided for the cell efficiency estimated (nu)
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
#' @return
#' @export
#' @examples
#' obj = new("scCOTAN",raw=read.csv("tests/testthat/raw.csv", header = T,row.names = 1 )) # with raw the raw data matrix with cells on colums and genes on rows
#'
setClass("scCOTAN", slots =
                        c(raw="ANY",raw.norm="ANY", coex="ANY",
                          nu="vector",lambda="vector",a="vector", hk = "vector", n_cells = "numeric",
                          meta="data.frame",yes_yes="ANY", clusters="vector",cluster_data="data.frame")) -> scCOTAN

#' initRaw
#'
#' It starts to fill some fields of cotan object.
#' @param object the dataframe containing the raw data
#' @param GEO a code reporting the GEO identification or other specific dataset code
#' @param sc.method a string reporting the method used for the sequencing
#' @param cond a string reporting the specific sampe condition or time point
#'
#' @return A COTAN object to store all information
#' @export
#' @examples
#' obj = new("scCOTAN",raw=read.csv("tests/testthat/raw.csv", header = T,row.names = 1 ))
#' obj = initRaw(obj,GEO="GSE95315" ,sc.method="10X",cond = "mouse dentate gyrus P0 subclustering of cl.2 - 0,1,2")
#'
setGeneric("initRaw", function(object,GEO,sc.method="10X", cond) standardGeneric("initRaw"))
setMethod("initRaw","scCOTAN",
          function(object,GEO,sc.method,cond) {
              print("Initializing S4 object")
              if(!attr(object@raw, "class")[1] == "dgCMatrix" | is.null(attr(object@raw, "class")) ){
                  object@raw = as(as.matrix(object@raw), "sparseMatrix")
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
#' Main function that can be used to check and clean the dataset. It also produce (using the function fun_linear) and store
#' the estimators for nu and lambda. It also fill the raw.norm (raw / nu) and n_cell (the initial number of cells in the dataset)
#' @param object COTAN object
#' @return a list of objects containing: "cl1" is the first cell cluster, "cl2" is the second cell cluster,
#' "pca.cell.2" is a ggplot2 cell pca plot, "object" is the COTAN object with saved the estimated lambda and mu, "mu_estimator", "D"
#' "pca_cells" pca numeric data.
#' @export
#'
#' @import ggplot2
#' dplyr
#'
#' @importFrom tibble rownames_to_column
#' @importFrom stats hclust
#' @importFrom stats cutree
#' @importFrom parallel mclapply
#' @importFrom Matrix rowMeans
#' @importFrom Matrix colMeans
#'
#' @examples
#' \dontrun{
#' ttm = clean(obj, cells)
#' obj = ttm$object
#' ttm$pca.cell.2
#' }
setGeneric("clean", function(object) standardGeneric("clean"))
setMethod("clean","scCOTAN",
          #clean_function <-
              function(object) {

                  mycolours <- c("A" = "#8491B4B2","B"="#E64B35FF")

                  my_theme = theme(axis.text.x = element_text(size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF" ),
                                   axis.text.y = element_text( size = 14, angle = 0, hjust = 0, vjust = .5, face = "plain", colour ="#3C5488FF"),
                                   axis.title.x = element_text( size = 14, angle = 0, hjust = .5, vjust = 0, face = "plain", colour ="#3C5488FF"),
                                   axis.title.y = element_text( size = 14, angle = 90, hjust = .5, vjust = .5, face = "plain", colour ="#3C5488FF"))


                cells = get.zero_one.cells(object)
              #cells = cells[rowSums(cells)> round((length(colnames(raw))/1000*1), digits = 0),] # 3 per mille delle cellule
              object@raw = object@raw[rownames(cells), colnames(cells)]
              #fwrite(cells,paste(out_dir,"cells_",t,".csv", sep = ""))
              #save(n_it,file = paste(out_dir,"n_it_",t, sep = ""))
              #if(iter){
               #   list1 = fun_linear_iter(object, mean_type)
              #}else{
                  list1 = fun_linear(object)
              #}

              #object@meta[(nrow(object@meta)+1),1:2] = c("mean type:",mean_type)
              #object@meta[(nrow(object@meta)+1),1:2] = c("Iteration:",iter)
              #object@meta[(nrow(object@meta)+1),1:2] = c("Type of estimation:","linear")
              # list1 = fun_iter_no_nu(cells = cells, raw = raw)
              gc()

              dist_cells = list1$dist_cells
              pca_cells = list1$pca_cells
              t_to_clust = as.matrix(list1$t_to_clust)
              mu_estimator = list1$mu_estimator
              object = list1$object
              to_clust = list1$to_clust


              raw_norm <- t(t(as.matrix(object@raw)) * (1/(as.vector(object@nu))))
              #raw_norm <- as.matrix(object@raw) %b/% matrix(as.vector(object@nu), nrow = 1)
              raw_norm <- as(as.matrix(raw_norm), "sparseMatrix")
              #end_time <- Sys.time()
              #print(paste("to clust; time",end_time - start_time, sep = " " ))
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
                  B = B[order(rowMeans(B[,2:length(colnames(B))]),decreasing = T), ]
              }else{
                  B = B[order(B[,2],decreasing = T), ]
              }

              #B = B[order(B[,2:length(colnames(B))],decreasing = T), ] #if just one column

              print(head(B, 15))

              C = arrange(B,rowMeans(B[2:length(colnames(B))]))
              rownames(C) = C$rowname
              D = data.frame("means"=rowMeans(C[2:length(colnames(C))]),"n"=NA )
              D = D[D$means>0,]
              D$n = c(1:length(D$means))

              # check if the pca plot is clean enought and from the printed genes, if the smalest group of cells are caratterised by particular genes

              pca_cells = cbind(pca_cells,"groups"=t_to_clust$groups)

              pca.cell.1 = ggplot(subset(pca_cells,groups == "A" ), aes(x=PC1, y=PC2,colour =groups)) +geom_point(alpha = 0.5, size=3)

              pca.cell.2=  pca.cell.1 + geom_point(data = subset(pca_cells, groups != "A" ), aes(x=PC1, y=PC2,colour =groups),alpha = 0.8, size=3)+
                  scale_color_manual("groups", values = mycolours)  +
                  my_theme + theme(legend.title = element_blank(),
                                   legend.text = element_text( size = 12,color = "#3C5488FF",face ="italic" ),
                                   legend.position="bottom")

              object@n_cells = length(colnames(object@raw))

              output = list("cl1"=cl1,"cl2"=cl2,"pca.cell.2"=pca.cell.2,"object"=object, "mu_estimator"=mu_estimator, "D"=D, "pca_cells"=pca_cells)
              return(output)
          }
)

#' cotan_analysis
#'
#' This is the main function that estimates the a vector to store all the negative binomial
#' dispersion factors. It need to be run after \code{\link{clean}}
#' @param object A COTAN object
#' @importFrom parallel mclapply
#' @return a COTAN object
#' @export
#' @examples
#' \dontrun{obj = cotan_analysis(obj)
#' }
setGeneric("cotan_analysis", function(object) standardGeneric("cotan_analysis"))
setMethod("cotan_analysis","scCOTAN",
          function(object) {

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
                                       fun_my_opt, ce.mat= cells, mu_est= mu_estimator  ,mc.cores = 11)

                  }else{
                      print("Final trance!")
                      tot1 = parallel::mclapply(rownames(mu_estimator)[p:length(rownames(mu_estimator))],
                                      fun_my_opt, ce.mat=cells, mu_est= mu_estimator  ,mc.cores = 11)
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
              print(paste("a min:",min(tot2$a) ,"| a max",max(tot2$a) , "| negative a %:",sum(tot2$a <0)/nrow(tot2)*100 ,sep=" "))

              gc()

              return(object)
          }
)


setGeneric("hk_genes", function(object) standardGeneric("hk_genes"))
setMethod("hk_genes","scCOTAN",
          function(object) {
              print("save effective constitutive genes")
              cells=as.matrix(object@raw)
              #---------------------------------------------------
              # Cells matrix : formed by row data matrix changed to 0-1 matrix
              cells[cells > 0] <- 1
              cells[cells <= 0] <- 0

              hk = names(which(rowSums(cells) == length(colnames(cells))))
              object@hk = hk

              return(object)
          }
)


#' get.coex
#'
#' This function estimates and stores the coex*sqrt(n_cells) matrix in the coex field.
#' It need to be run after \code{\link{cotan_analysis}}
#' @param object
#'
#' @return It returns a COTAN object
#' @export
#'
#' @examples
#' \dontrun{
#' obj = get.coex(obj)
#' }
setGeneric("get.coex", function(object) standardGeneric("get.coex"))
setMethod("get.coex","scCOTAN",
          function(object) {
              print("coex dataframe creation")
              hk = object@hk
              ll = obs_ct(object)

              object = ll$object

              ll$no_yes= ll$no_yes[!rownames(ll$no_yes) %in% hk,!colnames(ll$no_yes) %in% hk]
              ll$no_no = ll$no_no[!rownames(ll$no_no) %in% hk,!colnames(ll$no_no) %in% hk]
              si_si = object@yes_yes[!rownames(object@yes_yes) %in% hk,!colnames(object@yes_yes) %in% hk]
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

              coex = coex / sqrt(1/new_estimator_si_si + 1/new_estimator_no_no + 1/new_estimator_si_no + 1/new_estimator_no_si)

              #print("coex low values substituted with 0")
              #coex[abs(coex) <= 1]=0
              print("Diagonal coex values substituted with 0")
              diag(coex) = 0
              coex = coex / sqrt(object@n_cells)
              object@coex = coex

              return(object)
          }
)




setGeneric("obs_ct", function(object) standardGeneric("obs_ct"))
setMethod("obs_ct","scCOTAN",
          function(object) {
              cells=as.matrix(object@raw)
              #---------------------------------------------------
              # Cells matrix : formed by row data matrix changed to 0-1 matrix
              cells[cells > 0] <- 1
              cells[cells <= 0] <- 0
              cells = as.matrix(cells)
              print("Generating contingency tables for observed data")
              somma = rowSums(cells)
              somma = as.matrix(somma)
              si_any = do.call("cbind", replicate(length(rownames(somma)), somma, simplify = FALSE))
              #rm(somma) = object@yes_yes
              colnames(si_any) = rownames(si_any)
              if (is.null(object@yes_yes)) {
                  object = obs_yes_yes(object)
              }

              si_si = as.matrix(object@yes_yes)
              si_no = si_any - si_si
              #si_no = as(si_no, "sparseMatrix")
              si_any = t(si_any)
              no_si = si_any - si_si
              #rm(si_any)
              no_no = length(colnames(cells)) - (si_si + no_si + si_no)
              #no_no = as(no_no, "sparseMatrix")
              out = list("no_yes"=no_si,"yes_no"=si_no,"no_no"=no_no,"object"=object)
              return(out)

          }
)

setGeneric("obs_yes_yes", function(object) standardGeneric("obs_yes_yes"))
setMethod("obs_yes_yes","scCOTAN",
          function(object) {
              print("creation of observed yes/yes contingency table")

              cells=as.matrix(object@raw)
              #---------------------------------------------------
              # Cells matrix : formed by row data matrix changed to 0-1 matrix
              cells[cells > 0] <- 1
              cells[cells <= 0] <- 0

              si_si = as.matrix(cells) %*% t(as.matrix(cells))
              object@yes_yes = as(as.matrix(si_si), "sparseMatrix")
              return(object)
          }
)

setGeneric("expected_ct", function(object) standardGeneric("expected_ct"))
setMethod("expected_ct","scCOTAN",
          function(object) {


              cells=as.matrix(object@raw)
              #---------------------------------------------------
              # Cells matrix : formed by row data matrix changed to 0-1 matrix
              cells[cells > 0] <- 1
              cells[cells <= 0] <- 0

              mu_estimator = mu_est(object)
              mu_estimator = mu_estimator[!rownames(mu_estimator) %in% object@hk,]
              print("expected contingency tables creation")

              M = fun_pzero(object@a,mu_estimator[,colnames(cells)])
              N = 1-M #fun_pzero(as.numeric(tot2[,2]),mu_estimator[,colnames(cells)])

              n_zero_esti = rowSums(M) # estimated number of zeros for each genes
              n_zero_obs = rowSums(cells[!rownames(cells) %in% object@hk,] == 0) # observed number of zeros for each genes
              dist_zeros = sqrt(sum((n_zero_esti - n_zero_obs)^2))

              print(paste("The distance between estimated n of zeros and observed number of zero is", dist_zeros,"over", length(rownames(M)), sep = " "))

              if(any(is.na(M))){
                  print(paste("Errore: some Na in matrix M", which(is.na(M),arr.ind = T),sep = " "))
                  break()
              }

              gc()
              estimator_no_no = M %*% t(M)
              estimator_no_si = M %*% t(N)
              estimator_si_no = t(estimator_no_si)
              estimator_si_si = N %*% t(N)

              print("Done")
              #rm(M)
              #rm(N)

              out = list("estimator_no_no"=estimator_no_no, "estimator_no_yes"=estimator_no_si,"estimator_yes_no"=estimator_si_no,"estimator_yes_yes"=estimator_si_si)

              return(out)
          }
)

#' get.pval
#'
#' This function computes the p-values for genes in the COTAN oblect. It can be used genome-wide
#' or setting some specific genes of interest. By default it computes the p-values using the S
#' statistics (\eqn{\Chi^2})
#' @param object a COTAN object
#' @param gene.set.col an array of genes. It will be put in columns.
#' If left empty the function will do it genome-wide.
#' @param gene.set.row an array of genes. It will be put in rows.
#' If left empty the function will do it genome-wide.
#' @param type_stat By default it computes the S (\eqn{\Chi^{2}})
#'
#' @return a p-value matrix
#' @export
#'
#' @examples
#' \dontrun{
#' obj = get.pval(obj,type_stat="S")
#' }
setGeneric("get.pval", function(object, gene.set.col=c(),gene.set.row=c(), type_stat="S" ) standardGeneric("get.pval"))
setMethod("get.pval","scCOTAN",
          function(object, gene.set.col=c(),gene.set.row=c(), type_stat="S") {

              print(gene.set.col)

              if (!is.null(gene.set.row)) {

                  # a set for rows, not Genome Wide
                  cond.row = "on a set of genes on rows"
                  if (is.null(gene.set.col)) {
                      print("Error: can't have genome wide on columns and not rows! Use a subset on gene.set.col, not on rows.")
                      break()
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
            }else{
                print("Using function G")
                S = get.G(object)
            }


              if(cond.col == "on a set of genes on columns"){
                  S = S[,colnames(S) %in% gene.set.col]
              }
              if(cond.row == "on a set of genes on rows"){
                  S = S[rownames(S) %in% gene.set.row,]
              }


              p_value = pchisq(as.matrix(S), df=1, lower.tail=F)


              return(p_value)
          }
)


