#' get.gene.clusters
#'
#' This function perfor the gene clustering based on a pool of gene marker, using the gene
#' coexpression space
#'
#' @param obj a COTAN object
#' @param list.group.markers a named list with a name for each group a one or more marker gene for each group.
#' @param n.markers number of correlated genes to keep as other markers (default 25)
#' @param k.cuts number of estimated cluster (this define the high for the tree cut)
#' @param distance type of distance to use (default "cosine"... "euclidean" is also available)
#' @param hclust.method default is "ward.D2" but can be any method defined by hclust function
#' @import factoextra
#' @importFrom dendextend color_labels
#' @importFrom dendextend set
#' @importFrom stringr str_split
#' @importFrom dendextend color_branches
#' @importFrom stats hclust
#' @importFrom stats cutree
#' @importFrom RColorBrewer brewer.pal
#'
#' @return a list of
#' 1. the gene coexpression space dataframe,
#' 2. the eigenvalues plot (using eigenvalue from factoextra package)
#' 3. the pca component dataframe
#' 4. the tree plot for the coexpression space genes
#' @export
#'
#' @examples
setGeneric("get.gene.clusters", function(obj,list.group.markers,n.markers =25, k.cuts = 6,distance = "cosine",hclust.method="ward.D2")
    standardGeneric("get.gene.clusters"))
#' @rdname get.gene.clusters
setMethod("get.gene.clusters","scCOTAN",
 function(obj,list.group.markers,n.markers =25, k.cuts = 6,distance = "cosine",hclust.method="ward.D2"){

    g.space <- get.gene.coexpression.space(obj,
                                          n.genes.for.marker = n.markers,
                                          primary.markers = unlist(list.group.markers))
    g.space <- as.data.frame(g.space)

    coex.pca.genes <- prcomp(g.space,
                             center = TRUE,
                             scale. = F)
    #coex.pca.genes <- prcomp(t(g.space),
     #                        center = TRUE,
      #                       scale. = F)

    plot.eig <- fviz_eig(coex.pca.genes, addlabels=TRUE,ncp = 10)

    if (distance == "cosine") {
        hc.norm <- hclust(cosine.dissimilarity(as.matrix(t(g.space))),method =  hclust.method )
    }else if(distance == "euclidean"){
        hc.norm <- hclust(dist(g.space), method = hclust.method)
    }
    dend <- as.dendrogram(hc.norm)

    #pca_1 <- as.data.frame(coex.pca.genes$rotation[,1:10])
    pca_1 <- as.data.frame(coex.pca.genes$x[,1:10])
    pca_1 <- pca_1[order.dendrogram(dend),]

    cut <- stats::cutree(hc.norm, k = k.cuts)

    tmp <- calculatePValue(object = obj,
                           geneSubsetCol = unlist(list.group.markers),
                           geneSubsetRow = colnames(g.space))

    for (m in unlist(list.group.markers)) {
        tmp <- as.data.frame(tmp[order(tmp[,m]),])
        tmp$rank <- c(1:nrow(tmp))
        colnames(tmp)[ncol(tmp)] <- paste("rank", m, sep = ".")
    }
    rank.genes <- tmp[,(length(unlist(list.group.markers))+1):ncol(tmp)]
    #for (c in c(1:length(colnames(rank.genes)))) {
        colnames(rank.genes) <- stringr::str_split(colnames(rank.genes), pattern ="[.]",simplify = T)[,2]
    #}

    df.secondary.markers <- matrix(nrow = dim(rank.genes)[1],ncol = length(list.group.markers))
    rownames(df.secondary.markers) <- rownames(rank.genes)
    colnames(df.secondary.markers) <- names(list.group.markers)


    for (name in names(list.group.markers)) {
    #    print(name)
        df.secondary.markers[,name] <- rowSums(as.data.frame(rank.genes[,list.group.markers[[name]]]))
        df.secondary.markers[list.group.markers[[name]],name] <- 1
    }

    temp.coex = getGenesCoex(obj, genes = unlist(list.group.markers))
    temp.coex = temp.coex[rownames(temp.coex) %in% rownames(df.secondary.markers),]

    for (name in names(list.group.markers)) {
        to.change <- rowSums(as.data.frame(temp.coex[rownames(df.secondary.markers),list.group.markers[[name]]] < 0)) > 0
        to.change <- names(to.change[to.change > 0 ])
        df.secondary.markers[to.change,name] <- 100000
    }

    mylist.names <- colnames(df.secondary.markers)
    pos.link  <- set_names(vector("list", length(mylist.names)), mylist.names)
    for (g in rownames(df.secondary.markers)) {
        if(length( which(df.secondary.markers[g,] == min(df.secondary.markers[g,]))) == 1 ){
            pos.link[[which(df.secondary.markers[g,] == min(df.secondary.markers[g,])) ]] =
                c(pos.link[[which(df.secondary.markers[g,] == min(df.secondary.markers[g,])) ]], g)
        }
    }

    pca_1$highlight <- "not_marked"
    for (n in names(pos.link)) {
        pca_1[rownames(pca_1) %in% pos.link[[n]],]$highlight <- paste0("Genes related to ",n)
    }

    pca_1$hclust <- cut[rownames(pca_1)]
    if (k.cuts > 8) {
        col_vector <- RColorBrewer::brewer.pal(n = 8, name = "Set2")
        col_vector <- c(col_vector,RColorBrewer::brewer.pal(n = (k.cuts-8), name = "Set1"))
        col_vector <- col_vector[1:k.cuts]
    }else{
        col_vector = RColorBrewer::brewer.pal(n = k.cuts, name = "Set2")
    }

    pca_1$sec_markers <- 0
    pca_1[rownames(pca_1) %in% colnames(g.space),]$sec_markers <- 1

    pca_1$colors = "#B09C85FF"
    c = 1
    for (to.color in unique(pca_1$highlight)[!unique(pca_1$highlight) == "not_marked"]) {
        pca_1[pca_1$highlight == to.color, "colors"] = col_vector[c]
        c <- c+1
    }


    pca_1 <- pca_1[labels(dend),]
    col <- unique(pca_1$colors)
    pca_1$col_branches <- "#B09C85FF"
    pca_1$groupLabels <- pca_1$hclust
    for (c in unique(pca_1$hclust)) {
        for (m in unlist(list.group.markers)) {
            if (m %in% rownames(pca_1[pca_1$hclust == c,])) {
                pca_1[pca_1$hclust == c,"col_branches"] <- pca_1[m,"colors"]
                groupName <- names(which(sapply(list.group.markers, function(y) m %in% y)))
                pca_1[pca_1$hclust == c,"groupLabels"] <- groupName
            }
        }
    }

    dend <- dendextend::color_branches(dend, k = k.cuts,col = unique(pca_1[,c("hclust","col_branches")])[,2],
                           groupLabels = unique(pca_1[,c("hclust","groupLabels")])[,2])

    dend <- dendextend::color_labels(dend,labels = rownames(pca_1),col=pca_1$colors)

    return(list("g.space"=g.space, "plot.eig"=plot.eig,
                "pca_clusters"=pca_1,"tree_plot"=dend))
 }
)
