setGeneric("vec2mat_rfast", function(x,genes="all") standardGeneric("vec2mat_rfast"))
setMethod("vec2mat_rfast","numeric",
vec2mat_rfast <- function(x,genes="all"){
    if (genes == "all") {
        p <- sqrt(1 + 8 * length(x$coex.values))/ 2 - 0.5
        m <- matrix(0, p, p)
        m <- Rfast::lower_tri.assign(m, diag = TRUE,v = x$coex.values)
        m <- Rfast::upper_tri.assign(m,v = Rfast::upper_tri(Rfast::transpose(m)))
        #temp.names <- stringr::str_split(names(x),pattern = "[|]",simplify = T)
        rownames(m) <- x$genes
        colnames(m) <- x$genes
    }else{
        #all_genes = vector with gene names (it will create a coex matrix with the wanted
        # genes on columns and all genes in rows)
        m <- matrix(0, nrow = length(x$genes), ncol = length(genes))
        rownames(m) <- x$genes
        colnames(m) <- genes
        m
        for (pos.gene in match(genes, x$genes)) {
            temp.array <- x$coex.values[pos.gene]
            p =1
            l <- length(x$genes)
            s <- pos.gene
            while (p < (pos.gene-1)) {
                l <- l -1
                s <- s+l
                temp.array <- c(temp.array,x$coex.values[s])
                p <- p +1
                print(temp.array)
            }


            #linear part
            start.reading.position <- 0
            for (i in c(0:(pos.gene-2))) {
                start.reading.position <- start.reading.position + (length(x$genes)-i)
            }
            start.reading.position

            end.reading.position <- 0
            for (i in c(0:(pos.gene-1))) {
                end.reading.position <- end.reading.position + (length(x$genes)-i)
            }
            end.reading.position

            temp.array <- c(temp.array,x$coex.values[(start.reading.position+1):end.reading.position])
            m[,x$genes[pos.gene]] <- temp.array
        }


    }

    m
}

)

setGeneric("mat2vec_rfast", function(mat) standardGeneric("mat2vec_rfast"))
setMethod("mat2vec_rfast","numeric",
          mat2vec_rfast <- function(mat){
             v <- Rfast::lower_tri(mat,diag = T)
             names.v <- rownames(mat)
             v=list("genes"=rownames(mat),"coex.values"=v)
             v
          }
)

