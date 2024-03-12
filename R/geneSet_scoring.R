
#' Calculate score across genes and samples using PCA for a matrix with genes in a defined geneset
#'
#' @param gm normalized count matrix; rows are all genes in the signature that shall be summarized into one score; columns are samples
#' @param summarizationFunction character vector defining whether signature scores shall be based on principal component 1 ("PCA", default) or z-scores (other value)
#'
#' @return numeric vector of length ncol(gm); a score summarizing the rows of gm
#' @export
gsScore <- function(gm, summarizationFunction="PCA") {
    if (summarizationFunction == "PCA") {
        pc <- prcomp(t(gm),
                     retx=TRUE)
        gss <- pc$x[,1] * sign(cor(pc$x[,1], colMeans(gm)))
    } else {
        gss <- colMeans(gm)
    }

    return(gss)
}
