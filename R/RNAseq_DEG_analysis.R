

#' pca visulization for RNAseq grouping
#'
#' @param dat raw count matrix (genes are rows, columns are samples)
#' @param group group of samples; consistent with samples as column of the matrix
#' @param pal color of the groups;number of the colors should be consistent with number of levels of the group
#'
#' @return pca plot
#' @export
#'
#' @examples my_pca(dat, group)
my_pca <- function(dat, group,
                   pal=ggsci::pal_jco(alpha = 0.7)(9)) {

    dat=as.data.frame(t(dat))
    group = as.factor(group)

    ## before PCA analysis
    dat.pca <- FactoMineR::PCA(dat, graph = FALSE)
    ## pca visual
    pca.plot = factoextra::fviz_pca_ind(dat.pca,
                            geom.ind = "point",
                            col.ind = group,
                            palette = pal[1:nlevels(group)],
                            addEllipses = TRUE,
                            legend.title = "Groups")
    return(pca.plot)
}








#' clusting visulization for RNAseq grouping
#'
#' @param matrix raw count matrix (genes are rows, columns are samples)
#'
#' @return hclust plot
#' @export
#'
#' @examples my_hclust(exprSet)
my_hclust <- function(matrix) {
    # Define nodePar
    nodePar <- list(lab.cex = 0.6, pch = c(NA, 19),
                    cex = 0.7, col = "blue")
    hc=hclust(dist(t(matrix)))
    clust.plot = plot(as.dendrogram(hc), nodePar = nodePar,  horiz = TRUE)

    return(clust.plot)
}











#' Filtering a count matrix to remove low read genes
#' @description This function takes a count matrix and removes genes that don't reach a minimum namer of reads in a certain number of samples.
#'
#' @param mmatrix raw count matrix (genes are rows, columns are samples)
#' @param threshold minimum persentage of samples in which CPM threshold has to be reached; defaults to ten percent (0.1)
#' @param minCpm the minimum counts per million (CPM) that has to be reached; defaults to ten percent (0.25)
#'
#' @return filtered count matrix
#' @export
#'
#' @examples filter_count_matrix(matrix, threshold = 0.1, minCpm=0.25)
filter_count_matrix <- function(matrix, threshold = 0.1, minCpm=0.25) {
    tmp <- edgeR::cpm(matrix)
    sel <- apply(tmp, 1, function(x) {
        sum(x > minCpm) >= ncol(matrix)*threshold
    })
    return(matrix[sel,])
}







#' Adding differatiation column
#'
#' @param df  DEG data frame analysed with deseq2 limma or edgeR
#' @param logFC log2(fold change)
#' @param pV pvalue
#' @param method Deseq2 limma or edgeR (default Deseq2)
#'
#' @return DEG dataframe adding a column named Dif
#' @export
#'
#' @examples DEGsig(df, logFC =1, pV = 0.05, method = 'Deseq2')
DEGsig <- function(df, logFC, pV, method='Deseq2'){

    if(method =='Deseq2'){
        df = as.data.frame(df)

    }else if(method=='limma'){
        df = df %>% dplyr::rename(log2FoldChange=logFC, padj=`adj.P.Val`) %>% as.data.frame()

    }else if(method=='edgeR'){
        df = df %>% dplyr::rename(log2FoldChange=`log2(fc)`, padj=pval) %>% as.data.frame()

    }else{
        stop( "Only suitable for Deseq2, limma, edgeR!" )
    }

    ## increase a column named Dif, indicating up ,down or stable based on logFC and pV
    df$Dif=ifelse(df$padj > pV, 'stable',
                  ifelse( df$log2FoldChange > logFC, 'up',
                          ifelse( df$log2FoldChange < -logFC,'down','stable') )
    )
    return(df)

}





