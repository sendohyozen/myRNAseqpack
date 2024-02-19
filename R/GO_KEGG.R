

#' GO analysis for human gene symbols; human only!
#'
#' @param DEG.genes gene symbols of DEGs
#' @param background.genes background gene symbols ,default not provide
#'
#' @return a list consists of GO result objects, inculuding BP, MF, CC
#' @export
#'
#' @examples my_GO(deg.sig)
my_GO <- function(DEG.genes, background.genes){

    ## gene symobl to entrnzID
    df <- bitr(geneID = unique(DEG.genes),
               fromType = "SYMBOL",
               toType = c( "ENTREZID"),
               OrgDb = org.Hs.eg.db)
    gene = df[['ENTREZID']]

    ## GO analysis
    if(missing(background.genes)){

        ## defalut, background genes not provided
        go_enrich_results <- lapply( c('BP','MF','CC') , function(ont) {
            cat(paste('Now process',ont, '\n', collapse = ''))
            ego <- enrichGO(gene          = gene,
                            OrgDb         = org.Hs.eg.db,
                            ont           = ont ,
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.9,
                            qvalueCutoff  = 0.9,
                            readable      = TRUE)
            return(ego)
        })
    }else{

        ## background genes provided
        df_bg <- bitr(geneID = unique(background.genes),
                   fromType = "SYMBOL",
                   toType = c( "ENTREZID"),
                   OrgDb = org.Hs.eg.db)

        genes_background = df_bg[['ENTREZID']]

        go_enrich_results <- lapply( c('BP','MF','CC') , function(ont) {
            cat(paste('Now process',ont, '\n', collapse = ''))
            ego <- enrichGO(gene          = gene,
                            universe      = genes_background ,
                            OrgDb         = org.Hs.eg.db,
                            ont           = ont ,
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.9,
                            qvalueCutoff  = 0.9,
                            readable      = TRUE)
            return(ego)
        })
    }

    return(go_enrich_results)

}






#' KEGG analysis for human gene symbols; human only!
#'
#' @param DEG.genes gene symbols of DEGs
#'
#' @return kegg result object
#' @export
#'
#' @examples my_kegg(deg)
my_kegg <- function(DEG.genes){

    ## gene symobl to entrnzID
    df <- bitr(geneID = unique(DEG.genes),
               fromType = "SYMBOL",
               toType = c( "ENTREZID"),
               OrgDb = org.Hs.eg.db)
    gene = df[['ENTREZID']]

    ## KEGG analysis
    kegg_result <- enrichKEGG(gene = gene,
                        organism = 'hsa',
                        keyType = 'kegg',
                        pvalueCutoff = 0.9,
                        pAdjustMethod = "BH")


    return(kegg_result)
}






