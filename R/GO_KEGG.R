

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






#' barplot (ordered by count/ratio) and group by category of ontology (BP, MF, CC)
#'
#' @param GO.res GO result object  generated from my_GO() function including (BP, MF, CC)
#' @param top number of top pathways picked, default 10
#' @param y_value ordered by count or ratio, default "count", if need ratio, use  y_value = 'ratio'
#' @param decreasing.Y ordered by decreasing,default T
#' @param color.text if using colored labels, (go terms generated by geom_text, better give y.limit value) / if not show go terms on axis text
#' @param lab.size text size of colored labels, default 4
#' @param y.limit a vector with 2 numbers as the ymin and ymax: e.g y.limit = c(-30, 40), if not using default value
#' @param pal color palette represent for up or down regulation; down as the first color; up as the second color ;default in 'lanonc' pal
#' @param bar.width bar width ,default 0.6
#' @param title plot title ,default NULL
#' @param xlab plot title ,default NULL
#' @param ylab plot title ,default count / ratio
#' @param base.size base text size, default 15
#' @param legend.size base legend size, default 10
#' @param legend.position legend positon, default top
#' @param cord.flip if flip the plot default F
#'
#' @return a sorted ggplot bar plot grouped by GO terms
#' @export
#'
#' @examples GO_barplot_catergroy(GO.res = GOResults_list[[1]], color.text = T, y.limit = c(-30, 45))
#'
GO_barplot_catergroy <- function(GO.res, top=10, y_value="count", decreasing.Y = TRUE,
                                 color.text=F, lab.size = 4, y.limit,
                                 pal=ggsci::pal_lancet(palette = c("lanonc"), alpha = 0.6)(9),
                                 bar.width = 0.6, title =NULL, xlab =NULL, ylab = "Gene counts",
                                 base.size =15, legend.size=10, legend.position = 'top',
                                 cord.flip=F){

    BP = GO.res[[1]]@result[1:top, ]; BP$ontology = 'BP'
    MF = GO.res[[2]]@result[1:top, ]; MF$ontology = 'MF'
    CC = GO.res[[3]]@result[1:top, ]; CC$ontology = 'CC'

    data = rbind(BP, MF, CC) %>% dplyr::select(Description, Count, ontology, GeneRatio) %>%
        mutate(ontology = factor(ontology, levels=c("BP", "MF", "CC")) ) %>% as.data.frame()

    # arrange by count
    plot_data = data %>%
        arrange(ontology, desc(Count)) %>%
        mutate(Description = factor(.$Description, levels = .$Description) )  %>%
        dplyr::rename(y = Count)  %>%
        as.data.frame()


    if(y_value=="ratio"){
        data$bg_count = stringr::str_split(data$GeneRatio, pattern = '/', simplify = T)[,2] %>% as.numeric()

        plot_data = data %>%
            mutate(ratio = as.numeric(Count/bg_count)*100) %>%
            arrange(ontology, desc(ratio)) %>%
            mutate(Description = factor(.$Description, levels = .$Description) )  %>%
            dplyr::rename(y = ratio)  %>%
            as.data.frame()

        ylab = 'Gene ratio (%)'
    }

    #  setting color for category
    category_colors <- c('BP'=pal[1], 'MF'=pal[2], 'CC'= pal[3])


    # setting Y axis limits
    if(color.text==T & missing(y.limit)){
        ymin = -30
        ymax = max(plot_data$y) + 1
    }else if(color.text==F){
        cat('default y limits\n')
    }else{
        cat('y.limit need a vector with 2 numbers as the ymin and ymax: e.g y.limit = c(-30, 40)')
        ymin = y.limit[1]
        ymax = y.limit[2]
    }

    ## plotting
    p = ggplot(plot_data, mapping = aes(x = reorder(Description, order(ontology, y, decreasing = c(FALSE, decreasing.Y))), # works only data previously arranged
                                        y = y)) +
        geom_bar(aes(fill= ontology), stat = "identity", colour = 'black', position="dodge", width = bar.width) +
        scale_fill_manual(values = category_colors )

    # title
    p = p + ggtitle(title) +  xlab(xlab) +  ylab(ylab)

    # adding text as x lab
    if(color.text==T){

        p = p +  scale_y_continuous(limits=c(ymin, ymax))

        if(cord.flip==T){
            p = p + coord_flip() +
                geom_text(aes(label = Description, color = ontology, y = -0.5), vjust = 0.5, hjust = 1,  size = lab.size, angle = 0,
                          # position = position_dodge(width = bar.width)
                          show.legend = F) +
                scale_color_manual(values = category_colors)

            # theme
            p = p + theme_classic(base_size = base.size) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
                                                                 axis.ticks.x = element_blank(),
                                                                 legend.position = legend.position, legend.title = element_blank(), legend.text = element_text(size = legend.size))

        }else{
            p = p +  geom_text(aes(label = Description, color = ontology, y = -0.5), vjust = 0.5, hjust = 1,  size = lab.size, angle = 90,
                               # position = position_dodge(width = bar.width)
                               show.legend = F) +
                scale_color_manual(values = category_colors)

            # theme
            p = p + theme_classic(base_size = base.size) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
                                                                 axis.ticks.y = element_blank(),
                                                                 legend.position = legend.position, legend.title = element_blank(), legend.text = element_text(size = legend.size))
        }


    }else{

        if(cord.flip==T){
            p = p + coord_flip() +
                theme_classic(base_size = base.size) + theme(axis.text.x = element_text(colour = 'black', angle = 0, vjust = 0.5, hjust = 0.5), axis.ticks.x = element_blank(),
                                                             axis.ticks.y = element_blank(),
                                                             legend.position = legend.position, legend.title = element_blank(), legend.text = element_text(size = legend.size))
        }else{
            # if no color lab needed use default x axis instead
            p = p + theme_classic(base_size = base.size) + theme(axis.text.x = element_text(colour = 'black', angle = 90, vjust = 0.5, hjust = 1), axis.ticks.x = element_blank(),
                                                                 axis.ticks.y = element_blank(),
                                                                 legend.position = legend.position, legend.title = element_blank(), legend.text = element_text(size = legend.size))
        }
    }


    return(p)
}

