
#' GSEA analysis for human gene symbols
#'
#' @param DEG DEG dataframe from RNAseq results such as Deseq2
#' @param fc_column column name of fold change in the DEG
#' @param gene_column column name of genesymbol in the DEG
#' @param GMTfile gmt object of picked geneset gmt file
#'
#' @return GSEA object
#' @export
#'
#' @examples my_GSEA_analysis(DEG, 'log2FoldChange', 'geneID', GMT.hallmark)
my_GSEA_analysis <- function(DEG, fc_column, gene_column, GMTfile){

    ## prepare DEG dataframe
    DEG_df = data.frame(log2FC = DEG[[fc_column]],
                        gene_symbol = DEG[[gene_column]])
    ## prepare genelist
    DEG_df  = DEG_df  %>% arrange(desc(log2FC)) %>%
        distinct(gene_symbol, .keep_all = T)

    geneList = DEG_df[['log2FC']]
    names(geneList) = DEG_df[['gene_symbol']]
    geneList = sort(geneList,decreasing = T)

    ## GSEA
    egmt = GSEA(geneList = geneList,
                TERM2GENE = GMTfile,
                verbose = F,
                pvalueCutoff = 1)

    return(egmt)
}




#' GSEA plot a simple version; given one pathway name interested
#'
#' @param pathway pathway name , as shown in the gsea analysis
#' @param GSEA.res gsea result object, not a data frame!
#' @param base_size font size
#'
#' @return  a simple version of gsea plot of a single pathway
#' @export
#'
#' @examples GSEA.plot_simple(pathway, GSEA.res)
GSEA.plot_simple <- function(pathway, GSEA.res,  base_size = 4){
    p = enrichplot::gseaplot2(x = GSEA.res,
                              color = "red",
                              geneSetID = pathway,
                              pvalue_table = F,
                              base_size = base_size,
                              title = str_split(pathway, pattern = '_', n = 2, simplify = T)[2])

    return(p)
}





#' Optimized GSEA plot version; given one pathway name interested
#'
#' @param path pathway name , as shown in the gsea analysis
#' @param GSEA_Object gsea result object, not a data frame!
#' @param leading_genes leading edge gene symbols
#'
#' @return  a simple version of gsea plot of a single pathway
#' @export
#'
#' @examples my_GSEAPlot(path = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", GSEA_Object = gseaRes)
my_GSEAPlot <- function(path, GSEA_Object, leading_genes){

    if(!missing(leading_genes)){
        gseaNb(object = GSEA_Object,
               geneSetID = path,
               subPlot = 2,
               addPval = T, pvalX = 0.6, pvalY = 0.8,
               pCol = 'black',
               pHjust = 0,
               addGene = leading_genes,
               termWidth = 30)
    }else{
        gseaNb(object = GSEA_Object,
               geneSetID = path,
               subPlot = 2,
               addPval = T, pvalX = 0.6, pvalY = 0.8,
               pCol = 'black',
               pHjust = 0,
               termWidth = 30)
    }

}




#' GSEA presented as barplot ; given top pathways arranged by pvalue
#'
#' @param GESA.res gsea result object, not a data frame!
#' @param top  number of top pathways picked, default 30
#' @param pal color palette represent for up or down regulation; down as the first color; up as the second color ;default in 'lanonc' pal
#' @param title title of the plot , default GSEA
#' @param legend.position position of the legend, default right
#'
#' @return bar plot of GSEA result
#' @export
#'
#' @examples my_GSEA_barplot1(GSEAres)
GSEA_barplot1 <- function(GESA.res, top=30,
                          pal=ggsci::pal_lancet(palette = c("lanonc"), alpha = 0.6)(9),
                          title = 'GSEA', legend.position = 'right' ){

    dat = GESA.res@result
    colnames(dat)

    dat = dat %>% mutate(pathID = str_split(Description, pattern = '_', n = 2, simplify = T)[ ,2]) %>%
        mutate(Description = str_replace_all(pathID, pattern = '_', replacement = ' '))

    ## color for up or down regulation
    dat$color = ifelse(dat$NES > 0, 'up', 'down')
    dat$color = factor(dat$color,  levels = c('up', 'down'))

    ## top 30 pathways, arranged by pvalue,than ranged by NES
    dat = dat %>% arrange(pvalue) %>% .[1:top, ] %>%
        arrange(NES)

    ## plotting
    p = ggplot(data = dat, mapping = aes(x = reorder(Description, order(NES, decreasing = F)),
                                         y = NES,
                                         fill = color)) +
        geom_col() +
        coord_flip() +
        scale_fill_manual(values = c('up'=pal[2], 'down'=pal[1]))

    p = p + ggtitle(title) + xlab("") +
        theme_bw() +
        theme(legend.position= legend.position,
            legend.title = element_blank(),
            axis.title.y = element_blank())


    return(p)

}





#' GSEA presented as barplot ; given top pathways arranged by pvalue, pathway name beside the bar
#'
#' @param GESA.res gsea result object, not a data frame!
#' @param top number of top pathways picked, default 30
#' @param pal color palette represent for up or down regulation; down as the first color; up as the second color ;default in 'lanonc' pal
#' @param text.size text size of pathway names,default 4
#' @param xlab.size text size of xlab,default 15
#' @param legendtext.size  text size of legend,default 12
#'
#' @return bar plot of GSEA result
#' @export
#'
#' @examples my_GSEA_barplot2(GSEAres)
GSEA_barplot2 <- function(GESA.res, top=30,
                          pal=ggsci::pal_lancet(palette = c("lanonc"), alpha = 0.6)(9),
                          text.size=4, xlab.size=15, legendtext.size=12){
    ## data
    dat = GESA.res@result
    # colnames(dat)

    dat = dat %>% mutate(pathID = str_split(Description, pattern = '_', n = 2, simplify = T)[ ,2]) %>%
        mutate(Description = str_replace_all(pathID, pattern = '_', replacement = ' '))

    ## color for up or down regulation
    dat$color = ifelse(dat$NES > 0, 'up', 'down')
    dat$color = factor(dat$color,  levels = c('up', 'down'))

    ## top 30 pathways, arranged by pvalue,than ranged by NES
    dat = dat %>% arrange(pvalue) %>% .[1:top, ] %>%
        arrange(NES)

    ## plotting
    p = ggplot(dat, mapping = aes(x = reorder(Description, order(NES, decreasing = F)), # works only data previously arranged
                                  y = NES)) +
        geom_col(aes(fill= color), colour = 'black', position="dodge") +
        scale_fill_manual(values = c('up'=pal[2], 'down'=pal[1])) +
        coord_flip() +
        theme_void(base_size = 15, base_family = "Arial", ) +
        ggtitle('GSEA') +  xlab('') +  ylab('NES\n') +
        theme(plot.title = element_text(colour = "black", face = "bold", hjust = 0.5),
              legend.position = "bottom",
              panel.grid=element_blank(),
              legend.title = element_blank(),
              legend.text=element_text(size=legendtext.size),
              axis.text.x = element_text(color="black", size =xlab.size),
              axis.title.x = element_text(color="black"),
              axis.line.x = element_line(colour = "black"),
              axis.ticks.x = element_line(colour = "black", size=1.5, # ticks set up didn't work
                                          linetype=1, lineend=1))

    ## adding pathway names # see depmap project file 01_
    p <- p + geom_text(data = subset(dat, NES < 0),
                       mapping = aes(label = Description, x = Description, y = 0.1, colour = color),
                       angle = 0, hjust = 0, size = text.size, show.legend  = F) +

             geom_text(data = subset(dat, NES >= 0),
                       mapping = aes(label = Description, x = Description, y =  -0.1, colour = color),
                       angle = 0, hjust = 1,  size = text.size, show.legend  = F) +
        scale_color_manual(values = c('up'=pal[2], 'down'=pal[1]))


    return(p)
}







#' GSEA presented as chicklet-barplot ; given top pathways arranged by pvalue, pathway name beside the bar
#'
#' @param GESA.res gsea result object, not a data frame!
#' @param top number of top pathways picked, default 30
#' @param pal color palette represent for up or down regulation; down as the first color; up as the second color ;default pal_simpsons
#' @param text.size text size of pathway names,default 3
#' @param xlab.size text size of xlab,default 15
#' @param legendtext.size text size of legend,default 12
#'
#' @return chicklet-bar plot of GSEA result
#' @export
#'
#' @examples GSEA_chickletBar(GESA.res, top=30)
GSEA_chickletBar <- function(GESA.res, top=30,
                          pal=ggsci::pal_simpsons()(2),
                          text.size=3, xlab.size=15, legendtext.size=12){
    ## data
    dat = GESA.res@result

    dat = dat %>% mutate(pathID = str_split(Description, pattern = '_', n = 2, simplify = T)[ ,2]) %>%
        mutate(Description = str_replace_all(pathID, pattern = '_', replacement = ' '))

    ## color for up or down regulation
    dat$color = ifelse(dat$NES > 0, 'up', 'down')
    dat$color = factor(dat$color,  levels = c('up', 'down'))

    ## top 30 pathways, arranged by pvalue,than ranged by NES
    dat = dat %>% arrange(pvalue) %>% .[1:top, ] %>%
        arrange(NES) %>% as.data.frame()
    # factor levels for lab order
    dat$pathID = factor(dat$pathID, levels = dat$pathID)
    # Max lab for border of the plot axis
    max.value = ceiling(max(abs(dat$NES)))

    ## plotting
    p = ggplot(dat, mapping = aes(x = reorder(pathID, order(NES, decreasing = F)), # works only data previously arranged
                                  y = NES)) +
        # base color
        ggchicklet::geom_chicklet(aes(y = ifelse(NES>0, max.value, -max.value),
                          fill= color), alpha = 0.2) +
        scale_fill_manual(values = c('up'=pal[2], 'down'=pal[1])) +
        # bar color
        ggchicklet::geom_chicklet(aes(fill= color), colour = 'black', alpha=0.8) +
        scale_fill_manual(values = c('up'=pal[2], 'down'=pal[1])) +
        # text , decided by NES
        # a start position 0.1 for negative text,; -0.1 for positvie text
        # allign for hjust : 0 , left allign for negative; 1 right allign for positive text
        geom_text(aes(x = pathID,
                      y = ifelse(NES > 0,  -0.1, 0.1),
                      label = pathID,
                      colour = color,
                      hjust= ifelse(NES>0, 1, 0)),
                  size = text.size,
                  show.legend  = F) +  # remove a in the legend
        scale_color_manual(values = c('up'=pal[2], 'down'=pal[1])) +

        # theme
        theme_classic(base_size = xlab.size) +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
              legend.position = "bottom", legend.text=element_text(size=legendtext.size), legend.title = element_blank()) +
        labs(x = NULL, y='NES') +
        coord_flip()

    return(p)
}


