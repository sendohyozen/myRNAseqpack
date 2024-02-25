

#' GSVA analysis
#'
#' @param exprSet A normlized expression matrix (e.g. TPM or FPKM, not raw count); genes are rows, columns are samples
#' @param GMTfile gmt object of picked geneset gmt file
#'
#' @return  GSVA transformed geneset scored matrix (pathways are rows, columns are samples)
#' @export
#'
#' @examples my_GSVA_analysis(exprSet, GMTfile)
my_GSVA_analysis <- function(exprSet, GMTfile){

    ## prepareing a genelist from gmt object
    gs <- lapply(GMTfile@.Data, function(x){
        gene.list =  x@geneIds
    })

    names(gs) = unlist(lapply(GMTfile@.Data, function(x){
        x@setName
    }))

    cat('gmt geneset transformed successfully!')

    ## GSVA analysis
    GSVA_ES <- gsva(expr=exprSet,
                    gset.idx.list=gs,
                    parallel.sz=5)


}





#' limma result of DEG from GSVA scoring matrix  presented as barplot
#'
#' @param gsvaScore.limma.DEG limma result of DEG generated from GSVA scoring matrix
#' @param logFC log2(fold change) , deflaut 0.2
#' @param pV pvalue deflaut 0.05
#' @param addline if adding a helping line
#' @param axis_need_color if adding color on the axis lab of pathway name; showing up or down
#' @param palettes color palettes as down, stable, up
#'
#' @return bar plot of DEG using limma from a GSVA score matrix
#' @export
#'
#' @examples GSVA_barplot1(gsvaScore.limma.DEG)
GSVA_barplot1 <- function(gsvaScore.limma.DEG, logFC=0.2, pV=0.05,
                        addline = F,  axis_need_color = T,
                        palettes = c("#008020","#808080", "#08519C") ){

    ## use the function in the same package, :: needed or not?
    df = DEGsig(df = gsvaScore.limma.DEG,
                logFC = logFC, pV = pV,
                method='limma')
    ## simplified the pathway name in gmt file
    ID = str_split(rownames(df), pattern = '_', n = 2, simplify = T)[ ,2]
    ID = str_replace(ID, pattern = '_', replacement = '')

    ## data for plotting
    df = data.frame(ID = ID,
                    FC = df[['log2FoldChange']],
                    group = df[['Dif']])

    df.sort = df[order(df$FC), ]
    df.sort$ID = factor(df.sort$ID, levels = df.sort$ID)  # sorting as reordered names

    ## ploting
    if(axis_need_color==T){
        color_axis = ifelse(df.sort$group == 'down', palettes[1],
                            ifelse(df.sort$group == 'up', palettes[3], palettes[2]))
    }else{
        color_axis = "black"
    }


    p = ggplot(df.sort, aes(x = ID, y = FC, fill = group)) +
        geom_bar(stat = 'identity',alpha = 0.7) +
        theme_bw() +
        theme(legend.position="top",
              legend.title = element_blank(),
              panel.grid =element_blank(),
              legend.text= element_text(color="black", size=12),
              plot.title = element_text(hjust = 0.5, size=12),
              axis.text.x = element_text(color="black", size=12),
              axis.text.y = element_text(color=color_axis, size=8),
              axis.title.x = element_text(color="black", size=12),
              axis.title.y = element_text(color="black", size=12)) +
        theme(panel.border = element_rect(size = 0.6)) +
        ggtitle('GSVA') +  xlab('log2FC') +  ylab('') +
        scale_fill_manual(values = palettes) +
        coord_flip()


    ## adding helping line
    if(addline == T){
        p = p + geom_hline(yintercept = c(-logFC, logFC) ,lty=4,col="black", lwd=0.4)
    }

    return(p)
}




#' limma result of DEG from GSVA scoring matrix  presented as barplot, pathway name beside the bar
#'
#' @param gsvaScore.limma.DEG limma result of DEG generated from GSVA scoring matrix
#' @param top number of top pathways picked, default 20
#' @param pal color palette represent for up or down regulation; down as the first color; up as the second color ;default in 'lanonc' pal
#' @param text.size text size of pathway names,default 4
#' @param xlab.size text size of xlab,default 15
#' @param legendtext.size text size of legend,default 12
#'
#' @return bar plot of GSVA result
#' @export
#'
#' @examples GSVA_barplot2(gsvaScore.limma.DEG)
GSVA_barplot2 <- function(gsvaScore.limma.DEG, top=20,
                          pal=ggsci::pal_lancet(palette = c("lanonc"), alpha = 0.6)(9),
                          text.size=4, xlab.size=15, legendtext.size=12){
    ## data
    df = gsvaScore.limma.DEG %>% mutate(Dif= ifelse(logFC < 0, 'down', 'up'))
    ## simplified the pathway name in gmt file
    ID = str_split(rownames(df), pattern = '_', n = 2, simplify = T)[ ,2]
    ID = str_replace(ID, pattern = '_', replacement = '')

    ## data for plotting,colunames from limma
    dat = data.frame(ID = ID,
                    FC = df[['logFC']],
                    pvalue = df[['adj.P.Val']],
                    color = df[['Dif']])


    ## color for up or down regulation
    dat$color = factor(dat$color,  levels = c('up', 'down'))

    ## top 30 pathways, arranged by pvalue,than ranged by NES
    dat = dat %>% arrange(pvalue) %>% .[1:top, ] %>%
        arrange(FC)

    ## plotting
    p = ggplot(dat, mapping = aes(x = reorder(ID, order(FC, decreasing = F)), # works only data previously arranged
                                  y = FC)) +
        geom_col(aes(fill= color), colour = 'black', position="dodge") +
        scale_fill_manual(values = c('up'=pal[2], 'down'=pal[1])) +
        coord_flip() +
        theme_void(base_size = 15, base_family = "Arial", ) +
        ggtitle('GSVA') +  xlab('') +  ylab('log2FoldChange\n') +
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
    p <- p + geom_text(data = subset(dat, FC < 0),
                       mapping = aes(label = ID, x = ID, y = 0.05, colour = color),
                       angle = 0, hjust = "outward", size = text.size, show.legend  = F) +

        geom_text(data = subset(dat, FC >= 0),
                  mapping = aes(label = ID, x = ID, y =  -0.05, colour = color),
                  angle = 0, hjust = "outward",  size = text.size, show.legend  = F) +
        scale_color_manual(values = c('up'=pal[2], 'down'=pal[1]))


    return(p)
}








#' limma result of DEG from GSVA scoring matrix  presented as vacano plot
#'
#' @param gsvaScore.limma.DEG limma result of DEG generated from GSVA scoring matrix
#' @param logFC log2(fold change) , deflaut 0.2
#' @param pV pvalue deflaut 0.05
#' @param labeled.pathway pathway name to be labeled
#'
#' @return volcano plot of DEG using limma from a GSVA score matrix
#' @export
#'
#' @examples GSVA_volcano(gsvaScore.limma.DEG)
GSVA_volcano = function(gsvaScore.limma.DEG, logFC=0.2, pV=0.05,
                        pal  = c("#DC143C", "#808080","#00008B" ),
                        title='GSVA', x.name='log2 (Fold Change)', y.name='-log10 (P-value)',
                        labeled.pathway){

    ## use the function in the same package, :: needed or not?
    deg_df = DEGsig(df = gsvaScore.limma.DEG,
                logFC = logFC, pV = pV,
                method='limma')

    ## plotting data
    df = data.frame(Symbol = rownames(deg_df),
                    logFC = deg_df[['log2FoldChange']],
                    log10P = -log10(deg_df[['padj']]),
                    Dif = factor(deg_df[['Dif']], levels = c('up', 'stable', 'down')),
                    row.names = rownames(deg_df),
                    stringsAsFactors=FALSE)


    ## Adding lab text，ggrepl consitant with the dataframe
    df$labs = df$Symbol
    ## don't labeled ones give a '' value
    if(missing(labeled.pathway)){
        df$labs = ''
    }else{
        df$labs[!(df$labs %in% labeled.pathway)] = ''
    }


    ## plotting
    p = ggplot(data=df, aes(x=logFC, y=log10P, colour=Dif, fill=Dif)) +
        scale_color_manual(values = pal[1:3]) +
        geom_point(alpha=0.4, size=1) +
        ggrepel::geom_text_repel(aes(logFC, log10P, label= labs), size=3,
                                 show.legend  = F) +  # remove a in the legend
        ggtitle(title) +  xlab(x.name) +  ylab(y.name)

    p = p + geom_vline(xintercept=c(-logFC,  logFC),lty=4,col="black",lwd=0.4) +
        geom_hline(yintercept =  -log10(pV),lty=4,col="black",lwd=0.4) +
        theme_bw()

    return(p)
}

