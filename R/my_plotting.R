#' theme set with 4 axis line around, based on theme_bw
#'
#' @param legend.position legend position
#' @param baseSize base text size
#' @param legendSize legend size
#' @param xSize x lab size
#' @param ySize y lab size
#' @param x.angle x lab angle
#' @param x.vjust x lab vertical position
#'
#' @return a ggplot theme object
#' @export
#'
#' @examples my_themeFULL(); my_themeFULL(legend.position = 'top')
my_themeFULL = function(legend.position="right",
                        baseSize = 8, legendSize = 8, xSize = 10, ySize =10,
                        x.angle= 0, x.vjust = 0.5){
    # base size
    themeFULL = theme_bw(base_size = baseSize, base_family = "Arial") +
        # legend
        theme(legend.position=legend.position,
              panel.grid=element_blank(),
              legend.title = element_blank(),
              legend.text= element_text(color="black", size=legendSize),
              # title
              plot.title = element_text(hjust = 0.5),
              # axis
              axis.text.x = element_text(color="black", size= xSize, angle = x.angle, vjust = x.vjust), # x lab
              axis.text.y = element_text(color="black", size= ySize), # y lab
              axis.title.x = element_text(color="black"),
              axis.title.y = element_text(color="black"))

    return(themeFULL)
}



#' theme set with 2 axis line around ,bottom and left
#'
#' @param legend.position legend position
#' @param titleSize title text size
#' @param legendSize legend size
#' @param xlab.Size x lab size
#' @param ylab.Size y lab size
#' @param x.title.Size x axis title size
#' @param y.title.Size y axis title size
#' @param xlab.angle x lab angle
#' @param x.vjust x lab vertical position
#'
#' @return a ggplot theme object
#' @export
#'
#' @examples my_themeHALF(); my_themeHALF(legend.position = 'top')
my_themeHALF = function(legend.position = "bottom",
                        titleSize = 15, legendSize = 8,
                        xlab.Size = 10, ylab.Size =10,
                        x.title.Size = 12, y.title.Size =12,
                        xlab.angle= 0, x.vjust = 0.5){

    # font config for x,y lab, title
    axis_theme <- theme(plot.title = element_text(colour = "black", face = "bold", size = titleSize, hjust = 0.5), # hjust horizon position
                        axis.text.x = element_text(angle = xlab.angle, hjust = 0.5, vjust = x.vjust, color = "black", size=xlab.Size),
                        axis.text.y = element_text(color = "black",size = ylab.Size),
                        axis.title.x = element_text(color = "black",size = x.title.Size),
                        axis.title.y = element_text(color = "black",size = y.title.Size)
    )
    # background
    bg_theme <- theme(panel.background = element_rect(fill='transparent'),  # transparent background
                      panel.grid = element_blank(),                         # no gird
                      panel.border = element_blank(),                     # no border line
                      axis.line = element_line(colour = "black", size=0.8))   # axis line
    # lengend
    leg_theme <- theme(legend.text= element_text(color="black", size=legendSize),
                       legend.key = element_blank(),  # remove the gray background on legend fig
                       legend.title = element_blank(),
                       legend.position = legend.position )

    # return theme
    return(axis_theme + bg_theme + leg_theme)
}








#' vocano plot for DEG dataframe
#'
#' @param DEG.data  DEG dataframe (rownanmes genesymbol; colnames including log2FoldChange(need log2 transformed), pvalue raw value)
#' @param x column names meaning log2FoldChange
#' @param y column names meaning pvalue (raw value); transform for -log10P in the function
#' @param group factors of differentiation, including:'up','stable','down')
#' @param title plot title; deflaut blank
#' @param x.name plot x lab; deflaut log2 (Fold Change)
#' @param y.name plot y lab; deflaut -log10 (P-value)
#' @param pal color, default jco
#' @param plot_theme theme function returned object (e.g., my_themeFULL(), or theme_bw())
#' @param lab_text Gene names need to be labled
#' @param x_intercept  vline for lableing mininal log2foldchang3
#' @param y_intercept hline for lableing mininal -log10(pvalue)
#'
#' @return ggplot2 object
#' @export
#'
#' @examples  my_vocano1(DEG.data  = df, x = 'log2FoldChange', y = 'padj', group = 'Dif', x.name = 'log2(Foldchange)', y.name = '-log10(P.value)', lab_text = 'ASNS')
#'     my_vocano1(DEG.data  = df, x = 'log2FoldChange', y = 'padj', group = 'Dif',pal = c("#DC143C", "#808080","#00008B"), plot_theme = my_themeFULL())
my_vocano1 = function(DEG.data, x, y, group,
                      title='', x.name='log2 (Fold Change)', y.name='-log10 (P-value)',
                      lab_text='', x_intercept = 1, y_intercept = -log10(0.05),
                      pal = ggsci::pal_jco(alpha = 0.7)(9),
                      plot_theme ){
    ## DEG data processing
    Symbol = rownames(DEG.data)
    logFC = DEG.data[[x]]
    log10P = -log10(DEG.data[[y]])
    Dif = factor(DEG.data[[group]], levels = c('up','stable','down'))
    labs = Symbol
    labs[!(Symbol %in% lab_text)] = ''

    ## ploting data
    plot_df = data.frame(logFC, log10P, Dif, labs, stringsAsFactors=FALSE)

    ## ploting
    p =  ggplot(data=plot_df, aes(x=logFC, y=log10P, colour=Dif, fill=Dif)) +
        scale_color_manual(values = pal[1:3]) +
        geom_point(alpha=0.4, size=1) +
        ggrepel::geom_text_repel(aes(logFC, log10P, label= labs), size=3,
                                 show.legend  = F) + # remove a in the legend
        ggtitle(title) +  xlab(x.name) +  ylab(y.name)

    p = p + geom_vline(xintercept=c(-x_intercept,  x_intercept),lty=4,col="black",lwd=0.4) +
        geom_hline(yintercept =  y_intercept,lty=4,col="black",lwd=0.4)


    # adding theme if given
    if (!missing(plot_theme)){
        p = p + plot_theme
    }else{
        p = p + myRNAseq::my_themeHALF()
        }

    return(p)
    }





#' barplot for long dataframe
#'
#' @param df dataframe of long format
#' @param x column name of x axis data
#' @param y column name of y axis data
#' @param group column name of grouping; presenting bar color
#' @param title title, default blank
#' @param x.name x lab name, default blank
#' @param y.name y lab name, default blank
#' @param text_adding if need adding number of the colnum; default true
#' @param text_size size of the text adding
#' @param pal color, default jco; if given notice the number of factors
#' @param plot_theme theme function returned object (e.g., my_themeFULL(), or theme_bw())
#' @param x_level x lab levels, if x presenting as given order
#' @param group_level group lab levels, if group presenting as given order (show in legend figure)
#'
#' @return ggplot2 barplot
#' @export
#'
#' @examples my_barplot_long(df, x, y, group)
my_barplot_long = function(df, x, y, group, title='', x.name='', y.name='',
                           x_level, group_level,
                           text_adding = T, text_size = 5,
                           pal = ggsci::pal_jco(alpha = 0.7)(9),
                           plot_theme ){
    ## plot data
    df = as.data.frame(df)
    x_value = df[[x]]
    y_value = df[[y]]
    fillgroup = df[[group]]

    ## factor level given
    if (!missing(x_level)){x_value = factor(x_value, levels = x_level)}
    if (!missing(group_level)){fillgroup = factor(fillgroup, levels = group_level)}

    plot_df = data.frame(x_value,  y_value, fillgroup, stringsAsFactors=FALSE)


    ## ploting
    p = ggplot(data = plot_df, aes(x=x_value, y=y_value, group = fillgroup)) +
        geom_col(aes(fill=fillgroup), colour = 'black', position="dodge",width = 0.6) +
        theme(text=element_text(family="STHeitiSC-Medium"))

    p = p + scale_fill_manual(values = pal) +
        ggtitle(title) +  xlab(x.name) +  ylab(y.name)

    # adding text
    if(text_adding == T){
        p = p + geom_text(aes(label = y_value, y = y_value + max(y_value)*0.01),
                          colour = '#B0228C',position=position_dodge(0.6), size = text_size, vjust = 0)
    }

    # adding theme if given
    if (!missing(plot_theme)) {
        p = p + plot_theme
    }else{
        p = p + myRNAseq::my_themeFULL()
        }

    return(p)
}






#' boxplot for long dataframe
#'
#' @param df dataframe of long format
#' @param x column name of x axis data
#' @param y column name of y axis data
#' @param group column name of grouping; presenting bar color
#' @param title title, default blank
#' @param x.name x lab name, default blank
#' @param y.name y lab name, default blank
#' @param x_level x lab levels, if x presenting as given order
#' @param group_level group lab levels, if group presenting as given order (show in legend figure)
#' @param jitter.point if jitter points added; default true.
#' @param jitter.size size of jitter point
#' @param pal color, default jco; if given notice the number of factors
#' @param plot_theme theme function returned object (e.g., my_themeFULL(), or theme_bw())
#' @param pair_lab column name of the paried information
#'
#' @return ggplot2 box plot object
#' @export
#'
#' @examples my_boxplot_long(df, x, y, group)
my_boxplot_long = function(df, x, y, group, title='', x.name='', y.name='',
                           x_level, group_level,
                           pair_lab,
                           jitter.point = T, jitter.size = 1,
                           pal = ggsci::pal_jco(alpha = 0.7)(9),
                           plot_theme ){
    ## plot data
    df = as.data.frame(df)
    x_value = df[[x]]
    y_value = df[[y]]
    fillgroup = df[[group]]

    ## factor level given
    if (!missing(x_level)){x_value = factor(x_value, levels = x_level)}
    if (!missing(group_level)){fillgroup = factor(fillgroup, levels = group_level)}

    ## if paried data
    if(missing(pair_lab)){
        plot_df = data.frame(x_value,  y_value, fillgroup, stringsAsFactors=FALSE)
    }else{
        pair = df[[pair_lab]]
        plot_df = data.frame(x_value,  y_value, fillgroup, pair, stringsAsFactors=FALSE)
    }


    ## ploting
    p = ggplot(plot_df, aes(x=x_value, y=y_value, fill=fillgroup)) +
        geom_boxplot()

    p = p + scale_fill_manual(values = pal) +
        ggtitle(title) +  xlab(x.name) +  ylab(y.name)

    ## if pari line added
    if(!missing(pair_lab)){
        p = p + geom_point(alpha=0.6, size=1.5) +
            geom_line(aes(group = pair), color = 'gray', lwd = 0.5)
    }

    # adding jitter point
    if(jitter.point == T & missing(pair_lab)){
        p = p + geom_point(alpha=0.6, size=jitter.size, position=position_jitterdodge())
    }

    # adding theme if given
    if (!missing(plot_theme)) {
        p = p + plot_theme
    }else{
        p = p + myRNAseq::my_themeFULL()
        }

    return(p)
}





