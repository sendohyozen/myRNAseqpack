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
                      plot_theme){
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
    if (!missing(plot_theme)) {
        p = p + plot_theme}

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
                           plot_theme){
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
        p = p + plot_theme}

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
                           plot_theme){
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
        p = p + plot_theme}

    return(p)
}





#' bubleplot for go,kegg analysis
#'
#' @param df dataframe of go, kegg results
#' @param x column name of x axis
#' @param y column name of y axis
#' @param size column name showing size of the buble
#' @param color column name showing  color of the buble
#' @param title title of the plot , default "KEGG"
#' @param pal color palette of the bulbe color; a vector consists of two numbers, the fisrt for low , the second for high
#'
#' @return ggplot2 object of a buble plt
#' @export
#'
#' @examples my_bublePlot(df = data, x = 'pvalue', y = 'Description', size = 'Count' , color = 'p.adjust', title = 'KEGG')
my_bublePlot = function(df, x, y, size, color,
                        title='KEGG', pal=c('blue', 'red')){

    ## plot data
    df = as.data.frame(df)
    x_value = df[[x]]
    y_value = df[[y]]
    buble_Size = df[[size]]
    buble_color = df[[color]]

    plot_df = data.frame(x_value,  y_value,  buble_Size, buble_color, stringsAsFactors=FALSE)

    ## ploting
    p = ggplot(data= plot_df, aes(x=x_value, y=y_value)) +
        geom_point(aes(size=buble_Size, color=-1*log(buble_color))) +
        scale_color_gradient(low = pal[1], high = pal[2])

    p = p + labs(color = paste0('-log10(', color,')'),
                 size = size,
                 x = x, y = y,
                 title = title) +
        theme_bw() +
        theme(axis.text.x = element_text(size=10,angle = 90, hjust = 0.5, vjust = 0.5),
              axis.text.y = element_text(color = "black",size=12),
              axis.title.x = element_text(color = "black",size=14),
              axis.title.y = element_blank())


    return(p)
}


