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




#' matchstick_plot for go,kegg,gsea... analysis
#'
#' @param df dataframe of go, kegg,... results
#' @param top top term picked for plotting , arranged by count and pvalue , default top 20
#' @param y.decreasing if arranged the sticked with derceasing order, default T
#' @param logPV10 if use the -log10pvalue for filling colors, default T
#' @param x column name of x axis data, (the plot is coordfliped) default Description
#' @param y column name of y axis data, (the plot is coordfliped) default count
#' @param size column name of point size, default count
#' @param fill_col column name of filling color, default pvalue
#' @param fill_col_pal color palette of the bulbe color; a vector consists of two numbers, the fisrt for low , the second for high;
#' @param stick_col color of the stick
#' @param stick_width width of the stick, default 2
#' @param title title of the plot, default NULL
#' @param xlab x lab of the plot, default NULL
#' @param ylab y lab of the plot, default 'Gene Counts'
#' @param x.lab.size text size of the x labs, default 14
#' @param x.lab.angle text angle of the x labs, default 0
#' @param y.lab.size text size of the y labs, default 14
#' @param y.lab.angle text angle of the y labs, default 0
#' @param x.title.size text size of the x lab title, default 18
#' @param y.title.size text size of the y lab title, default 18
#' @param title.size text size of the plot title, default 20
#' @param col.limits limit number of the color bar, default NULL, if given, a vector consists of two numbers, the fisrt for min , the second for max
#' @param col.breaks break numbers of the color bar, if given, a vector consists of numbers, notes the numbers in the limits
#' @param col.labels labs of the break numbers , be consistent with the col.breaks
#' @param size.limits limit number of the size bar, default NULL, if given, a vector consists of two numbers, the fisrt for min , the second for max
#' @param size.breaks break numbers of the size bar, if given, a vector consists of numbers, notes the numbers in the limits
#' @param size.labels labs of the break numbers , be consistent with the size.breaks
#'
#' @return a ggplot2 object of matchstick plot
#' @export
#'
#' @examples my_matchstick_plot(df = GOdatframe, top = 20, y.decreasing = T,
#' x = 'Description', y = 'Count', size = 'Count', fill_col = 'pvalue',
#' col.breaks = c(3,4,5) , col.labels = c('low', '', 'high'),size.limits = c(20,40), size.breaks = c(20,30,40))

my_matchstick_plot <- function(df,  top=20, y.decreasing = T, logPV10 = T,
                               x='Description', y="count", size="count", fill_col="pvalue",
                               fill_col_pal = c("#BDD8E8", "#134B6C"), stick_col = "#FAAC90", stick_width=2,
                               title =NULL, xlab = NULL, ylab = "Gene counts",
                               x.lab.size=14, x.lab.angle=0, y.lab.size =14, y.lab.angle=0, x.title.size=18, y.title.size=18, title.size=20,
                               col.limits=NULL, col.breaks, col.labels,
                               size.limits=NULL, size.breaks, size.labels){

    ## plot data
    plot_df = data.frame(x_value = df[[x]],
                         y_value = df[[y]],
                         size_value =  df[[size]],
                         fill_col_value = df[[fill_col]])


    if(logPV10 == T){
        plot_df = plot_df %>% as.data.frame() %>% mutate(fill_col_value = -log10(fill_col_value))
        color.title = '-log10(pvalue)'

    }else{
        color.title = fill_col
        warning('Do not need to calculate the -log10(Pvalue) ? \n')
    }

    # arranged by decreasing y value (count) and pvalue
    plot_df = plot_df %>% arrange(desc(y_value), desc(fill_col_value)) %>% head(top) %>% as.data.frame()

    if(y.decreasing == T){
        # ordered x names
        plot_df$x_value = factor(plot_df$x_value, levels = rev(plot_df$x_value))
    }else{
        plot_df$x_value = factor(plot_df$x_value, levels = plot_df$x_value)
    }




    # plotting
    p <- ggplot(plot_df, aes(x=x_value, y=y_value)) +
        geom_segment(aes(x = x_value, xend = x_value,
                         y = 0, yend = y_value),
                     linewidth = stick_width,
                     color = stick_col,
                     linetype="solid") +
        geom_point(aes(color = fill_col_value,
                       size = size_value)) +
        scale_color_continuous(low = fill_col_pal[1], high = fill_col_pal[2]) +
        labs(x = xlab, y = ylab, title = title, color=color.title, size = size) +
        coord_flip()

    # set break points for continuous color bar
    if(!missing(col.breaks)){
        # defualt break = labs
        if(missing(col.labels)){
            col.labels = col.breaks
        }

        p = p + scale_color_gradient(low = fill_col_pal[1], high = fill_col_pal[2],
                                     limits = col.limits,
                                     breaks = col.breaks,
                                     labels = col.labels)
    }

    # set break points for continuous size bar
    if(!missing(size.breaks)){
        # defualt break = labs
        if(missing(size.labels)){
            size.labels = size.breaks
        }
        p = p + scale_size(limits = size.limits,
                           breaks = size.breaks,
                           labels = size.labels)
    }

    # theme
    p = p + theme_classic() +
        theme(plot.title = element_text(size = title.size, hjust=0.5),
              axis.text.x = element_text(size= x.lab.size, angle = x.lab.angle, vjust = 1),
              axis.text.y = element_text(size= y.lab.size, angle = y.lab.angle),
              axis.title.x = element_text(size = x.title.size),
              axis.title.y = element_text(size = y.title.size))

    return(p)
}

