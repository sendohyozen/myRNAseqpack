% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GO_KEGG.R
\name{my_bublePlot}
\alias{my_bublePlot}
\title{bubleplot for go,kegg analysis}
\usage{
my_bublePlot(df, x, y, size, color, title = "KEGG", pal = c("blue", "red"))
}
\arguments{
\item{df}{dataframe of go, kegg results}

\item{x}{column name of x axis}

\item{y}{column name of y axis}

\item{size}{column name showing size of the buble}

\item{color}{column name showing  color of the buble}

\item{title}{title of the plot , default "KEGG"}

\item{pal}{color palette of the bulbe color; a vector consists of two numbers, the fisrt for low , the second for high}
}
\value{
ggplot2 object of a buble plt
}
\description{
bubleplot for go,kegg analysis
}
\examples{
my_bublePlot(df = data, x = 'pvalue', y = 'Description', size = 'Count' , color = 'p.adjust', title = 'KEGG')
}
