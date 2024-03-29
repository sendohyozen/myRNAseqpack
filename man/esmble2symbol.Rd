% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expr_normalization.R
\name{esmble2symbol}
\alias{esmble2symbol}
\title{get a dataframe of esembleID to gene symbol with gene type  from a GTF file}
\usage{
esmble2symbol(gtf.file.path)
}
\arguments{
\item{gtf.file.path}{the file path of target GTF file}
}
\value{
a dataframe of esembleID, gene symbol, gene type
}
\description{
get a dataframe of esembleID to gene symbol with gene type  from a GTF file
}
\examples{
esmble2symbol(gtf.file.path = '../../!公共资源/gencode.v36.annotation.gtf.gz')
}
