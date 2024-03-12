

#' Normalize a count matrix
#'
#'@description This function takes a count matrix, calculates either DESeq sizeFactors or EdgeR normFactors and normalizes counts using limma::voom
#'
#' @param count.matrix raw count matrix (genes are rows, columns are samples)
#' @param DESeq logical; if set to FALSE (default), edgeR normFactors will used for normalization; otherwise, DEseq sizeFactors will be used
#'
#' @return normalized count matrix
#' @export
#'
#' @examples countMatrix_norm(expr)
countMatrix_norm <- function(count.matrix, DESeq=FALSE){
    if (DESeq) {
        nf <- DESeq2::estimateSizeFactorsForMatrix(count.matrix)   # DESeq2 package of estimateSizeFactors()
        v <- limma::voom(count.matrix,
                         lib.size = median(colSums(count.matrix)) * nf)  # lib.size: for the library size (sequencing depth) for each sample,如果不自定义， the library sizes will be computed from the column sums of the counts
    } else {
        dge <- edgeR::DGEList(count.matrix)
        dge <- edgeR::calcNormFactors(dge)   # edgeR  using TMM (trimmed mean of M-values) defaultly
        v <- limma::voom(dge)
    }

    v = as.data.frame(v$E)  # E presents the data in the list generated from voom function

    return(v)
}




#' Limit max or min value in the matrix
#' @d@description
#' Cap values in a numeric matrix at a specific value for heatmap visulization
#'
#'
#' @param m normalized or scalized matrix of exprSet
#' @param max.z max z score
#'
#' @return  limit max value to max.z p; min value to -max.z
#' @export
#'
#' @examples limitRange(matirx, max.z = 3)
limitRange <- function( m, max.z = 3){
    m[ m > max.z] <- max.z
    m[ m < -max.z] <- -max.z
    return(m)
}






#' get a dataframe of gene length from a GTF file
#'
#' @param gtf.file.path the file path of target GTF file
#'
#' @return a dataframe of geneID(esmbleID) and gene length
#' @export
#'
#' @examples gene_length(gtf.file.path = '../../!公共资源/gencode.v36.annotation.gtf.gz')
gene_length <- function(gtf.file.path){

    ## gtf object
    txdb <- GenomicFeatures::makeTxDbFromGFF(gtf.file.path, format="gtf")
    ## get start and end site
    exons_gene <- GenomicFeatures::exonsBy(txdb, by = "gene")
    ## avoid overlap
    exons_gene_lens <- lapply(exons_gene,
                              function(x){sum(IRanges::width(IRanges::reduce(x)))})

    ## transform to dataframe
    gene.length = do.call(rbind, lapply(exons_gene_lens, data.frame))
    gene.length = data.frame(gene_id = rownames(gene.length), effLenth = gene.length[,1])

    ## gene_id as esemble ID formation
    gene.esemble = stringr::str_split(gene.length$gene_id, pattern = '\\.', simplify = T)[  ,1] # remove the version
    gene.length$esmID =  gene.esemble
    gene.length = gene.length[!duplicated(gene.length$esmID),  ]  # remove duplicated
    rownames(gene.length) = gene.length$esmID

    return(gene.length)
}




#' Transform count matrix to tmp matrix
#'
#' @param count.matrix raw count matrix (genes are rows, columns are samples); rownames are esemble ID
#' @param gene.length.ob object get from gene_length function; a dataframe of geneID(esmbleID) and gene length; rownames are esemble ID
#'
#' @return a tpm matrix (genes are rows, columns are samples); rownames are esemble ID
#' @export
#'
#' @examples counts2TPM(count.matrix, gene.length.ob)
counts2TPM <- function(count.matrix, gene.length.ob){

    ## common genes
    gene = intersect(rownames(count.matrix), rownames(gene.length.ob))

    ## consistant dataframe using common genes
    tmp.count =  data.frame(count.matrix[gene, ], check.names = F)
    tmp.length = gene.length.ob[gene, ]

    ## if the gene orders is consistant
    if(identical(rownames(tmp.count),  rownames(tmp.length))){
        cat('count matrix and gene lenth dataframe is in the same order! Must in esemble ID!')
    }else{
        stop( "not consistant data!" )
    }

    ## function for single gene:count transformed to TPM matrix
    transform_math = function(counts, GeneLength){
        rate = log(counts) - log(GeneLength)
        denom = log(sum(exp(rate)))
        exp(rate - denom + log(1e6))
    }

    ## count matrix to TPM
    expr.TPM = apply(tmp.count, 2, transform_math, GeneLength=tmp.length$effLenth)

    return(expr.TPM)
}




#' get a dataframe of esembleID to gene symbol with gene type  from a GTF file
#'
#' @param gtf.file.path the file path of target GTF file
#'
#' @return  a dataframe of esembleID, gene symbol, gene type
#' @export
#'
#' @examples esmble2symbol(gtf.file.path = '../../!公共资源/gencode.v36.annotation.gtf.gz')
esmble2symbol <- function(gtf.file.path){

    ## get gtf file
    GTF.info = data.table::fread(gtf.file.path, header = F)

    ## get gene symbol id
    input = GTF.info[GTF.info[[3]]=='gene',  ]  # [[3]] gene a vector ，[,3] a dataframe

    ## get gene information
    Esemble = "gene_id \"([^\\s]+)\";"  # \\s  blank string; [^\\s] non blank string
    Genesymbol = "gene_name \"([^\\s]+)\";"
    geneType = "gene_type \"([^\\s]+)\";"

    ## Get Esemble_id
    Esemble_id = stringr::str_extract(input[[9]], Esemble)
    Esemble_id = stringr::str_extract(Esemble_id, "\"[^\\s]+\"")  # get strings in the ""
    Esemble_id = stringr::str_replace_all(Esemble_id, pattern = "\"", replacement = "") # remove ""
    Esemble_id = stringr::str_split(Esemble_id, pattern = '\\.', simplify = T)[ ,1]  # remove version

    ## Get gene Symbol
    Genesymbol_id = stringr::str_extract(input[[9]], Genesymbol)
    Genesymbol_id = stringr::str_extract(Genesymbol_id , "\"[^\\s]+\"")
    Genesymbol_id = stringr::str_replace_all(Genesymbol_id , pattern = "\"", replacement = "")

    ## Get Gene tpye
    genetype = stringr::str_extract(input[[9]], geneType)
    genetype = stringr::str_extract(genetype, "\"[^\\s]+\"")
    genetype = stringr::str_replace_all(genetype , pattern = "\"", replacement = "")

    ## form data frame
    ids = data.frame(Esemble = Esemble_id,
                     genSymbol = Genesymbol_id,
                     Type = genetype)
    return(ids)
}



