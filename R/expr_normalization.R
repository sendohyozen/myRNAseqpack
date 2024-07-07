

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




#' Transform count matrix to TPM matrix
#'
#' @param esemble.length.df data frame object get from gene_length function; a dataframe of geneID(esmbleID) and gene length; rownames are esemble ID
#' @param RownameIsSymbol if the count matrix given have rownames as gene symble, default F
#' @param gtf.file.path Only needed if gene symbol is given as rownames; the file path of target GTF file; default NULL
#' @param count.matrix raw count matrix (genes are rows, columns are samples); rownames are esemble ID/gene symbol
#'
#' @return a tpm matrix (genes are rows, columns are samples); rownames are esemble ID
#' @export
#'
#' @examples  counts2TPM(count.matrix = exprSet, esemble.length.df = gene_length, RownameIsSymbol = T, gtf.file.path = '../../!公共资源/gencode.v36.annotation.gtf.gz')
counts2TPM <- function(count.matrix, esemble.length.df, RownameIsSymbol=F, gtf.file.path){

    # count 矩阵必须全为数字，行为基因名或id，列为样本
    cat('the count matrix must be numeric, gene symbol/esemble ID in the rows, and sample names in the column! \n')
    is_numeric_df <- function(df) {
        all_numeric <- sapply(df, is.numeric)
        return(all(all_numeric))
    }

    stopifnot(is_numeric_df(count.matrix))


    # count矩阵行名是否为gene symbol需要对长度数据框进行调整
    if(RownameIsSymbol==T){

        # 强制对表达矩阵的基因行名取大写（genesymbol正式gtf文件中都是大写）
        rownames(count.matrix) = toupper(rownames(count.matrix))


        cat('the count matrix give rownames as gene symbol! must give gtf file path for id transformation! \n')

        ids = esmble2symbol(gtf.file.path = gtf.file.path)

        gene.length.ob = esemble.length.df %>% left_join(ids, by = c('esmID'='Esemble'))  %>%
            select(genSymbol, effLenth) %>%
            distinct(genSymbol, .keep_all = T) %>%
            column_to_rownames('genSymbol') %>%
            mutate(genesymbol = rownames(.)) %>%
            as.data.frame()

    }else{
        cat('the count matrix give rownames as esemble ID!  \n')
        gene.length.ob = as.data.frame(esemble.length.df)
    }


    ## common genes
    gene = intersect(rownames(count.matrix), rownames(gene.length.ob))
    cat(paste0(length(gene)), ' genes could be transformed from count to tpm!!\n')

    ## consistant dataframe using common genes
    tmp.count =  data.frame(count.matrix[gene, ], check.names = F)
    tmp.length = data.frame(gene.length.ob[gene, ], check.names = F)

    ## if the gene orders is consistant (gene.length.ob 对象不能是单列，否则会被转成向量，丢失行名)
    if(identical(rownames(tmp.count),  rownames(tmp.length))){
        cat('count matrix and gene lenth dataframe is in the same order!')
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



