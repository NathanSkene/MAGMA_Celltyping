
get_gene_info <- function(genelist,
                          attributes=c("size"),
                          one_row_per_gene=F){
    genelist <- unique(genelist)
    attributes <- unique(c("hgnc_symbol","ensembl_gene_id",attributes, "start_position","end_position"))
    message("Retrieving info from biomaRt for ", formatC(length(genelist),big.mark = ",")," genes.")
    ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
    # filters = biomaRt::listFilters(ensembl)
    gene_info = biomaRt::getBM(attributes=attributes[!attributes %in% "size"],
                               filters="hgnc_symbol", 
                               values=genelist, 
                               mart=ensembl,
                               uniqueRows = T)
    gene_info$size <- gene_info$end_position - gene_info$start_position
    if(one_row_per_gene){
        gene_info <- gene_info %>%
            dplyr::group_by(hgnc_symbol) %>%
            dplyr::slice_head(n = 1) %>%
            data.frame()
        rownames(gene_info) <- gene_info$hgnc_symbol
    }
    message("Done.")
    return(gene_info)
}
