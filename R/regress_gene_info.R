#' Regress out the effect of gene attributes
#' 
#' Queries \code{biomaRt} for \code{gene_attributes} and then regresses them out of \code{xmat}.
#' 
#' @param xmat gene x sample matrix.
#' @param attributes Gene attributes to extract from \code{biomaRt}
#'  and then regress from \code{xmat}.
#' @inheritParams iterate_lm
#'  
#' @examples
#' library(MAGMA.Celltyping)
#' xmat <- MAGMA.Celltyping::ctd_DRONC_human[[1]]$mean_exp[1:100,]
#' adjusted_df <- regress_gene_info(xmat) 
regress_gene_info <- function(xmat,
                              attributes=c("size"),
                              correction_method="BH"){
    gene_info <- get_gene_info(genelist = rownames(xmat),
                               attributes=attributes,
                               one_row_per_gene = T) 
    gene_intersect <- intersect(rownames(xmat), gene_info$hgnc_symbol)
    message(length(gene_intersect)," intersecting genes between xmat and bioMart query.") 
    xdat <- as.matrix(xmat[gene_intersect,]) 
    ydat <- as.matrix(gene_info[gene_intersect,attributes]) 
    #### Run model
    message("+ Training model")
    mod <- stats::lm(xdat ~ ydat) 
    message("+ Generating predictions")
    adjusted_df <- stats::predict(mod, data.frame(xdat))
    return(adjusted_df) 
}
