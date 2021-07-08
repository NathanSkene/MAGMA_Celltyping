#' Filter genes in a CellTypeDataset
#' 
#' 
#' 
#' @param ctd CellTypeDataset
#' @param gene_subset Genes to subset to.
#' 
#' @examples
#' library(MAGMA.Celltyping)
#' library(ewceData)
#' set.seed(1)
#' ctd <- ewceData::ctd()
#' gene_subset <- sample(rownames(ctd[[1]]$mean_exp), 1000)
#' ctd_subset <- filter_ctd_genes(ctd = ctd, gene_subset = gene_subset)
#' 
#' @export
filter_ctd_genes <- function(ctd, 
                             gene_subset){
    message("+ Filtering CTD to ",length(gene_subset)," genes.")
    new_ctd <- lapply(1:length(ctd), function(lvl){
        message("level ",lvl)
        ctd_lvl <- ctd[[lvl]]
        mat_names <- names(ctd_lvl)[!names(ctd_lvl) %in% c("annot","plotting")]
        other_names <- names(ctd_lvl)[names(ctd_lvl) %in% c("annot","plotting")]
        new_ctd_lvl <- lapply(mat_names, function(mat_nm){
            message("   - ",mat_nm)
            ctd_lvl[[mat_nm]][rownames(ctd_lvl[[mat_nm]]) %in% gene_subset, ]
        }) %>% `names<-`(mat_names) 
        for(nm in other_names){
            new_ctd_lvl[nm] <- ctd_lvl[nm]
        }
        return(new_ctd_lvl)
    }) %>% `names<-`(paste0("level",1:length(ctd)))
    return(new_ctd)
}
