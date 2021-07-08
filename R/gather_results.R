

#' Gather enrichment results from \code{celltype_associations_pipeline}
#' 
#' @examples 
#' \dontrun{
#' library(MAGMA.Celltyping)
#' local_files <- import_magma_files(download_dir=".")
#' #' magma_dirs <- unique(dirname(local_files))
#' res <- celltype_associations_pipeline(ctd=ewceData::ctd(), 
#'                                       ctd_name="Zeisel2018", 
#'                                       magma_dirs=magma_dirs, 
#'                                       genome_ref_path="~/Downloads/g1000_eur/g1000_eur") 
#' merged_res <- gather_results(res)
#' }
#' @export
gather_results <- function(MAGMA_results,
                          level=2, 
                          dataset_name=NULL,
                          species="mouse",
                          filetype="ctAssocMerged", 
                          fdr_thresh=NULL,
                          save_dir=NULL){ 
    merged_results <- lapply(names(MAGMA_results), function(nm){
        cbind(GWAS=stringr::str_split(nm,".annotated")[[1]][1],
              MAGMA_results[[nm]][[filetype]][[level]]$results
              # MAGMA_results[[nm]]$ctAssocsLinear[[level]]$results
        )
    }) %>% data.table::rbindlist() %>% 
        # Remove unlabeled cell clusters
        subset(!startsWith(Celltype,"X")) %>%
        dplyr::mutate(FDR=p.adjust(p = P, method="fdr"),
                      Celltype_id=Celltype,
                      Celltype=gsub(paste(dataset_name,species,"[.]",sep="|")," ", Celltype)) %>% 
        dplyr::arrange(FDR)
    ### Save 
    if(!is.null(save_dir)){
        save_path <- paste0(save_dir,"/MAGMA_celltyping.",dataset_name,".lvl",level,".csv")
        message("+ Saving full merged results to ==>",save_path)
        data.table::fwrite(merged_results, save_path) 
    } 
    ## Filter to only FDR-sig results
    if(!is.null(fdr_thresh)){
        message("+ Filtering results @ FDR<",fdr_thresh)
        merged_results <- subset(merged_results, FDR<fdr_thresh)
        message("+ ",nrow(merged_results)," significant results remaining.")
    }
    return(merged_results)
}
