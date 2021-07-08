


#' Calculate celltype associations using MAGMA
#'
#' Has the option of running multiple analyses with a single function.  
#' Assumes that you already have all MAGMA GWAS files precomputed.
#' Precomputed MAGMA GWAS files can be downlaoded via the \code{import_magma_files} fucnction. 
#'
#' @param ctd Cell type data structure containing $specificity_quantiles
#' @param ctd_name Used in file names 
#' @param magma_dirs Names of folders containing the pre-computed MAGMA GWAS files. 
#' \bold{NOTE}: Files within these folders must have the same naming scheme as the folders themselves.  
#' @param run_linear Run \code{calculate_celltype_associations} in linear mode.
#' @param run_top10 Run \code{calculate_celltype_associations} in top 10% mode.
#' @param run_condition Run \code{calculate_celltype_associations} in conditional mode.
#' @param suffix_linear This will be added to the linear results file name.
#' @param suffix_top10 This will be added to the top 10% results file name.
#' @param suffix_condition This will be added to the conditional results file name.
#' @inheritParams calculate_celltype_associations
#' @inheritParams calculate_conditional_celltype_associations
#' @param save_dir Folder to save results in (\code{save_dir=NULL} to not save any results).
#' 
#' @return A list containing the results of each selected celltype associations analysis.
#'
#' @examples
#' \dontrun{
#' library(MAGMA.Celltyping)
#' local_files <- import_magma_files(download_dir=".")
#' magma_dirs <- unique(dirname(local_files))
#' res <- celltype_associations_pipeline(ctd=ewceData::ctd(), ctd_name="Zeisel2018", magma_dirs=magma_dirs, genome_ref_path="~/Downloads/g1000_eur/g1000_eur")
#' }
#' @export
celltype_associations_pipeline <- function(ctd, 
                                          ctd_name,
                                          magma_dirs,
                                          genome_ref_path,
                                          specificity_species="mouse",
                                          run_linear=T,
                                          run_top10=T,
                                          run_condition=F,
                                          upstream_kb=35,
                                          downstream_kb=10, 
                                          suffix_linear="linear",
                                          suffix_top10="top10",
                                          suffix_condition="condition",
                                          controlTopNcells=1,
                                          force_new=F,
                                          save_dir=tempdir()){ 
    ### Example 
    # magma_dir<- "/Users/schilder/Desktop/model_celltype_conservation/raw_data/MAGMA/MAGMA_Files/ADHD_iPSYCH.all.annotated.35UP.10DOWN";upstream_kb=35;downstream_kb=10; specificity_species="human"; genome_ref_path=path.expand("~/Desktop/model_celltype_conservation/raw_data/MAGMA/g1000_eur/g1000_eur"); suffix_linear="linear";suffix_top10="top10";suffix_condition="condition";controlTopNcells=1; ctd_name="Aerts2021"
    
    
    ## Establish vars in case some are not computed.
    ctAssocsLinear <- NULL; ctAssocsTop <- NULL; ctCondAssocs <- NULL; ctAssocMerged <- NULL;
    
    MAGMA_results <- lapply(magma_dirs, function(magma_dir){
        message(basename(magma_dir))  
        fake_gwas_ss <- file.path(gsub("/MAGMA_Files","",dirname(magma_dir)),
                                  gsub(paste0(".",upstream_kb,"UP.",downstream_kb,"DOWN"),"",basename(magma_dir)))   
        #### Linear mode ####
        if(run_linear){
            message("+ Calculating celltype associations: linear mode") 
            ctAssocsLinear <- tryCatch(expr = {
                calculate_celltype_associations(ctd = ctd,
                                                gwas_sumstats_path = fake_gwas_ss,
                                                genome_ref_path=genome_ref_path,
                                                upstream_kb = upstream_kb,
                                                downstream_kb = downstream_kb,
                                                analysis_name = paste(ctd_name,suffix_linear,sep="_"),
                                                specificity_species = specificity_species,
                                                force_new = force_new)
            }, error = function(e){message(e); return(NULL)})
        }
        
        #### Top 10% mode ####
        if(run_top10){
            message("+ Calculating celltype associations: top10% mode") 
            ctAssocsTop <- tryCatch(expr = { 
                calculate_celltype_associations(ctd = ctd,
                                                gwas_sumstats_path=fake_gwas_ss,
                                                genome_ref_path=genome_ref_path,
                                                EnrichmentMode="Top 10%",
                                                upstream_kb = upstream_kb,
                                                downstream_kb = downstream_kb,
                                                analysis_name = paste(ctd_name,suffix_top10,sep="_"),
                                                specificity_species = specificity_species,
                                                force_new = force_new)
            },  error = function(e){message(e); return(NULL)})
        }
        
        
        #### Merge results ####
        if(all(!is.null(ctAssocsLinear), !is.null(ctAssocsTop))){
            message("+ Merging linear and top10% results")
            ctAssocMerged = merge_magma_results(ctAssoc1=ctAssocsLinear,
                                                ctAssoc2=ctAssocsTop) 
        }
        
        #### Conditional mode ####
        if(run_condition){ 
            message("+ Calculating celltype associations: conditional mode") 
            ctCondAssocs <- tryCatch({
                calculate_conditional_celltype_associations(ctd = ctd,
                                                            gwas_sumstats_path = fake_gwas_ss,
                                                            genome_ref_path=genome_ref_path,
                                                            analysis_name = paste(ctd_name,suffix_condition,sep="."),
                                                            upstream_kb = upstream_kb,
                                                            downstream_kb = downstream_kb,
                                                            controlTopNcells=controlTopNcells,
                                                            specificity_species = specificity_species)
            }, error=function(e){message(e); return(NULL)})
        }
        
        return(list(magma_dir=magma_dir,
                    ctAssocsLinear=ctAssocsLinear,
                    ctAssocsTop=ctAssocsTop,
                    ctAssocMerged=ctAssocMerged,
                    ctCondAssocs=ctCondAssocs))
    }) %>% `names<-`(basename(magma_dirs)) 
    
    
    if(!is.null(save_dir)){
        save_path <- file.path(save_dir,ctd_name,paste0("MAGMA_celltyping.",ctd_name,".rds"))
        print(paste("+ Saving results ==>",save_path))
        dir.create(dirname(save_path), showWarnings = F, recursive = T)
        saveRDS(MAGMA_results, save_path) 
    }
    return(MAGMA_results)
}
