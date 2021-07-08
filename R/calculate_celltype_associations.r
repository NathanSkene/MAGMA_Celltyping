#' Calculate celltype associations using MAGMA
#'
#' Assumes that you have already run \code{map.snps.to.genes()}
#'
#' @param ctd Cell type data structure containing $specificity_quantiles
#' @param gwas_sumstats_path File path of the summary statistics file
#' @param analysis_name Used in file names which area created
#' @param upstream_kb How many kb upstream of the gene should SNPs be included?
#' @param downstream_kb How many kb downstream of the gene should SNPs be included?
#' @param genome_ref_path Path to the folder containing the 1000 genomes .bed files (which can be downloaded 
#' from https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip)
#' @param specificity_species Species name relevant to the cell type data, i.e. "mouse" or "human"
#' @param genesOutCOND [Optional] Path to a genes.out file to condition on. Used if you want to condition on a different GWAS.
#' @param EnrichmentMode [Optional] Should either 'Linear' or 'Top 10\%' mode be used for testing enrichment?
#' @param force_new [Optional] Force new MAGMA analyses even if the pre-existing results files are detected.
#' @return File path for the genes.out file
#'
#' @examples
#' \dontrun{
#' ctAssocs = calculate_celltype_associations(ewceData::ctd(), gwas_sumstats_path, gwas_sumstats_path = "<PATH_TO_SUMSTATS>", upstream_kb=10, downstream_kb=1.5, genome_ref_path="<PATH_TO_g1000_eur>", specificity_species="mouse", EnrichmentMode="Top 10%")
#' }
#'
#' @export
calculate_celltype_associations <- function(ctd,
                                             gwas_sumstats_path, 
                                             analysis_name="MainRun",
                                             upstream_kb=10,
                                             downstream_kb=1.5,
                                             genome_ref_path,
                                             specificity_species="mouse",
                                             genesOutCOND=NA,
                                             EnrichmentMode="Linear",
                                             force_new=F){ 
    # Check EnrichmentMode has correct values
    if(!EnrichmentMode %in% c("Linear","Top 10%")){stop("EnrichmentMode argument must be set to either 'Linear' or 'Top 10%")}
    
    gwas_sumstats_path = path.expand(gwas_sumstats_path)
    magmaPaths = get.magma.paths(gwas_sumstats_path,upstream_kb,downstream_kb)
    
    # Check for errors in arguments
    check_inputs_to_magma_celltype_analysis(ctd,gwas_sumstats_path,analysis_name,upstream_kb,downstream_kb,genome_ref_path)
    
    output = list()
    for(annotLevel in 1:length(ctd)){
        # Prepare output list
        tmp = list() 
        sumstatsPrefix2 = sprintf("%s.level%s",magmaPaths$filePathPrefix, annotLevel) 
        path = sprintf("%s.%s.gsa.out", sumstatsPrefix2, analysis_name)
        #### Need to make sure colnames still match after theyre edited by prepare.quantile.groups step
        ctd <- standardise_ctd(ctd = ctd, lvl = annotLevel) 
        geneCovarFile <- NULL 
        print(path)
        
        if((!file.exists(path)) | force_new){ 
            message("+ Running MAGMA") 
            
            if(EnrichmentMode=="Linear"){
                # First match quantiles to the genes in the genes.out file... then write as the genesCovar file (the input to MAGMA)
                geneCovarFile = create_gene_covar_file(genesOutFile = sprintf("%s.genes.out",magmaPaths$filePathPrefix),ctd,annotLevel,specificity_species=specificity_species,genesOutCOND)
                
                if(is.na(genesOutCOND[1])){
                    magma_cmd = sprintf("magma --gene-results '%s.genes.raw' --gene-covar '%s' --model direction=pos --out '%s.%s'",magmaPaths$filePathPrefix,geneCovarFile,sumstatsPrefix2,analysis_name)
                }else{
                    conditionOn = paste(sprintf("ZSTAT%s",1:length(genesOutCOND)),collapse=",")
                    magma_cmd = sprintf("magma --gene-results '%s.genes.raw' --gene-covar '%s' --model direction=pos  condition-residualize='%s' --out '%s.%s'",magmaPaths$filePathPrefix,geneCovarFile,conditionOn,sumstatsPrefix2,analysis_name)
                }
            }else if(EnrichmentMode=="Top 10%"){
                # First match quantiles to the genes in the genes.out file... then write as the genesCovar file (the input to MAGMA)
                geneCovarFile = create_top10percent_genesets_file(genesOutFile = sprintf("%s.genes.out",magmaPaths$filePathPrefix),
                                                                  ctd = ctd,
                                                                  annotLevel = annotLevel,
                                                                  specificity_species=specificity_species)
                
                if(is.na(genesOutCOND[1])){
                    magma_cmd = sprintf("magma --gene-results '%s.genes.raw' --set-annot '%s' --out '%s.%s'",magmaPaths$filePathPrefix,geneCovarFile,sumstatsPrefix2,analysis_name)
                }else{
                    geneCovarFile2 = create_gene_covar_file(genesOutFile = sprintf("%s.genes.out",magmaPaths$filePathPrefix),ctd,annotLevel,specificity_species=specificity_species,genesOutCOND)
                    conditionOn = paste(sprintf("ZSTAT%s",1:length(genesOutCOND)),collapse=",")
                    magma_cmd = sprintf("magma --gene-results '%s.genes.raw' --set-annot '%s' twosided --gene-covar '%s' condition-only='%s' --out '%s.%s'",magmaPaths$filePathPrefix,geneCovarFile,geneCovarFile2,conditionOn,sumstatsPrefix2,analysis_name)
                }
            }
            print(magma_cmd)
            system(magma_cmd)
        } else {
            message("+ Importing precomputed MAGMA results...")
        }
        tmp$geneCovarFile = geneCovarFile
        tmp$results = load.magma.results.file(path = path,
                                              annotLevel = annotLevel,
                                              ctd = ctd,
                                              genesOutCOND=genesOutCOND,
                                              EnrichmentMode=EnrichmentMode,
                                              ControlForCT="BASELINE")
        output[[length(output)+1]] = tmp
    } # //End for loop
    
    # Calculate total number of tests performed
    totalTests = 0
    for(annotLevel in 1:length(output)){
        totalTests = totalTests + dim(output[[annotLevel]]$results)[1]
    }
    output$total_baseline_tests_performed = totalTests
    
    output$gwas_sumstats_path = gwas_sumstats_path
    output$analysis_name = analysis_name
    output$upstream_kb = upstream_kb
    output$downstream_kb = downstream_kb
    output$genome_ref_path = genome_ref_path
    
    return(output)
}


standardise_ctd <- function(ctd, lvl){
    #### Need to make sure colnames still match after theyre edited by prepare.quantile.groups step
    message("+ Standardising CTD.")
    for(nm in names(ctd[[lvl]])){
        if((!nm %in% c("annot","plotting")) & (!is.null(dim(ctd[[lvl]][[nm]]))) ){ 
            ctd[[lvl]][[nm]] <- as.matrix(data.frame(as.matrix(ctd[[lvl]][[nm]])))
            colnames(ctd[[lvl]][[nm]]) <- make.unique(colnames(ctd[[lvl]][[nm]]))
        } 
    }
    return(ctd)
}
