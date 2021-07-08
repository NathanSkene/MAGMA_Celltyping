#' Get MAGMA Paths
#'
#' Gets paths to folders based on location of the GWAS summary statistics file. This avoids the GWAS sumstats folder becoming a mess.
#'
#' @param gwas_sumstats_path Path to the GWAS summary statistics file
#' @param upstream_kb How many kb upstream of the gene should SNPs be included?
#' @param downstream_kb How many kb downstream of the gene should SNPs be included?
#' @param output_path If being called from magma.tileplot() then gwas_sumstats_path cannot be specified, so instead specificify the base folder where MAGMA_Figures exists
#'
#' @return Filepath for where MAGMA files should be created
#'
#' @examples
#' magmaPaths = get.magma.paths(gwas_sumstats_path,upstream_kb,downstream_kb)
#'
#' @export
get.magma.paths <- function(gwas_sumstats_path=NA,
                            upstream_kb=NA,
                            downstream_kb=NA,
                            output_path=NA){
    # Error checking
    if(is.na(gwas_sumstats_path) & is.na(output_path)){stop("Either gwas_sumstats_path or output_path must be specified")}
    #if(!is.na(gwas_sumstats_path) & !is.na(output_path)){stop("Only gwas_sumstats_path or output_path should be specified, not both")}
    if(!is.na(gwas_sumstats_path) & (is.na(upstream_kb) | is.na(downstream_kb))){
        stop("If gwas_sumstats_path is specified then upstream_kb and downstream_kb must also be specified")
    }
    
    # If output_path is not specified, then use the folder containing gwas_sumstats
    if(is.na(output_path)){ output_path = dirname(gwas_sumstats_path) }
    
    if(!is.na(gwas_sumstats_path)){
        gwasFileName = basename(gwas_sumstats_path)
        
        # Set the paths
        pathMagmaFiles = sprintf("%s/MAGMA_Files/%s.%sUP.%sDOWN",output_path,gwasFileName,upstream_kb,downstream_kb)
        pathFigs = sprintf("%s/MAGMA_Figures/%s FIGS",output_path,gwasFileName)
        pathTiles = sprintf("%s/MAGMA_Figures/Tileplots",output_path)
        
        # Create the folders (in case they don't exist yet)
        dir.create(sprintf("%s/MAGMA_Files",output_path), showWarnings = FALSE)
        dir.create(pathMagmaFiles, showWarnings = FALSE)
        dir.create(sprintf("%s/MAGMA_Figures",output_path), showWarnings = FALSE)
        dir.create(pathTiles, showWarnings = FALSE)    
        dir.create(pathFigs, showWarnings = FALSE)    
        
        prefix = sprintf("%s.%sUP.%sDOWN",gwasFileName,upstream_kb,downstream_kb)
        filePathPrefix = sprintf("%s/%s.%sUP.%sDOWN",pathMagmaFiles,gwasFileName,upstream_kb,downstream_kb)
        
        # Return the paths
        magmaPaths = list(tiles=pathTiles,figs=pathFigs,files=pathMagmaFiles,prefix=prefix,filePathPrefix=filePathPrefix,gwasFileName=gwasFileName, gwasFilePath=gwas_sumstats_path)    
    }else{
        dir.create(sprintf("%s/MAGMA_Figures",output_path), showWarnings = FALSE)   
        pathTiles = sprintf("%s/MAGMA_Figures/Tileplots",output_path)
        dir.create(pathTiles, showWarnings = FALSE)   
        magmaPaths = list(tiles=pathTiles)
    }

    return(magmaPaths)
}
