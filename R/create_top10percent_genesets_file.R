#' Create gene covar file
#'
#' The gene covar file is the input to MAGMA for the celltype association analysis. This code was functonalised 
#' because it is called by both baseline and conditional analysis.
#'
#' @param genesOutFile The output of the second call to MAGMA (performed in the map.snps.to.genes function)
#' @param ctd Cell type data structure. Must contain quantiles.
#' @param annotLevel Annot level for which the gene covar file should be constructed
#' @param specificity_species Species name relevant to the cell type data, i.e. "mouse" or "human"
#'
#' @return Filepath for the gene covar file
#'
#' @examples
#' library(MAGMA.Celltyping)
#' myGenesOut = tempfile()
#' data.table::fwrite(MAGMA.Celltyping::genesOut,sep="\t",file=myGenesOut)
#' ctd = prepare.quantile.groups(ctd=ewceData::ctd(),specificity_species="mouse",numberOfBins=40)
#' geneSetsFilePath = create_top10percent_genesets_file(genesOutFile=myGenesOut,ctd=ctd,annotLevel=1, specificity_species="mouse")
#' @export
create_top10percent_genesets_file <- function(genesOutFile,
                                              ctd,
                                              annotLevel,
                                              specificity_species){
    #### Map genes first so that the deciles computed in the following step only include usable genes ####
    quantDat2 <- get_top10percent(ctd = ctd, 
                                  annotLevel = annotLevel, 
                                  specificity_species = specificity_species)
    cts = setdiff(colnames(quantDat2),"entrez")
    ctRows = rep("",length(cts))
    names(ctRows) = cts
    for(ct in cts){
        ctRows[ct] = paste(c(ct,quantDat2[quantDat2[,ct]==10,1]),collapse=" ")
    } 
    # Write genes covar file to disk
    geneCovarFile=tempfile()
    #write.table(quantDat2,file=geneCovarFile,quote=FALSE,row.names=FALSE,sep="\t")
    write.table(ctRows,file=geneCovarFile,quote=FALSE,row.names=FALSE,sep="\t",col.names=FALSE)
    return(geneCovarFile)
}




get_top10percent <- function(ctd,
                             annotLevel,
                             specificity_species){ 
    ctd2 = map_specificity_to_entrez(ctd=ctd,
                                     annotLevel=annotLevel,
                                     specificity_species=specificity_species, 
                                     return_ctd = T)
    ctd3 = prepare.quantile.groups(ctd2,
                                   specificity_species=specificity_species,
                                   numberOfBins=10)
    quantDat2 <- ctd3[[annotLevel]]$quantDat2
    print(paste("+",dim(quantDat2)[1],"genes extracted from top10% specificity."))
    if(dim(quantDat2)[1]<100){stop("Less than one hundred genes detected after mapping genes between species. Was specificity_species defined correctly?")}
    return(quantDat2)
}
