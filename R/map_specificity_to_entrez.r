#' Map specificity to entrez
#'
#' Convenience function used in 'create_gene_covar_file()'
#'
#' @param ctd Cell type data structure. Must contain quantiles.
#' @param annotLevel Annot level for which the gene covar file should be constructed
#' @param specificity_species Species name relevant to the cell type data, i.e. "mouse" or "human"
#'
#' @return Matrix in which the first column is 'entrez' and then the specificity decile for each cell type
#'
#' @examples
#' genesCovarFilePath = create_gene_covar_file(genesOutFile,ctd)
#'
#' @export
map_specificity_to_entrez <- function(ctd,annotLevel,specificity_species, return_ctd=F){
    # Check specificity_species
    if(!specificity_species %in% c("human","mouse")){stop("Specificity species must be either 'human' or 'mouse'")}
    
    int_all_hgnc_wtEntrez <- MAGMA.Celltyping::all_hgnc_wtEntrez
    colnames(int_all_hgnc_wtEntrez)[1] = "human.symbol"
    
    if(specificity_species=="mouse"){
        # Because sumstats use entrez genes & ctd uses gene symbols, match entrez-->symbols
        entrez_mgi = merge(int_all_hgnc_wtEntrez,
                           One2One::ortholog_data_Mouse_Human$orthologs_one2one[,2:3],
                           by="human.symbol")
        entrez_mgi = entrez_mgi[!is.na(entrez_mgi$entrezgene),]
        entrez_mgi = entrez_mgi[entrez_mgi$mouse.symbol %in% rownames(ctd[[annotLevel]]$specificity_quantiles),] 
        
        # Get the quantiles from ctd and put into correct format, using entrez symbols
        quantDat = as.matrix(ctd[[annotLevel]]$specificity_quantiles[entrez_mgi$mouse.symbol,])
        quantDat2 = suppressWarnings(data.frame(entrez=entrez_mgi$entrezgene,quantDat))
        quantDat2 = quantDat2[!duplicated(quantDat2$entrez),]
    }
    
    if(specificity_species=="human"){
        # Get the quantiles from ctd and put into correct format, using entrez symbols
        humanSymsPresent = as.character(int_all_hgnc_wtEntrez$human.symbol[int_all_hgnc_wtEntrez$human.symbol %in% rownames(ctd[[annotLevel]]$specificity_quantiles)])
        entrezTable = int_all_hgnc_wtEntrez[int_all_hgnc_wtEntrez$human.symbol %in% humanSymsPresent,]
        quantDat = as.matrix(ctd[[annotLevel]]$specificity_quantiles[as.character(entrezTable$human.symbol),])
        quantDat2 = suppressWarnings(data.frame(entrez=entrezTable$entrez,quantDat))
        quantDat2 = quantDat2[!duplicated(quantDat2$entrez),]
    }
    if(return_ctd){
        ctd[[annotLevel]]$quantDat2 <- quantDat2
        return(ctd)
    }else {return(quantDat2)} 
}
