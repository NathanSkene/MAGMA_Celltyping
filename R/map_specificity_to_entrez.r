#' Map specificity to entrez
#'
#' Convenience function used in 'create_gene_covar_file()'
#'
#' @param genesOutFile The output of the second call to MAGMA (performed in the map.snps.to.genes function)
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
map_specificity_to_entrez <- function(genesOutFile,ctd,annotLevel,specificity_species){
    # Check specificity_species
    if(specificity_species %in% c("human","mouse")){
        
    data(all_hgnc_wtEntrez)
    colnames(all_hgnc_wtEntrez)[1] = "human.symbol"
    
    if(specificity_species=="mouse"){
        data(ortholog_data_Mouse_Human)
        
        # Because sumstats use entrez genes & ctd uses gene symbols, match entrez-->symbols
        entrez_mgi = merge(all_hgnc_wtEntrez,ortholog_data_Mouse_Human$orthologs_one2one[,2:3],by="human.symbol")
        entrez_mgi = entrez_mgi[!is.na(entrez_mgi$entrezgene),]
        entrez_mgi = entrez_mgi[entrez_mgi$mouse.symbol %in% rownames(ctd[[annotLevel]]$specificity_quantiles),] 
        
        # Get the quantiles from ctd and put into correct format, using entrez symbols
        quantDat = ctd[[annotLevel]]$specificity_quantiles[entrez_mgi$mouse.symbol,]
        quantDat2 = suppressWarnings(data.frame(entrez=entrez_mgi$entrezgene,quantDat))
        quantDat2 = quantDat2[!duplicated(quantDat2$entrez),]
    }
    
    if(specificity_species=="human"){
        # Get the quantiles from ctd and put into correct format, using entrez symbols
        humanSymsPresent = as.character(all_hgnc_wtEntrez$human.symbol[all_hgnc_wtEntrez$human.symbol %in% rownames(ctd[[annotLevel]]$specificity_quantiles)])
        entrezTable = all_hgnc_wtEntrez[all_hgnc_wtEntrez$human.symbol %in% humanSymsPresent,]
        quantDat = ctd[[annotLevel]]$specificity_quantiles[as.character(entrezTable$human.symbol),]
        quantDat2 = suppressWarnings(data.frame(entrez=entrezTable$entrez,quantDat))
        quantDat2 = quantDat2[!duplicated(quantDat2$entrez),]
    }
    } else {
        
        print("generating ortholog data for specific specificity species.")
        data(all_hgnc_wtEntrez)
        colnames(all_hgnc_wtEntrez)[1] = "human.symbol"

        library(One2One)

        # Download and format the homolog data from MGI
        allHomologs = load.homologs()

        # Get data on orthology between the two species
        species1 = "human" #this needs to be the gwas_species
        species2 = specificity_species
        ortholog_data = analyse.orthology(species1,species2,allHomologs)

        # Because sumstats use entrez genes & ctd uses gene symbols, match entrez-->symbols
        entrez_mgi = merge(all_hgnc_wtEntrez,ortholog_data$orthologs_one2one[,2:3],by="human.symbol")
        entrez_mgi = entrez_mgi[!is.na(entrez_mgi$entrezgene),]
        #the third column in the entrez_mgi data.frame is 'specificity_species'.symbol
        entrez_mgi = entrez_mgi[entrez_mgi[[3]] %in% rownames(ctd[[annotLevel]]$specificity_quantiles),] 

        # Get the quantiles from ctd and put into correct format, using entrez symbols
        quantDat = ctd[[annotLevel]]$specificity_quantiles[entrez_mgi[[3]],]
        quantDat2 = suppressWarnings(data.frame(entrez=entrez_mgi$entrezgene,quantDat))
        quantDat2 = quantDat2[!duplicated(quantDat2$entrez),]
    }
    return(quantDat2)
}
