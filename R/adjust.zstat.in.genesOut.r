#' Adjust MAGMA Z-statistic from .genes.out files
#'
#' Used when you want to directly analyse the gene level z-scores corrected for gene length etc
#'
#' @param ctd Cell type data structure
#' @param magma_GenesOut_file A MAGMA .genes.out file
#' @param sctSpecies Either 'human' or 'mouse'
#'
#' @examples
#' myGenesOut = tempfile()
#' data.table::fwrite(x=MAGMA.Celltyping::genesOut,sep="\t",file=myGenesOut)
#' magmaGenesOut = adjust.zstat.in.genesOut(EWCE::ctd,magma_GenesOut_file=myGenesOut,
#' sctSpecies="mouse")
#'
#' @export
#' @importFrom stats lm
#' @importFrom stats p.adjust
#' @importFrom utils read.table
adjust.zstat.in.genesOut <- function(ctd,
                                     magma_GenesOut_file=NA,
                                     sctSpecies="mouse"){
    allGenes = rownames(ctd[[1]]$specificity)
    
    if(sctSpecies=="mouse"){
        # Get mouse-->human othologs with human entrez
        orth = One2One::ortholog_data_Mouse_Human$orthologs_one2one[,2:3]
        hgnc2entrez = MAGMA.Celltyping::all_hgnc_wtEntrez
        colnames(hgnc2entrez)=c("human.symbol","entrez")
        orth2 = merge(orth,hgnc2entrez,by="human.symbol")
    }
    
    # Load the MAGMA data
    magma = utils::read.table(magma_GenesOut_file,stringsAsFactors = FALSE,header=TRUE)
    magma$entrez = magma$GENE
    if(sctSpecies=="mouse"){
        magma = merge(magma, orth2,by="entrez")
    }else if(sctSpecies=="human"){
        magma = merge(magma, MAGMA.Celltyping::hgnc2entrez,by="entrez")
    }
    
    magma = magma[order(magma$P),]
    if(sctSpecies=="mouse"){
        magma = magma[magma$mouse.symbol %in% allGenes,]
        magma = magma[!duplicated(magma$mouse.symbol),]
    }else if(sctSpecies=="human"){
        magma = magma[magma$hgnc_symbol %in% allGenes,]
        magma = magma[!duplicated(magma$hgnc_symbol),]
    }
    magma$Q = stats::p.adjust(magma$P,method="bonferroni")
    magma$logNSNPS=log(magma$NSNPS)
    magma$logNPARAM=log(magma$NPARAM)
    magma = magma[!(magma$CHR==6 & magma$START>=25000000 & magma$STOP<=34000000),] # DROP MHC: chr6, 25-34 mb
    magma$GENELEN = abs(magma$STOP-magma$START)
    magma$logGENELEN = log(magma$GENELEN)
    
    # Regress out effects of NSNPS and NPARAM (see 'boxplots_by_decile.r' and the section on downsampling for info)
    #--- NSNPS only really has mjaor effects (i.e. zscore+2) when a gene has ~10000 SNPS
    mod = stats::lm(ZSTAT~NSNPS+logNSNPS+NPARAM+logNPARAM+GENELEN+logGENELEN,data=magma)
    magma$ADJ_ZSTAT = magma$ZSTAT - (magma$NSNPS*mod$coefficients[2] + magma$logNSNPS*mod$coefficients[3] + 
                                         magma$NPARAM*mod$coefficients[4] + magma$logNPARAM*mod$coefficients[5] + 
                                         magma$GENELEN*mod$coefficients[6] + magma$logGENELEN*mod$coefficients[7]    )
    #magma = arrange(magma,desc(magma$ADJ_ZSTAT))
    #magma = magma[order(magma$ADJ_ZSTAT,decreasing=TRUE),]
    return(magma)
}