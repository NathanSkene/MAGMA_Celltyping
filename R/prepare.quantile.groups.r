#' Prepare quantile groups for each celltype based on specificity
#'
#' Quantile groups are stored in an extra matrix ('quantiles') in the returned CTD. This function also removes any genes
#' from the CTD data which are not 1:1 orthologs with the GWAS species.
#'
#' @param ctd Cell type data
#' @param specificity_species Species name relevant to the cell type data, i.e. "mouse" or "human"
#' @param gwas_species Species name relevant to the GWAS data, in almost all cases this will be "human"
#' @param numberOfBins How many bins should specificity be divided into?
#'
#' @return The ctd with additional quantiles matrix
#'
#' @examples
#' \dontrun{
#' ctd <- prepare.quantile.groups(ewceData::ctd())
#' } 
#' @import tibble
#' @import dplyr
#' @import EWCE
#' @importFrom One2One load.homologs
#' @importFrom One2One analyse.orthology
#' @export
prepare.quantile.groups <- function(ctd,
                                    specificity_species="mouse",
                                    gwas_species="human",
                                    numberOfBins=40){
    ### Account for changes in naming conventions across EWCE versions ####
    bin_columns_func <- if(packageVersion("EWCE")>="1.0.0"){EWCE::bin_columns_into_quantiles}else EWCE::bin.columns.into.quantiles
    bin_specificity_func <- if(packageVersion("EWCE")>="1.0.0"){EWCE::bin_specificity_into_quantiles}else EWCE::bin.specificity.into.quantiles
    
    #library(tibble)
    # First drop all genes without 1:1 homologs
    if(specificity_species != gwas_species){
        print("Dropping all genes that do not have 1:1 homologs between the two species")
        allHomologs = One2One::load.homologs()
        ortholog_data = One2One::analyse.orthology(specificity_species,gwas_species,allHomologs)$orthologs_one2one
        ctd = lapply(ctd,filter_by_orthologs,one2one_ortholog_symbols = ortholog_data[,2])
    }
    # Quantiles will be stored within the CTD as 'quantiles'
    # - Function below 'bin.columns.into.quantiles' has been moved to EWCE
    # bin.columns.into.quantiles <- function(spcValues,numberOfBins){
    #     quantileValues = rep(0,length(spcValues))
    #     quantileValues[spcValues>0] = as.numeric(cut(spcValues[spcValues>0], 
    #                                                  breaks=unique(quantile(spcValues[spcValues>0], probs=seq(0,1, by=1/numberOfBins), na.rm=TRUE)), 
    #                                                  include.lowest=TRUE))
    #     return(quantileValues)
    # }
    normalise.mean.exp <- function(spcMatrix){
        spcMatrix$mean_exp <- as(spcMatrix$mean_exp,"matrix")
        spcMatrix$specificity_quantiles <- as(spcMatrix$specificity_quantiles,"matrix")
        spcMatrix$linear_normalised_mean_exp = t(t(spcMatrix$mean_exp) *
                                                     (1/colSums(spcMatrix$mean_exp)))
        return(spcMatrix)
    }
    
    bin.expression.into.quantiles <- function(spcMatrix,numberOfBins){
        spcMatrix$expr_quantiles = as.matrix(apply(spcMatrix$linear_normalised_mean_exp,2,FUN=bin_columns_func,numberOfBins=numberOfBins))
        rownames(spcMatrix$expr_quantiles) = rownames(spcMatrix$linear_normalised_mean_exp)
        return(spcMatrix)
    }
    use.distance.to.add.expression.level.info <- function(spcMatrix){
        spcMatrix$specificity <- as(spcMatrix$specificity,"matrix")
        spcMatrix$spec_dist = spcMatrix$specificity
        for (ct in colnames(spcMatrix$expr_quantiles)) {
            resTab = data.frame(spec = spcMatrix$specificity_quantiles[,
                                                                       ct], exp = spcMatrix$expr_quantiles[, ct], gene = rownames(spcMatrix$linear_normalised_mean_exp))
            resTab$dist = sqrt((max(resTab$spec) - resTab$spec)^2 +
                                   (max(resTab$exp) - resTab$exp)^2)
            spcMatrix$spec_dist[, ct] = resTab$dist
        }
        spcMatrix$spec_dist <- as(spcMatrix$spec_dist,"matrix")
        spcMatrix$spec_dist = max(spcMatrix$spec_dist) - spcMatrix$spec_dist
        return(spcMatrix)
    }
    bin.specificityDistance.into.quantiles <- function(spcMatrix){  
        spcMatrix$specDist_quantiles = apply(spcMatrix$spec_dist, 2, FUN=bin_columns_func)
        rownames(spcMatrix$specDist_quantiles) = rownames(spcMatrix$spec_dist)
        spcMatrix$specDist_quantiles <- as(spcMatrix$specDist_quantiles,"matrix")
        return(spcMatrix) 
    }     
    ctd = lapply(ctd, normalise.mean.exp)
    ctd = lapply(ctd, bin_specificity_func, numberOfBins=numberOfBins)
    #ctd = lapply(ctd,bin.expression.into.quantiles)
    ctd = lapply(ctd, use.distance.to.add.expression.level.info)
    ctd = lapply(ctd, bin.specificityDistance.into.quantiles)
    return(ctd)
}
