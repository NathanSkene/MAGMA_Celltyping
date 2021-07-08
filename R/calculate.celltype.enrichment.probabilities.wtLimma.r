#' Calculate celltype enrichments without MAGMA using adjusted MAGMA Z-statistic from .genes.out files
#'
#' @param magmaAdjZ Output from adjust.zstat.in.genesOut()
#' @param ctd CellTypeData object
#' @param thresh A threshold on low specificity values (used to drop genes)
#' @param sctSpecies Either 'human' or 'mouse'
#' @param annotLevel CellTypeData level to test
#' #'
#' @examples
#' library(MAGMA.Celltyping)
#' library(ewceData)
#' ctd <- ewceData::ctd()
#'
#'  # The package stores an example genesOut file, so save this to a tempfile
#' myGenesOut = tempfile()
#' data.table::fwrite(MAGMA.Celltyping::genesOut,sep="\t",file=myGenesOut)
#' magmaAdjZ = adjust.zstat.in.genesOut(ctd, magma_GenesOut_file=myGenesOut,sctSpecies="mouse")
#' ps <- calculate.celltype.enrichment.probabilities.wtLimma(magmaAdjZ=magmaAdjZ, ctd=ctd)
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr rename
#' @importFrom dplyr filter
#' @importFrom limma lmFit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' @importFrom stats quantile
#' @importFrom rlang .data
calculate.celltype.enrichment.probabilities.wtLimma <- function(magmaAdjZ,
                                                                ctd,
                                                                thresh=0.0001,
                                                                sctSpecies="mouse",
                                                                annotLevel=1,
                                                                celltypes=NULL,
                                                                return_all=F){ 
    if("hgnc_symbol" %in% colnames(magmaAdjZ)){
        magmaAdjZ = magmaAdjZ %>% dplyr::rename(human.symbol=.data$hgnc_symbol)
    }
    
    # First get names of all cell types
    #### Ensure cell-type names are processed the same as in the MAGMA.Celltyping pipeline
    message("+ Preparing specificity matrix")
    spec <- ctd[[annotLevel]]$specificity
    og_colnames <- colnames(spec)
    spec <- as.matrix(data.frame(spec))
    allCellTypes = colnames(spec)
    celltype_dict <- setNames(og_colnames, allCellTypes)
  
    if(!is.null(celltypes)){ 
        selected_celltypes <- unname(get_celltype_dict(all_df = data.frame(Celltype=celltypes, dummy=1))) 
        allCellTypes <- allCellTypes[tolower(allCellTypes) %in% tolower(selected_celltypes)]
        message("+ ",length(allCellTypes)," cell-types selected.") 
    }
    
    
    # Initialise variables
    ps = coef = rep(0,length(allCellTypes))
    names(ps) = names(coef) = allCellTypes
    
    # Loop over each celltype testing for enrichment
    count=0
    res_all <- list()
    input_all <- list()
    for(ct1 in allCellTypes){
        count=count+1
        print(ct1)
        if(sctSpecies=="mouse"){
            mgiS  = magmaAdjZ$mouse.symbol
        }else{
            mgiS = magmaAdjZ$human.symbol
        }
        props = spec[mgiS,ct1]
        notExp = rep(0,length(props))
        
        # Drop any genes with expression below threshold
        mgiS = mgiS[props>thresh]
        props = props[props>thresh]
        notExp[props<thresh]=-1
        
        # Determine which expression decile genes fall into
        quantiles = stats::quantile(props[props>thresh],probs=seq(from=0,to=1,by=0.1))
        perc  = as.numeric(cut(props,quantiles,include.lowest=TRUE))
        perc[is.na(perc)] = 0
        
        # Merge decile groups with MAGMA zscores
        if(sctSpecies=="mouse"){
            geneGroups = data.frame(mgi_symbol=mgiS,proportion=props,percentile=perc)
            magma_with_ct1 = geneGroups %>%
                dplyr::rename(mouse.symbol=.data$mgi_symbol) %>% merge(magmaAdjZ,by="mouse.symbol") %>%
                dplyr::filter(.data$percentile>=0)
        }else{
            geneGroups = data.frame(hgnc_symbol=mgiS,proportion=props,percentile=perc)
            magma_with_ct1 = geneGroups %>%
                dplyr::rename(human.symbol=.data$hgnc_symbol) %>% merge(magmaAdjZ,by="human.symbol") %>%
                dplyr::filter(.data$percentile>=0)            
        }
        input_all[[ct1]] <- magma_with_ct1
        
        # Fit a linear model and get p-value and coefficient (slope)
        expMat = matrix(0,nrow=1,ncol=dim(magma_with_ct1)[1])
        expMat[1,] = magma_with_ct1$ADJ_ZSTAT
        rownames(expMat) = c("actual")
        colnames(expMat) = magma_with_ct1$mouse.symbol
        dmat <- stats::model.matrix(~ magma_with_ct1$percentile)
        fit <- limma::lmFit(expMat, dmat)
        fit2 <- limma::eBayes(fit)
        res <- limma::topTable(fit2,number=1)
        res[res$logFC<0,"P.Value"]=1
        res[res$logFC>=0,"P.Value"]=res[res$logFC>=0,"P.Value"]/2   
        
        # Convert p-value to one-sided
        ps[count] = res$P.Value
        res_all[[ct1]] <- res
    }
    
    if(return_all){ 
        #### res df
        res_all_df <- data.table::rbindlist(res_all, idcol = "Celltype_id") 
        res_all_df <- res_all_df %>%
            tibble::add_column(res_all_df, 
                               Celltype=celltype_dict[res_all_df$Celltype_id],
                               .after = "Celltype_id")  %>%
            dplyr::mutate(ctd_level=annotLevel)
        #### input df
        input_all_df <- data.table::rbindlist(input_all,idcol = "Celltype_id") %>%
            dplyr::rename(specificity_proportion=proportion, 
                          specificity_decile=percentile) 
        input_all_df <- input_all_df %>%
            tibble::add_column(input_all_df, 
                               Celltype=celltype_dict[input_all_df$Celltype_id],
                               .after = "Celltype_id") %>%
            dplyr::mutate(ctd_level=annotLevel)
        
        return(list(input=input_all_df,
                    results=res_all_df
        ))
    }else {return(ps)}
}


get_celltype_dict <- function(all_df){
    ##### Celltype names get changed during MAGMA_Celltyping. Translate here
    tmp <- all_df
    tmp$dummy <- 1
    tmp <- t(unique(tmp[,c("Celltype","dummy")]))
    colnames(tmp) <- tmp[1,]
    tmp <- data.frame(tmp)
    celltype_dict <- setNames(colnames(tmp),tmp[1,])
    return(celltype_dict)
}
