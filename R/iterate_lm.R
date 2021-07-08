#' Iterate linear regression across all traits and celltypes
#' 
#' @param xmat gene x trait matrix.
#' @param ymat gene x celltype matrix.
#' @param correction_method Multiple-testing correction method to be passed to \code{stats::p.adjust}.
#' @param qvalue_thresh q.value threshold to use when report significant results summary.
#' @param y_quantiles The number of quantiles to bin \code{ymat} data into.
#' @param nCores Number of cores to use in parallel. Will optimize if \code{NULL}.  
#' @export
iterate_lm <- function(xmat, 
                       ymat,
                       correction_method="BH",
                       qvalue_thresh=.05,
                       y_quantiles=NULL,
                       nCores=NULL){
    if(is.null(nCores)) nCores <- EWCE:::assign_cores(worker_cores = nCores)$worker_cores
    gene_intersect <- intersect(rownames(xmat), rownames(ymat))
    message(length(gene_intersect)," intersecting genes between xmat and ymat")
    ### Run lm  for all celltypes against this trait
    message("Running ",formatC(ncol(xmat)*ncol(ymat), big.mark=","),
            " tests: ",formatC(ncol(xmat), big.mark=",")," traits x ",
            formatC(ncol(ymat), big.mark=","), " celltypes.") 
    
    lm_res <- parallel::mclapply(1:ncol(xmat), function(i){
        tt <- colnames(xmat)[i]
        EWCE:::message_parallel(" - ",tt,": (",i,"/",ncol(xmat),")")
        lapply(colnames(ymat), function(ct){ 
            # EWCE:::message_parallel(" - ",ct) 
            if(!is.null(y_quantiles)){
                lm_dat <- data.frame(trait=xmat[gene_intersect,tt],
                                     celltype=cut(ymat[gene_intersect,ct], breaks=y_quantiles, labels = 1:y_quantiles)) 
            } else {
                lm_dat <- data.frame(trait=xmat[gene_intersect,tt],
                                     celltype=ymat[gene_intersect,ct])
            }
            mod <- stats::lm(data=lm_dat, 
                             formula = trait ~ celltype)
            res_df <- subset(broom::tidy(mod), term!="(Intercept)")
            res_df$term <- ct
            return(res_df)
        }) %>% data.table::rbindlist() 
    }, mc.cores = nCores) %>% 
        `names<-`(colnames(xmat)) %>%
        data.table::rbindlist(idcol = "trait")
    
    ### Multiple-testing correction
    lm_res <- lm_res %>% 
        dplyr::mutate(q.value=p.adjust(p = p.value, method = correction_method))
    ### Filter only sig results
    sig_res <- lm_res %>%
        subset(q.value<qvalue_thresh)
    message("\n",formatC(nrow(sig_res), big.mark=",")," significant results @ ",correction_method,"<",qvalue_thresh)
    ### Return FULL results (not just sig)
    return(lm_res)
} 