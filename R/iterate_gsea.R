#' Iterate gene set enrichment analysis across all traits and celltypes
#' 
#' @param xmat gene x trait matrix.
#' @param ymat gene x celltype matrix.
#' @param correction_method Multiple-testing correction method to be passed to \code{stats::p.adjust}.
#' @param qvalue_thresh q.value threshold to use when report significant results summary.
#' @param x_quantiles The number of quantiles to bin \code{ymat} data into.
#' @param y_quantiles The number of quantiles to bin \code{ymat} data into.
#' @param nCores Number of cores to use in parallel. Will optimize if \code{NULL}.
#' @export
iterate_gsea <- function(xmat, 
                         ymat,
                         correction_method="BH",
                         qvalue_thresh=.05,
                         x_quantiles=10,
                         y_quantiles=10,
                         use_quantiles=10,
                         nCores=1){
    gene_intersect <- intersect(rownames(xmat), rownames(ymat))
    message(length(gene_intersect)," intersecting genes between GWAS and CTD matrices.")
    ### Run lm  for all celltypes against this trait
    message("Running ",formatC(ncol(xmat)*ncol(ymat), big.mark=","),
            " tests: ",formatC(ncol(xmat), big.mark=",")," traits x ",
            formatC(ncol(ymat), big.mark=","), " celltypes.") 
    
    gsea_res <- parallel::mclapply(1:ncol(xmat), function(i){
        tt <- colnames(xmat)[i]
        EWCE:::message_parallel(" - ",tt,": (",i,"/",ncol(xmat),")")
        lapply(colnames(ymat), function(ct){ 
            # EWCE:::message_parallel(" - ",ct)  
            dat <- data.frame(trait= cut(xmat[gene_intersect,tt], breaks=x_quantiles, labels = 1:x_quantiles),
                              celltype = cut(ymat[gene_intersect,ct], breaks=y_quantiles, labels = 1:y_quantiles),
                              row.names = gene_intersect)  
            res <- GeneOverlap::newGeneOverlap(listA =rownames(subset(dat, trait %in% use_quantiles)),
                                               listB = rownames(subset(dat, celltype %in% use_quantiles))) %>%
                GeneOverlap::testGeneOverlap() 
            res_df <- data.frame(term=ct,
                                 trait_genes=length(res@listA),
                                 celltype_genes=length(res@listB),
                                 intersection=length(res@intersection),
                                 union=length(res@union),
                                 genome.size=length(res@genome.size),
                                 odds.ratio=res@odds.ratio,
                                 Jaccard=res@Jaccard,
                                 p.value=res@pval) 
            return(res_df)
        }) %>% data.table::rbindlist() 
    }, mc.cores = nCores) %>% 
        `names<-`(colnames(xmat)) %>%
        data.table::rbindlist(idcol = "trait")
    
    ### Multiple-testing correction
    gsea_res <- gsea_res %>% 
        dplyr::mutate(q.value=p.adjust(p = p.value, method = correction_method))
    ### Filter only sig results
    sig_res <- gsea_res %>%
        subset(q.value<qvalue_thresh)
    message("\n",formatC(nrow(sig_res), big.mark=",")," significant results @ ",correction_method,"<",qvalue_thresh)
    ### Return FULL results (not just sig)
    return(gsea_res)
} 
