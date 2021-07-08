
#' Plot enrichment results from \code{celltype_associations_pipeline}
#' 
#' @examples 
#' \dontrun{
#' library(MAGMA.Celltyping)
#' local_files <- import_magma_files(download_dir=".")
#' #' magma_dirs <- unique(dirname(local_files))
#' res <- celltype_associations_pipeline(ctd=ewceData::ctd(), 
#'                                       ctd_name="Zeisel2018", 
#'                                       magma_dirs=magma_dirs, 
#'                                       genome_ref_path="~/Downloads/g1000_eur/g1000_eur") 
#' merged_results <- gather_results(res)
#' heat <- results_heatmap(merged_results)
#' }
#' @export
results_heatmap <- function(merged_results, 
                            title=NULL,
                            x_lab=NULL,
                            fdr_thresh=.05,
                            facet_formula="EnrichmentMode ~ .",
                            x_var="Celltype",
                            y_var="GWAS",
                            fill_var="-log1p(FDR)",
                            scales = "free_y",
                            space = "fixed",
                            show_plot=T){
    library(ggplot2)
    if(!"GWAS" %in% colnames(merged_results)) merged_results$GWAS <- merged_results$dataset
    if(is.null(fdr_thresh)) fdr_thresh <- 1
    heat <- ggplot(data = subset(merged_results, FDR<fdr_thresh),
                   aes_string(x=x_var, y=y_var, fill= fill_var) ) +
        geom_tile(color="white") +
        facet_grid(facets = as.formula(facet_formula),
                   scales = scales, 
                   space = space) + 
        labs(title=title,
             subtitle=if(!is.null(fdr_thresh)) paste0("FDR < ",fdr_thresh) else NULL,
             x=x_lab) +
        scale_fill_gradient(low = "blue", high = "red", na.value = "white") + 
        theme_bw() +
        theme(axis.text.x = element_text(angle=45, hjust=1, size=7),
              strip.background = element_rect(fill="grey20"),
              strip.text = element_text(color = "white"),  
              panel.grid.minor = element_blank()) 
    if(show_plot) print(heat)
    return(heat)
}
