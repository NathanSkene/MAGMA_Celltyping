
create_genesets <- function(res_input,
                            n_genes=100,
                            spec_deciles=10,
                            human.symbol=T){
    symbol_species <- if(human.symbol)"human.symbol"else"mouse.symbol"
    message("+ Returning ",symbol_species)
    res_input <- subset(res_input, specificity_decile %in% spec_deciles)
    genesets <- lapply(unique(res_input$Celltype_id), function(x){
        g <- (subset(res_input, Celltype_id==x) %>%
                  dplyr::slice_max(ADJ_ZSTAT, n = n_genes, with_ties = F))[[symbol_species]]
        message("   ",x," : ",length(g)," genes")
        return(g)
    } ) %>% `names<-`(unique(res_input$Celltype_id))
    return(genesets)
}