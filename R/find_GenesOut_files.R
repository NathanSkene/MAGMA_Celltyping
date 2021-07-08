find_GenesOut_files <- function(sig_res, 
                                root_dir="raw_data/MAGMA/MAGMA_Files"){
    matches <- apply(stringr::str_split(sig_res$GCOV_FILE,"[.]", n = 5, simplify = T)[,1:4], 1, 
                     paste, collapse=".") %>% unique()    
    magma_GenesOut_file <- list.files(root_dir,paste0(matches,".genes.out"), 
                                      recursive = T, full.names = T)
    message(length(magma_GenesOut_file)," .genes.out file(s) found.")
    return(magma_GenesOut_file)
}


