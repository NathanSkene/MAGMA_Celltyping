



#' Import remote MAGMA files
#' 
#' Import preprocessed MAGMA GWAS files from remote
#'  \href{https://github.com/neurogenomics/MAGMA_Files}{GitHub repository}.
#' Each MAGMA file will be put into its own folder of the same name as the file. 
#' \bold{WARNING}: 
#' 
#' @param download_dir Parent folder where you want to save the MAGMA files.
#' @param overwrite If the files have already been downloaded in the specified directory, 
#' these downloads will be skipped. Set \code{overwrite=T} to force them to be re-downloaded.
#' @param nThread Number of threads to parallelize downloading.
#' 
#' @return 
#' Paths to where the MAGMA files have been downloaded locally.
#' 
#' @examples 
#' local_files <- import_magma_files(download_dir=".")
#' @export 
import_magma_files <- function(download_dir=tempdir(),
                               file_types=c(".genes.raw",".genes.out"),
                               overwrite=F,
                               nThread=parallel::detectCores()-2){ 
    magma_files <- GITHUB.list_files(creator="neurogenomics",
                                     repo="MAGMA_Files",
                                     branch="main",
                                     query = paste(file_types, collapse = "|"),
                                     return_download_api = T)
    local_files <- GITHUB.download_files(filelist = magma_files, 
                                          download_dir = file.path(download_dir,"MAGMA_Files"),
                                          nThread = nThread, 
                                          overwrite = overwrite)
    return(local_files) 
}





GITHUB.list_files <- function(creator="neurogenomics",
                              repo="MAGMA_Files",
                              branch=c("main","master"),
                              query=NULL,
                              return_download_api=T,
                              verbose=T){
    repo_api <- file.path("https://api.github.com/repos",creator,repo,
                          paste0("git/trees/",branch[1],"?recursive=1"))  
    req <- httr::GET(repo_api)
    httr::stop_for_status(req)
    filelist <- unlist(lapply(httr::content(req)$tree, "[", "path"), use.names = F)
    printer(paste(length(filelist),"files found in GitHub repo:", file.path(creator,repo)),v=verbose)
    if(!is.null(query)){
        # query_string <- "*Nalls23andMe_2019.*UKB.multi_finemap.csv.gz"
        bool <- grepl(query, filelist)
        filelist <- filelist[bool]
        printer(paste(length(filelist),"files found matching query."),v=verbose)
    }
    if(return_download_api){
        filelist <- file.path("https://github.com",creator,repo,"raw",branch,filelist)
    }
    return(filelist)
}



GITHUB.download_files <- function(filelist,
                                  download_dir=tempdir(), 
                                  overwrite=F,
                                  nThread=parallel::detectCores()-2,
                                  verbose=T){ 
    printer("+ Downloading",length(filelist),"files...",v=verbose)
    local_files <- unlist(parallel::mclapply(filelist, function(x){
        print(paste("Downloading",x))
        branch <- stringr::str_split(string = x, pattern = "/")[[1]][7]
        folder_structure <- paste(stringr::str_split(string = x, pattern = "/")[[1]][-c(1:7)], collapse="/")
        destfile <- file.path(download_dir, folder_structure)
        dir.create(dirname(destfile), showWarnings = F, recursive = T)
        if(!file.exists(destfile) & overwrite==F) download.file(url = x, destfile=destfile)
        return(destfile)
    }, mc.cores = nThread))
    return(local_files)
}

