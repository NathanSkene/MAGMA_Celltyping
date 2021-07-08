





#' Install the MAGMA command line tool
#' 
#' @source [MAGMA website](https://ctg.cncr.nl/software/magma)
#' @source [MAGMA.celltyping documentation](https://github.com/NathanSkene/MAGMA_Celltyping)
install_magma <- function(dest_dir="~"){ 
    
    get_os <- function () {
        OS=""
        switch(Sys.info()[['sysname']],
               Windows= {OS="Windows"},
               Linux  = {OS="Linux"},
               Darwin = {OS="Mac"})
        return(OS) 
    }
    check_magma <- ""
    try({
        check_magma <- system("magma", intern=T)
    })
   
    
    if(check_magma[1]=="No arguments specified. Please consult manual for usage instructions."){
        message("MAGMA already installed.")
    }else {
        message("MAGMA not installed. Downloading...")
        magma_url <- switch(get_os(),  
                            Mac="https://ctg.cncr.nl/software/MAGMA/prog/magma_v1.08a_mac.zip",
                            Linux="https://ctg.cncr.nl/software/MAGMA/prog/magma_v1.08a.zip",
                            NULL="https://ctg.cncr.nl/software/MAGMA/prog/magma_v1.08a.zip")
        destfile <- file.path(dest_dir,basename(magma_url))
        destpath <- gsub(".zip","",destfile)
        download.file(magma_url,
                      destfile = destfile) 
        unzip(destfile, junkpaths = T,
              exdir = gsub(".zip","",destfile),
              overwrite = T)
        file.remove(destfile)
        message("MAGMA installation complete.")
        message("MAGMA path:\n",destpath)
        # Create a symlink to the actualy magma executable
        R.utils::createLink(link="/usr/local/bin/magma",
                            target=file.path(destpath,"magma"))
    }
}
