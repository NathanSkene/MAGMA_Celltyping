#' Check that sumstats has correct columns and that they are in the correct order for MAGMA and LDSC
#'
#' @return col_headers The new column headers for the sumstats file
#'
#' @examples
#' format_sumstats_for_magma(path)
#'
#' @import data.table
#' @import stringr
#' @export
format_sumstats_for_magma <- function(path, N){
  
  # Checking if the file exists should happen first
  if (!file.exists(path)) {stop("Path to GWAS sumstats is not valid")}
  
  # This almost surely modifies the file (since most sumstats from different studies are differently formatted), so it makes more sense to just make a temporary file <tmp>, and return the address of the temp
  sumstats_file <- readLines(path)
  tmp <- tempfile()
  writeLines(sumstats_file, con=tmp)
  path <- tmp
  
  # Ensure that tabs separate rows
  row_of_data <- strsplit(sumstats_file[2], "\t")[[1]]
  if (length(row_of_data) == 1) {
    if (grep(" ", row_of_data) == 1) {
      print("WARNING: This GWAS sumstat file has space field separators instead of tabs (unusual, not proper input for MAGMA). Temp file with corrected FS created and used instead.")
      sumstats_file <- gsub(pattern = " ", replace = "\t", x = before)
    }
  }
  
  sumstats_file[1] = standardise.sumstats.column.headers.crossplatform(sumstats_file[1])
  col_headers = sumstats_file[1]
  col_headers = strsplit(col_headers, "\t")[[1]]
  
  # Check if there are CHR and BP columns
  if(!sum(c("SNP","BP") %in% col_headers)==2){
    # If not, see if there is a column storing the data in a format like: CHR:BP:A2:A1 or just CHR:BP
    # - UKBB data from Ben Neale has a Variant column with CHR:POS:REF:ALT where ALT allele is the effect allele in the model [NB: the ALT allele is NOT always the minor allele]
    # -- For input to LDSC, A1 is effect allele, A2 is non-effect allele
    # - DIAGRAM diabetes data has a Chr:Position column with CHR:BP
    # - BMI adjusted for smoking has markername with CHR:BP (with the chromosome name having 'chr' preceeding)
    # - Agression [EAGLE] just doesn't have any CHR or BP data
    print("Summary statistics file does not have obvious CHR or BP columns. Checking to see if they are joined in another column")
    
    # Obtain a row of the actual data
    row_of_data <- strsplit(sumstats_file[2], "\t")[[1]]
    
    # Check if there is a column of data with CHR:BP:A2:A1 format
    fourStepCol = grep(".*:.*:\\w:\\w",row_of_data)
    if(length(fourStepCol)){
      # Convert the ':' into '\t'
      sumstats_file <- gsub(pattern = ":", replace = "\t", x = sumstats_file)
      # Replace the column name with four names
      curColName = col_headers[fourStepCol]
      # Write the new column headers to file
      first_line = paste(col_headers, collapse = "\t")
      new_first_line = gsub(sprintf("^%s\\t|\\t%s\\t|\\t%s$",curColName,curColName,curColName), "CHR\tBP\tA2\tA1\t", paste(col_headers, collapse = "\t"))
      sumstats_file[1]=new_first_line
      col_headers = strsplit(new_first_line, "\t")[[1]]
      print(sprintf("Column %s has been replaced with CHR BP A2 A1", curColName))
      print(col_headers)
      row_of_data <- strsplit(sumstats_file[2], "\t")[[1]]
    }
    
    # Check if there is a column of data with CHR:BP format
    twoStepCol = grep(".*:.*", row_of_data)
    if (length(twoStepCol)) {
      # Convert the ':' into '\t'
      sumstats_file <- gsub(pattern = ":", replace = "\t", x = sumstats_file)
      # Replace the column name with four names
      curColName = col_headers[twoStepCol]
      # Write the new column headers to file
      first_line = paste(col_headers,collapse = "\t")
      new_first_line = gsub(curColName,"CHR\tBP",paste(col_headers,collapse = "\t"))
      sumstats_file[1] <- new_first_line
      col_headers = strsplit(new_first_line,"\t")[[1]]
      print(sprintf("Column %s has been replaced with CHR BP",curColName))
      print(col_headers)
      row_of_data <- strsplit(sumstats_file[2], "\t")[[1]]
    }
    
    # Restandardise in case the joined column headers were unusual
    sumstats_file[1] = standardise.sumstats.column.headers.crossplatform(sumstats_file[1])
    col_headers = strsplit(sumstats_file[1], "\t")[[1]]
  }
  
  # If SNP is present... BUT not CHR or BP then need to find the relevant locations
  rows_of_data <- c(sumstats_file[1], sumstats_file[2]); col_headers = strsplit(rows_of_data[1], "\t")[[1]]; writeLines(sumstats_file, con = path)
  if(sum(c("CHR","BP") %in% col_headers)==0 & sum("SNP" %in% col_headers)==1){
    #library(data.table)
    #sumstats = fread(path)
    #SNP_LOC_DATA = load_snp_loc_data()
    #SNP_LOC_DATA_2 = SNP_LOC_DATA[SNP_LOC_DATA$Build=="GRCh37",1:3]
    #sumstats2 = merge(sumstats,SNP_LOC_DATA_2,by="SNP")
    #sumstats3 = data.frame(sumstats2)[,c("SNP","CHR","BP",setdiff(colnames(sumstats2),c("SNP","CHR","BP")))]
    #fwrite(sumstats3,file=path,sep="\t"); sumstats_file <- readLines(path)
    stop("I've blocked this function because I've not tested it since replacing SNP_LOC_DATA. You should test it works manually. Let me know if you test it!")
  }
  
  # If CHR and BP are present... BUT not SNP then need to find the relevant SNP ids
  rows_of_data <- c(sumstats_file[1], sumstats_file[2]); col_headers = strsplit(rows_of_data[1], "\t")[[1]]
  if(sum(c("CHR","BP") %in% col_headers)==2 & sum("SNP" %in% col_headers)==0){
    print("There is no SNP column found within the data. It must be inferred from CHR and BP information.")
    
    genomebuild <- as.numeric(readline("Which genome build is the data from? 1 for GRCh37, 2 for GRCh38... "))
    if(!genomebuild %in% c(1,2)){stop("Genome build must be entered as either 1 (for GRCh37) or 2 (for GRCh38)")}
    
    SNP_LOC_DATA = load_snp_loc_data()
    if(genomebuild==1){genomebuild="GRCh37"}else{genomebuild="GRCh38"}
    snpLocDat = SNP_LOC_DATA[SNP_LOC_DATA$Build==genomebuild,][,-4]
    library(data.table)
    sumstats = fread(path)
    sumstats$CHR = as.factor(sumstats$CHR)
    if(length(grep("chr",sumstats$CHR[1]))!=0){sumstats$CHR = gsub("chr","",sumstats$CHR)}
    sumstats2 = merge(sumstats,snpLocDat,by=c("CHR","BP"))
    # Remove any rows where P is NaN
    sumstats2 = sumstats2[!is.nan(sumstats2$P),]
    fwrite(sumstats2,file=path,sep="\t")
    sumstats_file <- readLines(path)
  }
  
  # Check that all the vital columns are present
  rows_of_data <- c(sumstats_file[1], sumstats_file[2]); col_headers = strsplit(rows_of_data[1], "\t")[[1]]
  for(key_column in c("SNP","CHR","BP","P","A1","A2")){
    code_example = "sed -i '' '1s/p_value/P/' IQ.Sniekers.2017.txt"
    if(!key_column %in% col_headers){
      print("Header of file:")
      print(rows_of_data)
      stop(sprintf("cannot find a %s column in GWAS sumstats file. \nUse code such as '%s' to fix",key_column,code_example))
    }
  }
  
  # Check there is at least one signed sumstats column
  print("Checking that there is at least one signed sumstats column (eg.: Z, OR, BETA, LOG_ODDS, SIGNED_SUMSTAT)")
  signed_stat_column_names = c("Z","OR","BETA","LOG_ODDS","SIGNED_SUMSTAT")
  if(sum(signed_stat_column_names %in% col_headers)<1 %in% col_headers){
    print("Header of file:")
    print(rows_of_data)
    stop("ERROR: cannot find a column name representing signed statistic in GWAS sumstats file. I.e. Z, OR, BETA")
  }
  
  # Check that first three column headers are SNP, CHR, BP (in that order)
  print("Checking that the first three column headers are SNP, CHR and BP in this order.")
  if(!sum(col_headers[1:3]==c("SNP","CHR","BP"))==3){
    whichSNP = which(col_headers=="SNP")[1]
    whichCHR = which(col_headers=="CHR")[1]
    whichBP = which(col_headers=="BP")[1]
    otherCols = setdiff(1:length(col_headers),c(whichSNP,whichCHR,whichBP))
    #x=read.table(path)
    x=fread(path)
    #write.table(x=x[,c(whichSNP,whichCHR,whichBP,otherCols)], file=path, sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
    xx = setcolorder(x, c(whichSNP,whichCHR,whichBP,otherCols))
    fwrite(x=xx, file=path, sep="\t")
    sumstats_file <- readLines(path)
  }
  
  # The formatting process can (rarely) result in duplicated columns, i.e. CHR, if CHR:BP is expanded and one already exists... delete duplicates
  print("Removing duplicated columns (if any).")
  rows_of_data <- c(sumstats_file[1], sumstats_file[2]); col_headers = strsplit(rows_of_data[1], "\t")[[1]]
  if(sum(duplicated(col_headers))>0){
    notDup = which(!duplicated(col_headers))
    write.table(x=read.table(path)[,notDup], file=path, sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE); sumstats_file <- readLines(path)
  }
  

  sumstats <- fread(path)
    # MAGMA cannot handle P-values as low as 3e-400... so convert them to zeros
    #if (as.logical(as.numeric(readline("MAGMA cannot handle P-values as low as 3e-400. Do you want MAGMA.celltyping to convert any (if) existing ones to zeroes? 0 for NO, 1 for YES")))) {
  rows_of_data <- c(sumstats_file[1], sumstats_file[2]); col_headers = strsplit(rows_of_data[1], "\t")[[1]]
  sumstats$P = as.numeric(as.character(sumstats$P)) # Note, I've not tested this since changing it from the original code... which was using gawk/sed
    #}
  fwrite(x=sumstats, file=path, sep="\t") #, quote=FALSE, row.names = FALSE, col.names = FALSE)
    
  sumstats_file <- readLines(path)
  
  # Sometimes the N column is not all integers... so round it up
  # - I edited this on 3rd April 2020 so it uses data.table... but don't have a dataset to check that it still works
  if("N" %in% col_headers) {
    #whichN = which(col_headers %in% "N")
    print("Sometimes the N column is not all integers. MAGMA.celltyping will round them up if such instances exist.")
      #rows_of_data <- c(sumstats_file[1], sumstats_file[2]); col_headers = strsplit(rows_of_data[1], "\t")[[1]]
      sumstats <- fread(path) # read.table(path)
      sumstats$N = round(as.numeric(as.character(sumstats$N))) # Note, I've not tested this since changing it from the original code... which was using gawk/sed
      fwrite(x=sumstats, file=path, sep="\t") #, quote=FALSE, row.names = FALSE, col.names = FALSE)
      #for (i in seq_along(sumstats[,which(col_headers=="N")])) {
      #  if (sumstats[i,which(col_headers=="N")]=="N") {next} # To skip the header.
      #  sumstats[i,which(col_headers=="N")] <- round(as.numeric(as.character(sumstats[i,which(col_headers=="N")]))) # This converts anything under 3e-400 to zeros.
      #}
      #write.table(x=sumstats, file=path, sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE); sumstats_file <- readLines(path)  
   } else {
     print("There is no N column in the summary statistics file. Adding N column based on N argument.")
     sumstats <- fread(path) # read.table(path)
     sumstats$N = N
     sumstats$N = round(as.numeric(as.character(sumstats$N))) # Note, I've not tested this since changing it from the original code... which was using gawk/sed
     fwrite(x=sumstats, file=path, sep="\t") #, quote=FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  # All rows should start with either SNP or rs... if they don't drop them
  print("Dropping all rows that don't start with 'rs'")
  sumstats_file <- readLines(path)
  sumstats_file <- c(sumstats_file[1],sumstats_file[grepl("^rs",sumstats_file)])
  writeLines(text=sumstats_file, con = path)
  
  # Keep only rows which have the number of columns expected (I've commmented this out because it takes forever to run... but presumably I wrote it for a reason!)
  # print("Keeping only rows which have the number of columns expected.")
  # sumstats_file <- readLines(path); expected_number_of_columns <- length(strsplit(sumstats_file[1],"\t")[[1]]); good_ones <- sumstats_file[1]
  # for (line in sumstats_file) {
  #   if (line == sumstats_file[1]) {next} # Skip header
  #   if ( length(strsplit(line,"\t")[[1]]) == expected_number_of_columns ) {good_ones <- c(good_ones, line)} # This adds every line that has expected number of columns to a temporary list
  # }
  # writeLines(text=good_ones, con = path)
  
  # Try to remove duplicated RSIDs
  print("Removing duplicated RSIDs.")
  #sumstats <- read.table(path)
  sumstats <- fread(path)
  if(sum(duplicated(sumstats[,1]))>0){
    notDup = which(!duplicated(sumstats[,1]))
    notDupLines=sumstats[notDup,]
    fwrite(notDupLines, file=path, sep="\t") #, quote=FALSE, row.names = FALSE, col.names = FALSE); sumstats_file <- readLines(path)
    rm(notDupLines);gc()
  }
  
  # Show how the data now looks
  print("Succesfully finished preparing sumstats file:")
  print("Header of file:")
  con <- file(path,"r") ; rows_of_data <- readLines(con,n=2) ; close(con)
  print(rows_of_data)
  col_headers = strsplit(rows_of_data[1],"\t")[[1]]
  
  return(tmp) # Returns address of modified file
}
