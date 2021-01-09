###
#   File name : Filtering_AIMs.R
#   Author    : Hyunjin Kim
#   Date      : Jan 8, 2021
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : Load the large Ancestry Informative Marker chip result file
#               and filter it to only contain genotypes for SNP rs12768894.
#
#   Instruction
#               1. Source("Filtering_AIMs.R")
#               2. Run the function "filter_AIMs" - specify the input paths and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Filtering_AIMs.R/Filtering_AIMs.R")
#               > filter_AIMs(Aims_result_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Aisha/GSA/GSA Results with Metadata.csv",
#                             chunk_size=10000,
#                             outputDir="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin_Results/Aisha/Filtered_AIMs/")
#
#   * chunk_size = the number of lines that each chunk will contain
#                  only those lines will be on the memory for each iteration
#
###

filter_AIMs <- function(Aims_result_path="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Aisha/GSA/GSA Results with Metadata.csv",
                        chunk_size=10000,
                        outputDir="Z:/ResearchHome/ResearchHomeDirs/thomagrp/common/Hyunjin_Results/Aisha/Filtered_AIMs/") {
  
  ### load libraries
  if(!require(LaF, quietly = TRUE)) {
    install.packages("LaF")
    require(LaF, quietly = TRUE)
  }
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  
  ### get the numer of rows in list2 file
  ### if you know the exact number already, just type it
  ### but if you do not know, use the code below:
  ### Aims_result_nrow <- determine_nlines(Aims_result_path)-1
  ### this will take 5 mins, so if you just type the number,
  ### you will save 5 mins
  Aims_result_nrow <- 205162026
  
  ### set the number of iteration
  iteration_num <- as.integer(Aims_result_nrow/chunk_size)
  
  ### define columns we need
  ### colnames(read.csv(Aims_result_path, nrows = 1, row.names = 1)) will let you know the column names of the file
  columns <- c("SNP.Name",
               "Chr",
               "Position",
               "RsID",
               "GC.Score",
               "Sample.Name.x",
               "Sample.Name.y",
               "Allele1...Forward",
               "Allele2...Forward",
               "SNP")
  
  ### set progress bar
  pb <- txtProgressBar(min = 0, max = iteration_num, style = 3)
  
  ### detect a data model for the file and open a connection
  model <- detect_dm_csv(Aims_result_path, sep=",", header=TRUE,
                         colClasses = "character",
                         stringsAsFactors = FALSE, check.names = FALSE)
  df.laf <- laf_open(model)
  
  ### process residuals first then iterate the others
  data <- next_block(df.laf, nrows=(Aims_result_nrow-(iteration_num*chunk_size)))
  idx <- which(data[,"RsID"] == "rs12768894")
  filtered_data <- data[idx,columns]
  
  ### start time
  start_time <- Sys.time()
  
  ### iteratively load the data and filter it
  ### if it was to handle the whole data, it would be faster to compute
  ### concordance while streaming, but in this case, what we are comparing
  ### is a smaller subset than the 54GB, so we will get the data we need
  ### and compute the concordance later with the downsized subset.
  ### getting the subset would take 40 mins but once we have it,
  ### it would be quick to compute the concordance.
  ### if we compute it in streaming, it would be slower in this case.
  for(i in 1:iteration_num) {
    ### load the next block
    data <- next_block(df.laf, nrows=chunk_size)
    
    ### get indicies of rows we want to keep
    idx <- which(data[,"RsID"] == "rs12768894")
    
    ### update the filtered data
    filtered_data <- rbind(filtered_data, data[idx,columns])
    
    ### update the progress bar
    setTxtProgressBar(pb, i)
    
    ### garbage collection
    gc()
  }
  
  ### close the connections
  close(pb)
  close(df.laf)
  
  ### end time
  end_time <- Sys.time()
  
  ### print out the running time
  cat(paste("Running Time:",
            signif(as.numeric(difftime(end_time, start_time, units = "mins")), digits = 3),
            "mins"))
  
  ### remove duplicates (the Aims_result originally has the duplicates)
  filtered_data <- filtered_data[which(!duplicated(filtered_data)),]
  
  ### create outputDir
  dir.create(outputDir, showWarnings = TRUE, recursive = TRUE)
  
  ### save the filtered data
  saveRDS(filtered_data, file = paste0(outputDir, "rs12768894_filtered_AIMs_Result.RDS"))
  write.xlsx2(filtered_data, file = paste0(outputDir, "rs12768894_filtered_AIMs_Result.xlsx"), row.names = FALSE)
  
}
