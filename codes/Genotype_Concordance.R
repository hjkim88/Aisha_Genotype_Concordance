###
#   File name : Genotype_Concordance.R
#   Author    : Hyunjin Kim
#   Date      : Apr 9, 2020
#   Email     : hyunjin.kim@stjude.org
#   Purpose   : There are two SNP list and we want to get/check the overlapped SNP names
#               and compute genotype concordance of them. If they were vcf files, then
#               I would just use PICARD tool for it but since they are all csv/excel files
#               I will implement it in R. One thing to be aware is, one file (csv) has
#               a large file size (54GB).
#
#   Instruction
#               1. Source("Genotype_Concordance.R")
#               2. Run the function "compute_concordance" - specify the input paths and the output directory
#               3. The results will be generated under the output directory
#
#   Example
#               > source("The_directory_of_Genotype_Concordance.R/Genotype_Concordance.R")
#               > compute_concordance(list1_path="C:/Users/hkim8/SJ/Aisha_SNP/Fluidigm Genotypes.xlsx",
#                                     list2_path="C:/Users/hkim8/SJ/Aisha_SNP/GSA Results with Metadata.csv",
#                                     overlap_path="C:/Users/hkim8/SJ/Aisha_SNP/Overlap SNPs.xlsx",
#                                     chunk_size=10000,
#                                     gc_score_threshold=0.15,
#                                     outputDir="./results/")
#
#   * chunk_size = the number of lines that each chunk will contain
#                  only those lines will be on the memory for each iteration
#   * gc_score_threshold = we will only use SNPs that have GC score > [gc_score_threshold] for the comparison
#                          if it is set to 0, that means there will be no filtering based on the GC score
#
###

compute_concordance <- function(list1_path="C:/Users/hkim8/SJ/Aisha_SNP/Fluidigm Genotypes.xlsx",
                                list2_path="C:/Users/hkim8/SJ/Aisha_SNP/GSA Results with Metadata.csv",
                                overlap_path="C:/Users/hkim8/SJ/Aisha_SNP/Overlap SNPs.xlsx",
                                chunk_size=10000,
                                gc_score_threshold=0.15,
                                outputDir="./results/") {
  
  ### load libraries
  if(!require(xlsx, quietly = TRUE)) {
    install.packages("xlsx")
    require(xlsx, quietly = TRUE)
  }
  if(!require(LaF, quietly = TRUE)) {
    install.packages("LaF")
    require(LaF, quietly = TRUE)
  }
  if(!require(gplots, quietly = TRUE)) {
    install.packages("gplots")
    require(gplots, quietly = TRUE)
  }
  
  ### load the data
  list1 <- read.xlsx2(list1_path, sheetIndex = 1,
                      stringsAsFactors = FALSE, check.names = FALSE)
  overlap <- read.xlsx2(overlap_path, sheetIndex = 1,
                        stringsAsFactors = FALSE, check.names = FALSE)[,"RsID"]
  
  ### get list1 sample IDs
  list1_sample_ids <- unique(list1$ID)
  list1_sample_ids <- list1_sample_ids[order(list1_sample_ids)]
  
  ### get the numer of rows in list2 file
  ### if you know the exact number already, just type it
  ### but if you do not know, use the code below:
  ### list2_nrow <- determine_nlines(list2_path)-1
  ### this will take 5 mins, so if you just type the number,
  ### you will save 5 mins
  list2_nrow <- 205162026
  
  ### set the number of iteration
  iteration_num <- as.integer(list2_nrow/chunk_size)
  
  ### define columns we need
  ### colnames(read.csv(list2_path, nrows = 1, row.names = 1)) will let you know the column names of the file
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
  model <- detect_dm_csv(list2_path, sep=",", header=TRUE,
                         colClasses = "character",
                         stringsAsFactors = FALSE, check.names = FALSE)
  df.laf <- laf_open(model)
  
  ### process residuals first then iterate the others
  data <- next_block(df.laf, nrows=(list2_nrow-(iteration_num*chunk_size)))
  idx <- intersect(which(data[,"RsID"] %in% overlap),
                   intersect(which(data[,"Sample.Name.x"] %in% list1_sample_ids),
                             which(data[,"Sample.Name.x"] == data[,"Sample.Name.y"])))
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
    idx <- intersect(which(data[,"RsID"] %in% overlap),
                     intersect(which(data[,"Sample.Name.x"] %in% list1_sample_ids),
                               which(data[,"Sample.Name.x"] == data[,"Sample.Name.y"])))
    
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
  
  ### remove duplicates (the list2 originally has the duplicates)
  filtered_data <- filtered_data[which(!duplicated(filtered_data)),]
  
  ### save the filtered data
  write.xlsx2(filtered_data, file = paste0(outputDir, "filtered_Illumina_data.xlsx"), row.names = FALSE)
  
  
  ### save the filtered data to the list2
  list2 <- filtered_data
  rm(filtered_data)
  gc()
  
  ### numerize some character columns in the list2
  list2[,c("Chr",
           "Position",
           "GC.Score")] <- sapply(list2[c("Chr",
                                          "Position",
                                          "GC.Score")], as.numeric)
  
  ### filter the list2 with the GC score threshold
  list2 <- list2[which(as.numeric(list2$GC.Score) > gc_score_threshold),]
  
  ### create a matrix (SNP x Sample)
  common_snps <- intersect(unique(list1$SNP), unique(list2$RsID))
  common_samples <- intersect(unique(list1$ID), unique(list2$Sample.Name.x))
  concordance_snp_sample <- matrix(NA, length(common_snps), length(common_samples))
  rownames(concordance_snp_sample) <- common_snps
  colnames(concordance_snp_sample) <- as.character(common_samples)
  
  ### compute the concordance of all the SNPs and all the samples
  for(snp in common_snps) {
    for(samp in common_samples) {
      ### get specific index per list
      list1_idx <- intersect(which(list1$SNP == snp), which(list1$ID == samp))  
      list2_idx <- intersect(which(list2$RsID == snp), which(list2$Sample.Name.x == samp))
      
      ### compute the concordance (0: no intersection, 1: half-same, 2: exactly the same)
      if((length(list1_idx) > 0) && (length(list2_idx) > 0)) {
        list1_genotypes <- c(substr(list1$Genotype[list1_idx], 1, 1), substr(list1$Genotype[list1_idx], 2, 2))
        list2_genotypes <- c(list2$Allele1...Forward[list2_idx], list2$Allele2...Forward[list2_idx])
        cnt <- 0
        for(i in 1:length(list1_genotypes)) {
          for(j in 1:length(list2_genotypes)) {
            if(list1_genotypes[i] == list2_genotypes[j]) {
              cnt <- cnt + 1
              break
            }
          }
        }
        concordance_snp_sample[snp, samp] <- cnt
      }
    }
  }
  
  ### write out the concordance result
  write.xlsx2(data.frame(RsID=rownames(concordance_snp_sample), concordance_snp_sample,
                         stringsAsFactors = FALSE, check.names = FALSE),
              file = paste0(outputDir, "table_genotype_concordance.xlsx"),
              row.names = FALSE)
  
  ### create a heatmap with the concordance
  png(paste0(outputDir, "heatmap_genotype_concordance.png"),
      width = 1200, height = 1000, res = 150)
  heatmap.2(concordance_snp_sample,
            main = "Genotype Concordance",
            xlab = "", ylab = "",
            col=c("orange", "gray", "green2", "darkgreen"),
            scale="none", key=FALSE,
            dendrogram = 'none', trace = 'none',
            labRow = rownames(concordance_snp_sample),
            labCol = FALSE,
            Rowv = FALSE, Colv = FALSE,
            adjRow = c(0, 0.5),
            adjCol = c(1, 0.5),
            cexRow = 1, cexCol = 0.8,
            na.rm = TRUE, na.color = "gray",
            margins = c(3, 6))
  legend("left", title = "Concordance",
         legend = c("NA", "0", "1", "2"),
         fill = c("gray", "orange", "green2", "darkgreen"), cex = 1, box.lty = 0)
  dev.off()
  
  ### create an empty data frame for the concordance (per snp and per sample)
  concordance_snp <- data.frame(RsID=common_snps,
                                Percentage=rep("", length(common_snps)),
                                Numbers=rep("", length(common_snps)),
                                stringsAsFactors = FALSE, check.names = FALSE)
  rownames(concordance_snp) <- common_snps
  concordance_sample <- data.frame(SampleID=common_samples,
                                   Percentage=rep("", length(common_samples)),
                                   Numbers=rep("", length(common_samples)),
                                   stringsAsFactors = FALSE, check.names = FALSE)
  rownames(concordance_sample) <- common_samples
  
  ### find the concordance per SNP
  for(snp in common_snps) {
    same_num <- length(which(concordance_snp_sample[snp,] == 2))
    total_num <- length(which(!is.na(concordance_snp_sample[snp,])))
    if(same_num == 0) {
      concordance_snp[snp,"Percentage"] <- 0
    } else {
      concordance_snp[snp,"Percentage"] <- signif((same_num/total_num), digits = 3)
    }
    concordance_snp[snp,"Numbers"] <- paste(same_num, "/", total_num)
  }
  
  ### find the concordance per sample
  for(samp in common_samples) {
    same_num <- length(which(concordance_snp_sample[,samp] == 2))
    total_num <- length(which(!is.na(concordance_snp_sample[,samp])))
    if(same_num == 0) {
      concordance_sample[samp,"Percentage"] <- 0
    } else {
      concordance_sample[samp,"Percentage"] <- signif((same_num/total_num), digits = 3)
    }
    concordance_sample[samp,"Numbers"] <- paste(same_num, "/", total_num)
  }
  
  ### save the results
  write.xlsx2(concordance_snp,
              file = paste0(outputDir, "table_genotype_concordance_per_snp.xlsx"),
              row.names = FALSE)
  write.xlsx2(concordance_sample,
              file = paste0(outputDir, "table_genotype_concordance_per_sample.xlsx"),
              row.names = FALSE)
  
}
