###################### this function is for the DNA methylation to inter

#' Extract training samples from the coverage file after Bismark call methylation
#'
#' @param cov_filename The file contain "chr", "start", "end", "methylation", "m.read" and "um.read"
#' @param flank_size Side chain length, total window width 2*flank_size
#' @param output_prefix prefix of output file
#' @param train_size Training sample size
#' @param ref_genome_dir Reference genome file directory, and individual chromosomes are stored separately in .fa file.
#' @param file_type cov_file type ".gz" or ".txt"
#'
#' @return A list containing Flank_methyl,Flank_DNA_onehot,Flank_DNA_2mer,Center_loci_table,Flank_m and Flank_um)
#' @export
#'
#' @examples
make_tain_Rdata <- function(cov_filename = "data/Methyl/ENCFF003FWN.0.5.deduplicated.bismark.cov.gz", flank_size = 100,
                                 output_prefix = "data/Train_Test_Data/ENCFF003FWN.0.5", train_size = 50000, ref_genome_dir,file_type = ".gz")
{
  source("R/DNA_coding.R")
  if (file_type == ".gz") {
    gf <- gzfile(cov_filename, "rt")
    cov_file <- read.delim(gf, header = FALSE, stringsAsFactors = FALSE)
    colnames(cov_file) <- c("chr", "start", "end", "methylation", "m.read", "um.read")
  } else {
    cov_file <- read.table(cov_filename, quote = "\"", comment.char = "", stringsAsFactors = FALSE)
    colnames(cov_file) <- c("chr", "start", "end", "methylation", "m.read", "um.read")
  }
  cov_file <- cov_file[cov_file[, "chr"] %in% c(paste0("chr", c(1:22)), "chrX", "chrY"), ]
  cov_file <-   cov_file[order(as.numeric(as.character(cov_file[,2])),decreasing=F),]
  cov_file <-   cov_file[order(as.character(cov_file[,1]),decreasing=F),]

  cov <- as.numeric(as.character(cov_file[, "m.read"])) + as.numeric(as.character(cov_file[, "um.read"]))

  train_methylation_start <- mean(cov) + 2 * sd(cov)
  train_methylation_end <- mean(cov) + 3 * sd(cov)

  ########## get the train data################
  train_index_all <- which(cov < train_methylation_end & cov > train_methylation_start)
  train_index_all_index <- sample(train_index_all, size = ceiling(train_size * 1.2),replace = F)  ####Sample more then the train size than del some bad sites
  train_index_all_index <- train_index_all_index[order(train_index_all_index, decreasing = F)]
  train_index_all_index <- train_index_all_index[which(train_index_all_index > (flank_size+1) & train_index_all_index < (nrow(cov_file)-flank_size-1))]

  Center_loci_table <- cov_file[train_index_all_index,]#cov_file is order

  Flank_methyl <- matrix(0, nrow = length(train_index_all_index), ncol = 2 * flank_size)
  Flank_m <- matrix(0, nrow = length(train_index_all_index), ncol = 2 * flank_size)
  Flank_um <- matrix(0, nrow = length(train_index_all_index), ncol = 2 * flank_size)

  Flank_DNA_onehot <- get_DNA_inter_mutiple_core(Center_loci_table[, c(1, 2)], flank_size, ref_genome_dir,OnehotOr2mer = 1)
  Flank_DNA_2mer <- get_DNA_inter_mutiple_core(Center_loci_table[, c(1, 2)], flank_size, ref_genome_dir,OnehotOr2mer = 2)
  for (i in 1:length(train_index_all_index)) {
    Flank_methyl[i, ] <- as.numeric(as.character(cov_file[c((train_index_all_index[i] - flank_size):
                                                              (train_index_all_index[i] - 1), (train_index_all_index[i] + 1):(train_index_all_index[i] + flank_size)), "methylation"]))
  }#to get the flank methylation, outmit the center site.
  for (i in 1:length(train_index_all_index)) {
    Flank_m[i, ] <- as.numeric(as.character(cov_file[c((train_index_all_index[i] - flank_size):
                                                         (train_index_all_index[i] - 1), (train_index_all_index[i] + 1):(train_index_all_index[i] + flank_size)), "m.read"]))
  }#to get the flank m reads count, outmit the center site.
  for (i in 1:length(train_index_all_index)) {
    Flank_um[i, ] <- as.numeric(as.character(cov_file[c((train_index_all_index[i] - flank_size):
                                                          (train_index_all_index[i] - 1), (train_index_all_index[i] + 1):(train_index_all_index[i] + flank_size)), "um.read"]))
  }#to get the flank um reads count, outmit the center site.
  # filt some bad loci#######
  out_filt <- c()
  for (i in 1:length(train_index_all_index)) {
    if (sum(apply(Flank_DNA_onehot[[i]], 1, sum)) < 0 || is.na(Flank_methyl[i, ])) {
      out_filt <- c(out_filt, i)
    }
  }
  if (length(out_filt) > 0) {
    Flank_DNA_onehot[out_filt] <- NULL
    Flank_DNA_onehot <- Flank_DNA_onehot[1:train_size]
    Flank_DNA_2mer[out_filt] <- NULL
    Flank_DNA_2mer <-  Flank_DNA_2mer[1:train_size]
    Center_loci_table <- Center_loci_table[-out_filt, ][1:train_size, ]
    Flank_m <- Flank_m[-out_filt, ][1:train_size, ]
    Flank_um <- Flank_um[-out_filt, ][1:train_size, ]
    Flank_methyl <- Flank_methyl[-out_filt, ][1:train_size, ]
  } else {
    Flank_DNA_onehot <- Flank_DNA_onehot[1:train_size]
    Flank_DNA_2mer <- Flank_DNA_2mer[1:train_size]
    Center_loci_table <- Center_loci_table[1:train_size, ]
    Flank_m <- Flank_m[1:train_size, ]
    Flank_um <- Flank_um[1:train_size, ]
    Flank_methyl <- Flank_methyl[1:train_size, ]
  }
  Train_data <- list(Flank_methyl,Flank_DNA_onehot,Flank_DNA_2mer,Center_loci_table,Flank_m,Flank_um)
  save(Train_data, file = paste0(output_prefix, ".train.Rdata"))
  rm(Flank_methyl,Flank_DNA_onehot,Flank_DNA_2mer,Center_loci_table,Flank_m,Flank_um)
}

#' Get the features of sites to be predicted from the coverage file after Bismark call methylation
#'
#' @param cov_filename The file contain "chr", "start", "end", "methylation", "m.read" and "um.read"
#' @param prediction_index Index of low coverage sites to be predicted，low coverage sites are contained in cov_file. Index satrt from 1.
#' @param flank_size Side chain length, total window width 2*flank_size
#' @param output_prefix prefix of output file
#' @param ref_genome_dir Reference genome file directory, and individual chromosomes are stored separately in .fa file.
#' @param file_type cov_file type ".gz" or ".txt"
#'
#' @return A list containing Flank_methyl,Flank_DNA_onehot,Flank_DNA_2mer,Center_loci_table,Flank_m and Flank_um)
#' @export
#'
#' @examples
make_test_Rdata <- function(cov_filename = "data/Methyl/ENCFF003FWN.0.5.deduplicated.bismark.cov.gz",prediction_index,  flank_size = 100,
                            output_prefix = "data/Train_Test_Data/ENCFF003FWN.0.5", ref_genome_dir ,file_type = ".gz") {
  source("R/DNA_coding.R")
  if (file_type == ".gz") {
    gf <- gzfile(cov_filename, "rt")
    cov_file <- read.delim(gf, header = FALSE, stringsAsFactors = FALSE)
    colnames(cov_file) <- c("chr", "start", "end", "methylation", "m.read", "um.read")
  } else {
    cov_file <- read.table(cov_filename, quote = "\"", comment.char = "", stringsAsFactors = FALSE)
    colnames(cov_file) <- c("chr", "start", "end", "methylation", "m.read", "um.read")
  }
  cov_file <- cov_file[cov_file[, "chr"] %in% c(paste0("chr", c(1:22)), "chrX", "chrY"), ]
  cov_file <-   cov_file[order(as.numeric(as.character(cov_file[,2])),decreasing=F),]
  cov_file <-   cov_file[order(as.character(cov_file[,1]),decreasing=F),]


  ########## get the test data################
  Center_loci_table <- cov_file[prediction_index,]#cov_file is order

  Flank_methyl <- matrix(0, nrow = nrow(Center_loci_table), ncol = 2 * flank_size)
  Flank_m <- matrix(0, nrow = nrow(Center_loci_table), ncol = 2 * flank_size)
  Flank_um <- matrix(0, nrow = nrow(Center_loci_table), ncol = 2 * flank_size)

  Flank_DNA_onehot <- get_DNA_inter_mutiple_core(Center_loci_table[, c(1, 2)], flank_size, ref_genome_dir,OnehotOr2mer = 1)
  Flank_DNA_2mer <- get_DNA_inter_mutiple_core(Center_loci_table[, c(1, 2)], flank_size, ref_genome_dir,OnehotOr2mer = 2)
  for (i in 1:nrow(Center_loci_table)) {
    Flank_methyl[i, ] <- as.numeric(as.character(cov_file[c((train_index_all_index[i] - flank_size):
                                                              (train_index_all_index[i] - 1), (train_index_all_index[i] + 1):(train_index_all_index[i] + flank_size)), "methylation"]))
  }#to get the flank methylation, outmit the center site.
  for (i in 1:nrow(Center_loci_table)) {
    Flank_m[i, ] <- as.numeric(as.character(cov_file[c((train_index_all_index[i] - flank_size):
                                                         (train_index_all_index[i] - 1), (train_index_all_index[i] + 1):(train_index_all_index[i] + flank_size)), "m.read"]))
  }#to get the flank m reads count, outmit the center site.
  for (i in 1:nrow(Center_loci_table)) {
    Flank_um[i, ] <- as.numeric(as.character(cov_file[c((train_index_all_index[i] - flank_size):
                                                          (train_index_all_index[i] - 1), (train_index_all_index[i] + 1):(train_index_all_index[i] + flank_size)), "um.read"]))
  }#to get the flank um reads count, outmit the center site.
  # filt some bad loci#######
  out_filt <- c()
  for (i in 1:nrow(Center_loci_table)) {
    if (sum(apply(Flank_DNA_onehot[[i]], 1, sum)) < 0 || is.na(Flank_methyl[i, ])) {
      out_filt <- c(out_filt, i)
    }
  }
  if (length(out_filt) > 0) {
    Flank_DNA_onehot[out_filt] <- NULL
    Flank_DNA_onehot <- Flank_DNA_onehot[1:train_size]
    Flank_DNA_2mer[out_filt] <- NULL
    Flank_DNA_2mer <-  Flank_DNA_2mer[1:train_size]
    Center_loci_table <- Center_loci_table[-out_filt, ][1:train_size, ]
    Flank_m <- Flank_m[-out_filt, ][1:train_size, ]
    Flank_um <- Flank_um[-out_filt, ][1:train_size, ]
    Flank_methyl <- Flank_methyl[-out_filt, ][1:train_size, ]
  } else {
    Flank_DNA_onehot <- Flank_DNA_onehot[1:train_size]
    Flank_DNA_2mer <- Flank_DNA_2mer[1:train_size]
    Center_loci_table <- Center_loci_table[1:train_size, ]
    Flank_m <- Flank_m[1:train_size, ]
    Flank_um <- Flank_um[1:train_size, ]
    Flank_methyl <- Flank_methyl[1:train_size, ]
  }
  Train_data <- list(Flank_methyl,Flank_DNA_onehot,Flank_DNA_2mer,Center_loci_table,Flank_m,Flank_um)
  save(Train_data, file = paste0(output_prefix, ".test.Rdata"))
  rm(Flank_methyl,Flank_DNA_onehot,Flank_DNA_2mer,Center_loci_table,Flank_m,Flank_um)
}

#' Extract specific length DNA sequence information of all CpG sites from cov file
#' @param flank_size Side chain length, total window width 2*flank_size
#' @param ref_genome_dir Reference genome file directory, and individual chromosomes are stored separately in .fa file.
#' @param CpG_report_file_name The file contains "chr" and "position" in ".txt" file
#' @param pre_DNA_inter_dir The folder of output file
#'
#'
#' @return Three ".Rdata" files in the pre_DNA_inter_dir.
#' @export
#'
#' @examples
get_all_onehot_2mer <- function(flank_size=50,ref_genome_dir="ref_hg38/hg38/",CpG_report_file_name="data/Methyl/H1_hESC.0.3.CpG_report.txt",pre_DNA_inter_dir="data/All_hg38_DNA/")
{
  cov_file <- read.delim(CpG_report_file_name, header = FALSE, stringsAsFactors = FALSE)
  cov_file <- cov_file[cov_file[,3]=="+",]
  cov_file <- cov_file[cov_file[, 1] %in% c(paste0("chr", c(1:22)), "chrX", "chrY"), ]
  cov_file <-   cov_file[order(as.numeric(as.character(cov_file[,2])),decreasing=F),]
  cov_file <-   cov_file[order(as.character(cov_file[,1]),decreasing=F),]
  for(temp_chr in c(paste0("chr", c(1:22)), "chrX", "chrY"))
  {
    cov_file_temp <- cov_file[cov_file[,1] == temp_chr,]
    Chr_Flank_DNA_2mer <- get_DNA_inter_mutiple_core(cov_file_temp [, c(1, 2)], flank_size, ref_genome_dir,OnehotOr2mer = 2)
    save(Chr_Flank_DNA_2mer,file = paste0(pre_DNA_inter_dir,temp_chr,"_Flank_DNA_2mer.Rdata"))
    save(cov_file_temp ,file = paste0(pre_DNA_inter_dir,temp_chr,"_CpG_plus_loci.Rdata"))
    Chr_Flank_DNA_onehot <- get_DNA_inter_mutiple_core(cov_file_temp[, c(1, 2)], flank_size, ref_genome_dir,OnehotOr2mer = 1)
    save(Chr_Flank_DNA_onehot,file = paste0(pre_DNA_inter_dir,temp_chr,"_Flank_DNA_onehot.Rdata"))
  }
}


#' Extract all sites' features with high coverage from cov file

#' @param cov_filename The file contain "chr", "start", "end", "methylation", "m.read" and "um.read"
#' @param flank_size Side chain length, total window width 2*flank_size
#' @param output_prefix prefix of output file
#' @param pre_DNA_inter_dir The folder of output file
#' @param file_type cov_file type ".gz" or ".txt"
#' @param cov_cutoff >cov_cutoff site were used as training data
#' @param cov_high_threshold There are extremely abnormal high coverage sites during sequencing.
#' A fixed value can be set according to prior knowledge to filter out these sites.
#' @return A list containing Flank_methyl,Flank_DNA_onehot,Flank_DNA_2mer,Center_loci_table,Flank_m and Flank_um)
#' @export
#'
#' @examples
make_lowcov_train_data_from_all <- function(cov_filename = "data/Methyl/GM12878.0.5.strand_plus.cov", flank_size = 50,
                                           output_prefix = "data/Train_Test_Data/GM12878.0.5.lowcov", pre_DNA_inter_dir="data/All_hg38_DNA/",
                                           file_type = ".gz",
                                           cov_cutoff=3,cov_high_threshold=1000)
{
  source("R/DNA_coding.R")
  library("foreach")
  library("doParallel")
  if (file_type == ".gz") {
    gf <- gzfile(cov_filename, "rt")
    cov_file <- read.delim(gf, header = FALSE, stringsAsFactors = FALSE)
    colnames(cov_file) <- c("chr", "start", "end", "methylation", "m.read", "um.read")
  } else {
    cov_file <- read.table(cov_filename, quote = "\"", comment.char = "", stringsAsFactors = FALSE)
    colnames(cov_file) <- c("chr", "start", "end", "methylation", "m.read", "um.read")
  }
  cov_file <- cov_file[cov_file[, "chr"] %in% c(paste0("chr", c(1:22)), "chrX", "chrY"), ]
  cov_file <-   cov_file[order(as.numeric(as.character(cov_file[,2])),decreasing=F),]
  ########## get the test data################

  for(temp_chr in c(paste0("chr", c(1:22)), "chrX", "chrY"))
  {
    print(temp_chr)
    temp_cov_file <- cov_file[cov_file[,"chr"]==temp_chr,]
    temp_cov <- as.numeric(as.character(temp_cov_file[, "m.read"])) + as.numeric(as.character(temp_cov_file[, "um.read"]))

    test_index_all_index <- which(temp_cov > cov_cutoff & temp_cov < cov_high_threshold)
    Center_loci_table <- temp_cov_file[test_index_all_index,]#cov_file may ordered in chr1 chr10...


    Flank_DNA_onehot <- list()
    Flank_DNA_2mer <- list()
    Center_loci_table_index <- paste0(Center_loci_table[,1],"_",Center_loci_table[,2])
    print("start make 2_mer")

    load(paste0(pre_DNA_inter_dir,temp_chr,"_Flank_DNA_2mer.Rdata"))
    load(paste0(pre_DNA_inter_dir,temp_chr,"_Flank_DNA_onehot.Rdata"))
    load(paste0(pre_DNA_inter_dir,temp_chr,"_CpG_plus_loci.Rdata"))
    Chr_CpG_plus_index <- paste0(cov_file_temp[,1],"_",cov_file_temp[,2])
    index_match <- match(Center_loci_table_index,Chr_CpG_plus_index)

    print("Convert the flanking DNA sequence into onehot and 2mer")
    Flank_DNA_onehot[which(!is.na(index_match))] <- Chr_Flank_DNA_onehot[index_match[which(!is.na(index_match))]]
    Flank_DNA_2mer[which(!is.na(index_match))] <- Chr_Flank_DNA_2mer[index_match[which(!is.na(index_match))]]

    print("start make parallel")
    cl <- makeCluster(8)
    registerDoParallel(cl)
    Flank_methyl <- foreach(i =c( 1:length(test_index_all_index)),.combine="rbind")  %dopar% {
      max_t <- nrow(temp_cov_file)
      tmep_methyl <- rep(0, 2 * flank_size)
      if(test_index_all_index[i] - flank_size < 1 || (test_index_all_index[i] + flank_size) > max_t)
      {
        tmep_methyl <- -1
      } else
      {
        tmep_methyl<- as.numeric(as.character(temp_cov_file[c((test_index_all_index[i] - flank_size):(test_index_all_index[i] - 1), (test_index_all_index[i] + 1):(test_index_all_index[i] + flank_size)), "methylation"]))
      }
      return(tmep_methyl)
    }
    stopCluster(cl)

    cl <- makeCluster(8)
    registerDoParallel(cl)
    Flank_m <- foreach(i =c( 1:length(test_index_all_index)),.combine="rbind")  %dopar% {
      max_t <- nrow(temp_cov_file)
      tmep_m <- rep(0, 2 * flank_size)
      if(test_index_all_index[i] - flank_size < 1 || (test_index_all_index[i] + flank_size) > max_t)
      {
        tmep_m  <- 0
      } else
      {
        tmep_m <- as.numeric(as.character(temp_cov_file[c((test_index_all_index[i] - flank_size): (test_index_all_index[i] - 1), (test_index_all_index[i] + 1):(test_index_all_index[i] + flank_size)), "m.read"]))
      }
      return(tmep_m)
    }
    stopCluster(cl)

    cl <- makeCluster(8)
    registerDoParallel(cl)
    Flank_um <- foreach(i =c( 1:length(test_index_all_index)),.combine="rbind")  %dopar% {
      max_t <- nrow(temp_cov_file)
      tmep_m <- rep(0, 2 * flank_size)
      if(test_index_all_index[i] - flank_size < 1 || (test_index_all_index[i] + flank_size) > max_t)
      {
        tmep_um <- 0
      } else
      {
        tmep_um <- as.numeric(as.character(temp_cov_file[c((test_index_all_index[i] - flank_size): (test_index_all_index[i] - 1), (test_index_all_index[i] + 1):(test_index_all_index[i] + flank_size)), "um.read"]))
      }
      return(tmep_um)
    }
    stopCluster(cl)

    print(paste("Stop parallel",length(test_index_all_index)))


    print("start filt some bad loci")
    out_filt <- c()
    for (i in 1:length(test_index_all_index)) {
      if(length(Flank_DNA_2mer[[i]])==0)
      {
        out_filt <- c(out_filt, i)
      }else{
        if (sum(apply(Flank_DNA_2mer[[i]], 1, sum)) < 0 || is.na(Flank_methyl[i, ]) || sum(Flank_methyl[i, ]) < 0 ) {
          out_filt <- c(out_filt, i)
        }
      }
    }

    if (length(out_filt) > 0) {
      keep <- c(1:length(test_index_all_index))[-out_filt]
      Flank_DNA_onehot <- Flank_DNA_onehot[keep]
      Flank_DNA_2mer <-  Flank_DNA_2mer[keep]
      Center_loci_table <- Center_loci_table[-out_filt, ]
      Flank_m <- Flank_m[-out_filt, ]
      Flank_um <- Flank_um[-out_filt, ]
      Flank_methyl <- Flank_methyl[-out_filt, ]
    }
    Test_data <- list(Flank_methyl,Flank_DNA_onehot,Flank_DNA_2mer,Center_loci_table,Flank_m,Flank_um)
    save(Test_data, file = paste0(output_prefix,"_" ,temp_chr,".train.Rdata"))
    rm(Flank_methyl,Flank_DNA_onehot,Flank_DNA_2mer,Center_loci_table,Flank_m,Flank_um)
  }
}

#' Extract all sites' features with low coverage from cov file

#' @param cov_filename The file contain "chr", "start", "end", "methylation", "m.read" and "um.read"
#' @param flank_size Side chain length, total window width 2*flank_size
#' @param output_prefix prefix of output file
#' @param pre_DNA_inter_dir The folder of output file
#' @param file_type cov_file type ".gz" or ".txt"
#' @param cov_high_threshold There are extremely abnormal high coverage sites during sequencing.
#' A fixed value can be set according to prior knowledge to filter out these sites.
#' @return A list containing Flank_methyl,Flank_DNA_onehot,Flank_DNA_2mer,Center_loci_table,Flank_m and Flank_um)
#' @export
#'
#' @examples
make_lowcov_test_data_from_all <- function(cov_filename = "data/Methyl/GM12878.0.5.strand_plus.cov", flank_size = 50,
                                             output_prefix = "data/Train_Test_Data/GM12878.0.5.lowcov", pre_DNA_inter_dir="data/All_hg38_DNA/",
                                             file_type = ".gz",
                                             cov_high_threshold=1000)
{
  source("R/DNA_coding.R")
  library("foreach")
  library("doParallel")
  if (file_type == ".gz") {
    gf <- gzfile(cov_filename, "rt")
    cov_file <- read.delim(gf, header = FALSE, stringsAsFactors = FALSE)
    colnames(cov_file) <- c("chr", "start", "end", "methylation", "m.read", "um.read")
  } else {
    cov_file <- read.table(cov_filename, quote = "\"", comment.char = "", stringsAsFactors = FALSE)
    colnames(cov_file) <- c("chr", "start", "end", "methylation", "m.read", "um.read")
  }
  cov_file <- cov_file[cov_file[, "chr"] %in% c(paste0("chr", c(1:22)), "chrX", "chrY"), ]
  cov_file <-   cov_file[order(as.numeric(as.character(cov_file[,2])),decreasing=F),]
  ########## get the test data################

  for(temp_chr in c(paste0("chr", c(1:22)), "chrX", "chrY"))
  {
    print(temp_chr)
    temp_cov_file <- cov_file[cov_file[,"chr"]==temp_chr,]
    temp_cov <- as.numeric(as.character(temp_cov_file[, "m.read"])) + as.numeric(as.character(temp_cov_file[, "um.read"]))

    test_index_all_index <- which(temp_cov <=3 | temp_cov > cov_high_threshold)
    Center_loci_table <- temp_cov_file[test_index_all_index,]#cov_file may ordered in chr1 chr10...


    Flank_DNA_onehot <- list()
    Flank_DNA_2mer <- list()
    Center_loci_table_index <- paste0(Center_loci_table[,1],"_",Center_loci_table[,2])
    print("start make 2_mer")

    load(paste0(pre_DNA_inter_dir,temp_chr,"_Flank_DNA_2mer.Rdata"))
    load(paste0(pre_DNA_inter_dir,temp_chr,"_Flank_DNA_onehot.Rdata"))
    load(paste0(pre_DNA_inter_dir,temp_chr,"_CpG_plus_loci.Rdata"))
    Chr_CpG_plus_index <- paste0(cov_file_temp[,1],"_",cov_file_temp[,2])
    index_match <- match(Center_loci_table_index,Chr_CpG_plus_index)

    print("Convert the flanking DNA sequence into onehot and 2mer")
    Flank_DNA_onehot[which(!is.na(index_match))] <- Chr_Flank_DNA_onehot[index_match[which(!is.na(index_match))]]
    Flank_DNA_2mer[which(!is.na(index_match))] <- Chr_Flank_DNA_2mer[index_match[which(!is.na(index_match))]]

    print("start make parallel")
    cl <- makeCluster(8)
    registerDoParallel(cl)
    Flank_methyl <- foreach(i =c( 1:length(test_index_all_index)),.combine="rbind")  %dopar% {
      max_t <- nrow(temp_cov_file)
      tmep_methyl <- rep(0, 2 * flank_size)
      if(test_index_all_index[i] - flank_size < 1 || (test_index_all_index[i] + flank_size) > max_t)
      {
        tmep_methyl <- -1
      } else
      {
        tmep_methyl<- as.numeric(as.character(temp_cov_file[c((test_index_all_index[i] - flank_size):(test_index_all_index[i] - 1), (test_index_all_index[i] + 1):(test_index_all_index[i] + flank_size)), "methylation"]))
      }
      return(tmep_methyl)
    }
    stopCluster(cl)

    cl <- makeCluster(8)
    registerDoParallel(cl)
    Flank_m <- foreach(i =c( 1:length(test_index_all_index)),.combine="rbind")  %dopar% {
      max_t <- nrow(temp_cov_file)
      tmep_m <- rep(0, 2 * flank_size)
       if(test_index_all_index[i] - flank_size < 1 || (test_index_all_index[i] + flank_size) > max_t)
      {
        tmep_m  <- 0
      } else
      {
        tmep_m <- as.numeric(as.character(temp_cov_file[c((test_index_all_index[i] - flank_size): (test_index_all_index[i] - 1), (test_index_all_index[i] + 1):(test_index_all_index[i] + flank_size)), "m.read"]))
      }
      return(tmep_m)
    }
    stopCluster(cl)

    cl <- makeCluster(8)
    registerDoParallel(cl)
    Flank_um <- foreach(i =c( 1:length(test_index_all_index)),.combine="rbind")  %dopar% {
      max_t <- nrow(temp_cov_file)
      tmep_m <- rep(0, 2 * flank_size)
      if(test_index_all_index[i] - flank_size < 1 || (test_index_all_index[i] + flank_size) > max_t)
      {
        tmep_um <- 0
      } else
      {
        tmep_um <- as.numeric(as.character(temp_cov_file[c((test_index_all_index[i] - flank_size): (test_index_all_index[i] - 1), (test_index_all_index[i] + 1):(test_index_all_index[i] + flank_size)), "um.read"]))
      }
      return(tmep_um)
    }
    stopCluster(cl)

    print(paste("Stop parallel",length(test_index_all_index)))


    print("start filt some bad loci")
    out_filt <- c()
    for (i in 1:length(test_index_all_index)) {
      if(length(Flank_DNA_2mer[[i]])==0)
      {
        out_filt <- c(out_filt, i)
      }else{
      if (sum(apply(Flank_DNA_2mer[[i]], 1, sum)) < 0 || is.na(Flank_methyl[i, ]) || sum(Flank_methyl[i, ]) < 0 ) {
        out_filt <- c(out_filt, i)
      }
      }
    }

    if (length(out_filt) > 0) {
      keep <- c(1:length(test_index_all_index))[-out_filt]
      Flank_DNA_onehot <- Flank_DNA_onehot[keep]
      Flank_DNA_2mer <-  Flank_DNA_2mer[keep]
      Center_loci_table <- Center_loci_table[-out_filt, ]
      Flank_m <- Flank_m[-out_filt, ]
      Flank_um <- Flank_um[-out_filt, ]
      Flank_methyl <- Flank_methyl[-out_filt, ]
    }
    Test_data <- list(Flank_methyl,Flank_DNA_onehot,Flank_DNA_2mer,Center_loci_table,Flank_m,Flank_um)
    save(Test_data, file = paste0(output_prefix,"_" ,temp_chr,".test.Rdata"))
    rm(Flank_methyl,Flank_DNA_onehot,Flank_DNA_2mer,Center_loci_table,Flank_m,Flank_um)
  }
}
