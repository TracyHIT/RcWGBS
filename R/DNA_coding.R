###################### this function is for the DNA to inter
#' Using onehot coding to transform DNA sequence into a data matrix containing only 0，1
#'
#' The input argument is a short DNA sequence
#' Output is a data matrix
#'
#' @param seq a short DNA sequence
DNA_To_inter_onehot <- function(seq) {
    dir <- matrix(c(c(1, 0, 0, 0), c(0, 1, 0, 0), c(0, 0, 1, 0), c(0, 0, 0, 1)), nrow = 4)
    colnames(dir) <- c("A", "C", "G", "T")
    tensorflow_DNA <- matrix(0, nrow = 4, ncol = nchar(seq))
    seq <- toupper(seq)
    for (i in 1:nchar(seq)) {
        fit <- try(tensorflow_DNA[, i] <- dir[, substr(seq, i, i)], silent = TRUE)
        if ("try-error" %in% class(fit)) {
            tensorflow_DNA[, i] <- c(-1, -1, -1, -1)
        }
    }
    return(tensorflow_DNA)
}

###################### this function is for the Dinucleotide binary encoding to inter
#'  Encoding 2mer to transform a short DNA sequence into a data matrix containing only 0，1
#'
#' The input argument is a short DNA sequence
#' Output is a data matrix
#'  @param seq a short DNA sequence
DNA_To_inter_2mer <- function(seq) {
    dir <- c(0, 1, 2, 3)
    names(dir) <- c("A", "C", "G", "T")
    tensorflow_DNA <- matrix(0, nrow = 4, ncol = nchar(seq) - 1)
    seq <- toupper(seq)
    for (i in 1:(nchar(seq) - 1)) {
        temp_2mer <- substr(seq, i, i + 1)
        fit <- try(temp_Num <- dir[substr(temp_2mer, 1, 1)] * 4 + dir[substr(temp_2mer, 2, 2)], silent = TRUE)
        if ("try-error" %in% class(fit) | is.na(fit)) {
            tensorflow_DNA[, i] <- c(-1, -1, -1, -1)
        } else {
            for (j in 1:4) {
                tensorflow_DNA[j, i] <- floor(temp_Num/(2^(4 - j)))
                temp_Num <- temp_Num%%(2^(4 - j))
            }
            #print(tensorflow_DNA[, i])
        }
    }
    return(tensorflow_DNA)
}

###################### this function is getting the center_chr_loci on the genome to binary encoding

#' Using the genomic location to get a short sequence then to transform a data matrix containing only 0，1 with single core computing
#'
#' @param center_chr_loci dataframe or matrix that containing genomic changes and locations
#' @param flank_size Side chain length, total window width 2*flank_size
#' @param fre_gen_dir Reference genome file directory, and individual chromosomes are stored separately in .fa file.
#' @param OnehotOr2mer Index of 1 or 2，1 for onehot, 2 for 2mer
#'
#' @return A matrix contain 0,1
#' @export
#'
#' @examples
#' center_chr_loci<-data.frame(chr=c("chr1","chr2"),Position=c(100000,200000))
#' get_DNA_inter(center_chr_loci, 100, "data/hg38/", 1)
get_DNA_inter <- function(center_chr_loci, flank_size = 100, fre_gen_dir = "data/hg38/", OnehotOr2mer = 1) {
    library("Biostrings")
    library("rlist")
    chr <- ""
    Seq <- rep("", nrow(center_chr_loci))
    DNAtens <- list()
    if (OnehotOr2mer == 1) {
        for (j in 1:nrow(center_chr_loci)) {
            if (center_chr_loci[j, 1] != chr) {
                chr <- as.character(center_chr_loci[j, 1])
                s <- readDNAStringSet(paste0(fre_gen_dir, chr, ".fa"))
            }
            Seq[j] <- toString(subseq(s, start = as.numeric(as.character(center_chr_loci[j, 2])) - flank_size, end = as.numeric(as.character(center_chr_loci[j,
                2])) + flank_size))
            DNAtens <- list.append(DNAtens, DNA_To_inter_onehot (Seq[j]))
            # cat ('Has get ',j, ' sequence \n');
            if (floor(j/100) == j/100) {
                cat("Has get ", j, " sequence \n")
            }
        }
    } else {
        if (OnehotOr2mer == 2)
            for (j in 1:nrow(center_chr_loci)) {
                if (center_chr_loci[j, 1] != chr) {
                  chr <- as.character(center_chr_loci[j, 1])
                  s <- readDNAStringSet(paste0(fre_gen_dir, chr, ".fa"))
                }
                Seq[j] <- toString(subseq(s, start = as.numeric(as.character(center_chr_loci[j, 2])) - flank_size, end = as.numeric(as.character(center_chr_loci[j,
                  2])) + flank_size))
                DNAtens <- list.append(DNAtens, DNA_To_inter_2mer(Seq[j]))
                # cat ('Has get ',j, ' sequence \n');
                if (floor(j/100) == j/100) {
                  cat("Has get ", j, " sequence \n")
                }
            }
    }
    return(DNAtens)
}

#' Using the genomic location to get a short sequence then to transform a data matrix containing only 0，1 with 8 concurrent threads parallel computing
#'
#' @param center_chr_loci dataframe or matrix that containing genomic changes and locations
#' @param flank_size Side chain length, total window width 2*flank_size
#' @param fre_gen_dir Reference genome file directory, and individual chromosomes are stored separately in .fa file.
#' @param OnehotOr2mer Index of 1 or 2，1 for onehot, 2 for 2mer
#'
#' @return A matrix contain 0,1
#' @export
#'
#' @examples
#' center_chr_loci<-data.frame(chr=c("chr1","chr2"),Position=c(100000,200000))
#' get_DNA_inter(center_chr_loci, 100, "data/hg38/", 1)
get_DNA_inter_mutiple_core <- function(center_chr_loci, flank_size = 100, dir = "data/hg38/", OnehotOr2mer = 1) {

  library("rlist")
  library("foreach")
  library("doParallel")
  s <- NULL
  chr <- ""
  Seq <- rep("", nrow(center_chr_loci))
  if (OnehotOr2mer == 1) {
    cl <- makeCluster(8)
    registerDoParallel(cl)
    DNAtens <- foreach(j = c(1:nrow(center_chr_loci))) %dopar% {
      library("Biostrings")
      source("R/DNA_coding.R")
      if (center_chr_loci[j, 1] != chr) {
        chr <- as.character(center_chr_loci[j, 1])
        s <- readDNAStringSet(paste0(dir, chr, ".fa"))
      }
      Seq[j] <- toString(subseq(s, start = as.numeric(as.character(center_chr_loci[j, 2])) - flank_size, end = as.numeric(as.character(center_chr_loci[j, 2])) + flank_size))
      return(DNA_To_inter_onehot (Seq[j]))
    }
    stopCluster(cl)
  }
 else {
    if (OnehotOr2mer == 2)
       cl <- makeCluster(8)
       registerDoParallel(cl)
       DNAtens <- foreach(j = c(1:nrow(center_chr_loci))) %dopar% {
        library("Biostrings")
        source("R/DNA_coding.R")
        if (center_chr_loci[j, 1] != chr) {
          chr <- as.character(center_chr_loci[j, 1])
          s <- readDNAStringSet(paste0(dir, chr, ".fa"))
        }
        Seq[j] <- toString(subseq(s, start = as.numeric(as.character(center_chr_loci[j, 2])) - flank_size, end = as.numeric(as.character(center_chr_loci[j,                                                                                                                                                   2])) + flank_size))
        return( DNA_To_inter_2mer(Seq[j]))
      }
      stopCluster(cl)
  }
  return(DNAtens)
}

