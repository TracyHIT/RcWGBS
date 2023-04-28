#' Here we  can use hg38 pre-extracted 2-mer
# cov_file :("chr", "start", "end", "methylation", "m.read", "um.read")
source("R/DNA_coding.R")
source("R/Methylation_DNA_coding.R")
make_lowcov_train_data_from_all(cov_filename = "data/Compare_sc/HepG2_1_hg38.cov", flank_size = 50,
                                            output_prefix = "data/Train_Test_Data/HepG2_1_hg38.lowcov", pre_DNA_inter_dir="data/All_DNA/",
                                            file_type = ".gz",
                                            cov_cutoff=3,cov_high_threshold=1000)

