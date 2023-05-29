#' Here we  can use hg38 pre-extracted 2-mer
# cov_file :("chr", "start", "end", "methylation", "m.read", "um.read")
source("R/DNA_coding.R")
source("R/Methylation_DNA_coding.R")
source("R/Method_model_lowcov.R")
#get training data
make_tain_Rdata (cov_filename = "data/Compare_sc/HepG2_1_hg38.cov", flank_size = 50,
                            output_prefix = "data/Train_Test_Data/HepG2_1_hg38", train_size = 100000,ref_genome_dir="/home/zqluoximei/ref_hg38/hg38/",file_type = ".txt",cov_cutoff_low=50,cov_cutoff_high=-1)
#get the sites need to predict
make_lowcov_test_data_from_all(cov_filename = "data/Compare_sc/HepG2_1_hg38.cov", flank_size = 50,
                                            output_prefix = "data/Train_Test_Data/HepG2_1_hg38.lowcov", pre_DNA_inter_dir="data/All_DNA/",
                                            file_type = ".gz",cov_cutoff=3,
                                            cov_high_threshold=1000)

#to predict Bulk WGBS data
model_2mer_methyl_train(train_data = "data/Train_Test_Data/HepG2_1_hg38.train.Rdata",
                                   figure_history_prefix = "pre_model2_2mer_methyl_GM12878",flank_size=50,
                                   save_model = "model_2mer_methyl_GM12878.0.3.h5",epochs_Num=20,
                                   batch_size_Num = 512,
                                   validation_split_Num = 0.2,
                                   poolingsize = 2)
#to predict single-cell WGBS data
model_2mer_methyl_train_single_cell(train_data = "data/Train_Test_Data/HepG2_1_hg38.train.Rdata",
                                    figure_history_prefix = "pre_model2_2mer_methyl_GM12878_SC",flank_size=50,
                                    save_model = "model_2mer_methyl_GM12878.0.3.h5",epochs_Num=20,
                                    batch_size_Num = 512,
                                    validation_split_Num = 0.2,
                                    poolingsize = 2)
