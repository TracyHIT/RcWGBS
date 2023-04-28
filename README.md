# RcWGBS
Recall DNA methylation levels at low coverage sites using CNN model in WGBS. 

Here, 50 base pairs (bp) of DNA sequences from upstream and downstream regions, along with the methylation levels of 50 CpG sites, were selected as features for the model.

Refer to R/example.R to make train dataset
>make_tain_Rdata()

We preprocessed the 50bp sequences on both sides of all CpG sites in hg38. Due to the large amount of data, it is not possible to upload to github. Aslo can run >get_all_onehot_2mer() to prodcuce these dateset. The download link is as follows:

Refer to R/example.R to get features for the site need to predict
>make_lowcov_test_data_from_all()
