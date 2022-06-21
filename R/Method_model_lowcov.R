
#' Training the CNN model
#'
#' @param train_data The train dataset produced by "make_tain_Rdata()". A list containing Flank_methyl,Flank_DNA_onehot,Flank_DNA_2mer,Center_loci_table,Flank_m and Flank_um)
#' @param figure_history_prefix The output of the history during training the model
#' @param flank_size Side chain length, total window width 2*flank_size 
#' @param save_model Save the taining model
#' @param epochs_Num Set the CNN model's epochs
#' @param batch_size_Num Set the CNN model's batch_size
#' @param validation_split_Num Set the CNN model's validation_split
#' @param poolingsize Set the poolingsize. Recommended setting is 2
#'
#' @return
#' @export
#'
#' @examples
model_2mer_methyl_train<- function(train_data = "data/Train_Test_Data/GM12878.0.3.train.Rdata",
                                  figure_history_prefix = "pre_model2_2mer_methyl_GM12878",flank_size=50,
                                  save_model = "model_2mer_methyl_GM12878.0.3.h5",epochs_Num=20,
                                  batch_size_Num = 512,
                                  validation_split_Num = 0.2,
                                  poolingsize = 2)
{
  library(keras)
  load(train_data)
  #Flank_methyl,Flank_DNA_onehot,Flank_DNA_2mer,Center_loci_table,Flank_m,Flank_um)
  table_loci <- Train_data[[4]][,c("chr","start","m.read","um.read","methylation")]
  train_array <- array(NA,dim=c(nrow(table_loci),flank_size*2,5,1))
  
  for(i in 1:nrow(table_loci))
  {
      for(j in 1:4)
      {
        train_array[i,,j,1] <- Train_data[[3]][[i]][j,]
      }
      train_array[i,,5,1] <- Train_data[[1]][i,]
  }
  
  train_y <- Train_data[[4]][,"methylation"]#downsample methyl
  
  model <- keras_model_sequential() 
  model %>%
    # Start with hidden 2D convolutional layer being fed nrow: 5,ncol: 200,1: chanel
    layer_conv_2d(filters = 32, kernel_size = c(5,5), activation = "relu",
                  input_shape = c(flank_size*2,5,1)) %>%
    layer_max_pooling_2d(pool_size = c(poolingsize, 1)) %>%
    
    layer_conv_2d(filters = 64, kernel_size = c(3, 1), activation = "relu") %>%
    layer_max_pooling_2d(pool_size = c(poolingsize, 1)) %>%
    
    layer_conv_2d(filters = 128, kernel_size = c(3, 1), activation = "relu") %>%
    layer_max_pooling_2d(pool_size = c(poolingsize, 1)) %>%
    
    layer_flatten() %>%
    layer_dense(units = 200, activation = "relu") %>%
    layer_dense(units = 1,  activation = 'linear')
  
  summary(model)
  
  
  model %>% compile(
    loss = 'mean_squared_error',
    optimizer = 'adam',
    metrics = c('mae')
  )

  history <- model %>% fit(
    train_array,train_y, shuffle=TRUE,
    epochs = epochs_Num, batch_size = batch_size_Num, 
    validation_split =  validation_split_Num,
    lrate=callback_reduce_lr_on_plateau(monitor = "val_loss",factor = 0.1,patience = 3,verbose = 2,mode = "auto",epsilon=0.0001,min_lr = 0),
    callbacks=callback_early_stopping(monitor = "val_loss",patience = 3,verbose = 2,mode = "auto")
  )
  #keras.callbacks.ReduceLROnPlateau(monitor='val_loss', factor=0.1, patience=10, verbose=0, mode='auto', epsilon=0.0001, cooldown=0, min_lr=0)
  pdf(file= paste0(figure_history_prefix,".pdf"))
  plot(history)
  Sys.sleep(2)
  dev.off()
  
  model %>% evaluate(train_array, train_y,verbose = 0)
  model %>% save_model_hdf5(save_model)
}

#' Prediction of low coverage sites with trained models.
#' 
#' @param load_model_hdf5_file The taining model file path
#' @param testfile File path of the site to be predicted. The testfile produced by "make_test_Rdata()" or "make_lowcov_test_data_from_all". A list containing Flank_methyl,Flank_DNA_onehot,Flank_DNA_2mer,Center_loci_table,Flank_m and Flank_um)
#' @param flank_size Side chain length, total window width 2*flank_size 
#' @param Predicted_file_name Output file after prediction. The data.frame contains "chr","start","m.read","um.read" and "estPredict". "estPredict" is the prediction value.
#' 
#'
#' @return
#' @export
#'
#' @examples
Predict_test_low_chr <- function(load_model_hdf5_file,testfile, flank_size,Predicted_file_name)
{
  library(keras)
  model<- load_model_hdf5(load_model_hdf5_file)
  load(testfile)
  table_loci <- Test_data[[4]][,c("chr","start","m.read","um.read")]
  test_array <- array(NA,dim=c(nrow(table_loci),flank_size*2,5,1))
  
  for(i in 1:nrow(table_loci))
  {
    for(j in 1:4)
    {
        test_array[i,,j,1] <- Test_data[[3]][[i]][j,]#bug!
      }
      test_array[i,,5,1] <- Test_data[[1]][i,]
    }
  
  estPredict <- model %>% predict(test_array)
  estPredict[estPredict < 0]  <- 0
  estPredict[estPredict > 1]  <- 1
  Prediction <- cbind(table_loci,estPredict) 
  save(Prediction,file =Predicted_file_name)
}

