### read in dependencies ----------------------------------
library(geomorph)
library(car)
library(Morpho)
require(MASS)     # contains lda() and qda() functions
require(sampling) # easier way to generate subsets of data
require(mda)      # flexible discriminant analysis and mixture discriminant analysis
require(class)    # k nearest neighbours (knn)
require(adabag)   # bagging()
require(ada)      # boosting function, ada()
require(nnet)     # neural net function, nnet()
require(e1071)    # svm()
require(randomForest) # randomForest() 
require(plyr)
#########################################################################################

### set working directory at images folder  ------------------------
setwd("Documents/DSD Shape Data/")  #work computer
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/DSD Shape Data/")   #home computer

### read in source files ----------------------------
source("scripts/defineSliders.R")   #contains semi-landmark and link configurations
source("scripts/NDT_functions.R")   #contains custom scripts for analysis
source("scripts/WRP_FUNCTIONS.R")
source("scripts/dataRead_40x.R")        #reads in raw data, performs GPA, and other useful data processing


### Common Functions ----------------
#creating testing and training sets, genotype
#########################################################################################
strata_5var <- function(dataset, variable_name, variable_table, percentage, variable_pos, dataframe_length = 73){
  #this function stratifies a dataset by five variables, 
  #and makes testing and training sets based on the strata
  #testing and training sets returned as a list
  start_position <- dataframe_length - 35
  my_strata <- strata(dataset, stratanames = c(variable_name), method = "srswor", 
                      size = c(round(variable_table[1] *(percentage)), 
                               round(variable_table[2] * (percentage)),
                               round(variable_table[3] * (percentage)),
                               round(variable_table[4] * (percentage)),
                               round(variable_table[5] * (percentage)),
                               round(variable_table[6] * (percentage))))
  my_data <- dataset[, c(variable_pos, start_position:dataframe_length)]
  my_training <- my_data[rownames(my_data) %in% my_strata$ID_unit, ]
  my_testing <- my_data[!rownames(my_data) %in% my_strata$ID_unit, ] 
  my_train_test_info <- list(my_training, my_testing)
  return(my_train_test_info)
}

strata_4var <- function(dataset, variable_name, variable_table, percentage, variable_pos, dataframe_length = 73){
  #this function stratifies a dataset by two variables, 
  #and makes testing and training sets based on the strata
  #testing and training sets returned as a list
  start_position <- dataframe_length - 35
  my_strata <- strata(dataset, stratanames = c(variable_name), method = "srswor", 
                      size = 
                        c(round(variable_table[1] *(percentage)), 
                          round(variable_table[2] * (percentage)),
                          round(variable_table[3] *(percentage)), 
                          round(variable_table[4] * (percentage))))
  my_data <- dataset[, c(variable_pos, start_position:dataframe_length)]
  my_training <- my_data[rownames(my_data) %in% my_strata$ID_unit, ]
  my_testing <- my_data[!rownames(my_data) %in% my_strata$ID_unit, ] 
  my_train_test_info <- list(my_training, my_testing)
  return(my_train_test_info)
}

strata_2var <- function(dataset, variable_name, variable_table, percentage, variable_pos, dataframe_length = 73){
  #this function stratifies a dataset by two variables, 
  #and makes testing and training sets based on the strata
  #testing and training sets returned as a list
  start_position <- dataframe_length - 35
  my_strata <- strata(dataset, stratanames = c(variable_name), method = "srswor", size = 
                        c(round(variable_table[1] *(percentage)), round(variable_table[2] * (percentage))))
  my_data <- dataset[, c(variable_pos, start_position:dataframe_length)]
  my_training <- my_data[rownames(my_data) %in% my_strata$ID_unit, ]
  my_testing <- my_data[!rownames(my_data) %in% my_strata$ID_unit, ] 
  my_train_test_info <- list(my_training, my_testing)
  return(my_train_test_info)
}

resample_default <- function(function_name, dataset, datatable, reps){
  out_list <- list()
  for(i in 1:reps){
    myresult <- function_name(dataset, datatable)
    out_list <- unlist(c(out_list, myresult))
  }
  mean_score <- mean(out_list)
  error_score <- sd(out_list)
  return_list <- list(mean_score,error_score,out_list) 
  return(return_list)
}

myHeatMap <- function( dim, names, data, xLab="Actual", yLab="Predicted", Main="" ){
  myCols <- colorRampPalette(c("#ffffff7f", "#ffff007f", "#ff00007f"), space="rgb")
  image( 1:dim, 1:dim, data[,dim:1], col=myCols(100), xaxt='n', yaxt='n', xlab=xLab, ylab=yLab, main=Main )
  abline(h=0:dim+0.5, col="grey")
  abline(v=0:dim+0.5, col="grey")
  box( lwd=1 )
  text( 1:dim, rep(dim:1, each=dim), sub('^0$', '',data) )
  axis(1, at=1:dim, labels=names, cex.axis=0.8 )
  axis(2, at=1:dim, labels=rev(names), cex.axis=0.8, las=1 )
}


### calculate variables / dataframes ---------------- 
wings <- data.frame("genotype"= paste(ID, devTemp, sep="-"), two.d.cell)
#wings <- data.frame("genotype"= ID, two.d.cell)   # use this when only looking at species
#wing_average <- aggregate( wings[,1:36], list( ID=wings$ID, genotype=wings$genotype, sex=wings$sex), mean) #this is to average L and R wings from each fly, so not useful here
wing_average <- wings
wings_mean_pc <- prcomp(wing_average[,2:37])
average_wings <- data.frame(wing_average, wings_mean_pc$x[,1:36])
tabw_g <- table(average_wings$genotype)


Names <- c( expression( paste( italic( "AF16-20" ))),
            expression( paste( italic( "AF16-30" ))),
            expression( paste( italic( "N2-20" ))),
            expression( paste( italic( "N2-30" ))) )



### Neural Net --------------

#neural network  | to look at just N2 vs AF16, change strata_5var to strata_2var below
nnet_repeat_function_genotype <- function(dataset, datatable){
  nnet_wings_genotype <- strata_4var(dataset = dataset, variable_name = "genotype", variable_table = datatable, 
                                     percentage = 2/3, variable_pos = 1, dataframe_length = 73)
  train_set <- nnet_wings_genotype[[1]]
  test_set <- nnet_wings_genotype[[2]]
  nnet_model <- nnet(formula=genotype ~., data = train_set, size=10, decay=0)
  nnet_training_table <- table(actual = train_set$genotype,
                               predicted = predict(nnet_model, type="class"))                         
  nnet_test_table <- table(actual = test_set$genotype,
                           predicted = predict(nnet_model, newdata=test_set, type="class"))
  prediction_success <- 100*sum(diag(nnet_test_table)/sum(nnet_test_table))
  return(prediction_success)
}

nnet_score_CI_genotype <- resample_default(function_name = nnet_repeat_function_genotype, 
                                           dataset = average_wings, 
                                           datatable = tabw_g, 
                                           reps = 100) 
#returns values for [[1]] mean score (mean model precision), [[2]] error score (standard deviation or precision), and [[3]] list of individual model precision values
nnet_score_CI_genotype


nnet_genotype_confusion <- function(dataset, datatable){
  nnet_wings_genotype <- strata_4var(dataset = dataset, variable_name = "genotype", variable_table = datatable, 
                                     percentage = 2/3, variable_pos = 1, dataframe_length = 73)
  train_set <- nnet_wings_genotype[[1]]
  test_set <- nnet_wings_genotype[[2]]
  nnet_model <- nnet(formula=genotype ~., data = train_set, size=10, decay=0)
  nnet_training_table <- table(actual = train_set$genotype,
                               predicted = predict(nnet_model, type="class"))                         
  nnet_test_table <- table(actual = test_set$genotype,
                           predicted = predict(nnet_model, newdata=test_set, type="class"))
  #prediction_success <- 100*sum(diag(nnet_test_table)/sum(nnet_test_table))
  return(nnet_test_table)
}
### landmark data
wings_info <- nnet_genotype_confusion( average_wings, tabw_g)

wings_info_perc <- 100*( wings_info / rowSums( wings_info ) )
myHeatMap( dim=4, names=Names, data=round(wings_info_perc, 0), Main="Neural Net" )


image( 0:5, 0:1, matrix( c(100,80,60,40,20,0), byrow=T ), xaxt='n', yaxt='n', xlab='', ylab='', main="% correct" )
text( 0:5, 0.5, c(0,20,40,60,80,100) )
abline( v=0.5:5.5, col="grey")
box( lwd=1 )


### Random Forest ------------

random_forest_genotype <- function(train_data, test_data){
  #random forest predictions
  #makes a training and testing table, outputs a list including
  #training table, test table, and overall prediction success percentage
  vector_of_trees <- c(10, 100, 500, 1000)
  random_forest_return <- list()
  final_output_val <- 0
  final_setting <- ""
  for (i in vector_of_trees)
  {
    random_forest_model <- randomForest(genotype ~., data = train_data, ntree = i)
    rf_training_table <- table(actual = train_data$genotype,
                               predicted = predict(random_forest_model, type="class"))
    rf_test_table <- table(actual = test_data$genotype,
                           predicted = predict(random_forest_model, newdata=test_data, type="class"))
    prediction_success <- 100 * sum(diag(rf_test_table)/sum(rf_test_table))
    if (prediction_success > final_output_val){final_output_val <- prediction_success; final_setting <- i}
    random_forest_return_part <- list(rf_training_table, rf_test_table, prediction_success, random_forest_model, i)
  }
  random_forest_return <- list(final_output_val, final_setting)
  return(random_forest_return)
}



#function for looping the random forest genotype function with different testing/training sets
rf_genotype_fun <- function(dataset, datatable){
  rf_wings_genotype <- strata_5var(dataset, "genotype", datatable, 2/3, 1, 73)
  rf_training <- rf_wings_genotype[[1]]
  rf_test <- rf_wings_genotype[[2]]
  mytest <- random_forest_genotype(rf_training, rf_test)
  myreturn <- mytest
  return(myreturn)
}


resample_settings_rf <- function(function_name, dataset, datatable, reps){
  mylist_10 <- list()
  mylist_100 <- list()
  mylist_500 <- list()
  mylist_1000 <- list()
  out_list <- list()
  for(i in 1:reps){
    myresult <- function_name(dataset, datatable)
    if (myresult[[2]] == 10) {mylist_10 <- unlist(c(mylist_10, myresult[[1]]))}
    if (myresult[[2]] == 100) {mylist_100 <- unlist(c(mylist_100, myresult[[1]]))}
    if (myresult[[2]] == 500) {mylist_500 <- unlist(c(mylist_500, myresult[[1]]))}
    if (myresult[[2]] == 1000) {mylist_1000 <- unlist(c(mylist_1000, myresult[[1]]))}}
  out_list <- list(mylist_10, mylist_100, mylist_500, mylist_1000)
  return(out_list)
}


rf_list_genotype <- resample_settings_rf(rf_genotype_fun, average_wings, tabw_g, 100)
rf_list_genotype #500 trees works the best here | lists results from 10, 100, 500, and 1000 random tree models


rf_genotype_confusion <- function(dataset, datatable){
  trees = list(1000)
  rf_wings_genotype <- strata_4var(dataset, "genotype", datatable, 2/3, 1, 73) #different dimensions, not average_wings dataframe
  rf_training <- rf_wings_genotype[[1]]
  rf_test <- rf_wings_genotype[[2]]
  results_list <- list()
  for (i in trees){
    random_forest_model <- randomForest(genotype ~., data = rf_training, ntree = i)
    rf_test_table <- table(actual = rf_test$genotype,
                           predicted = predict(random_forest_model, newdata=rf_test, type="class"))
  }
  return(rf_test_table)
}

wings_info <- rf_genotype_confusion( average_wings, tabw_g)

wings_info_perc <- 100*( wings_info / rowSums( wings_info ) )
myHeatMap( dim=4, names=Names, data=round(wings_info_perc, 0), Main="Random Forest" )



### Support Vector Machine (SVM) ------

svm_genotype_fun <- function(dataset, datatable, kernel_type = "radial"){
  svm_wings_genotype <- strata_4var(dataset, "genotype", datatable, 2/3, 1, 73)
  svm_training <- svm_wings_genotype[[1]]
  svm_test <- svm_wings_genotype[[2]]
  svm_model <- svm(genotype ~., data = svm_training, kernel = kernel_type)
  svm_training_table <- table(actual = svm_training$genotype,
                              predicted = predict(svm_model, type="class"))
  svm_test_table <- table(actual = svm_test$genotype,
                          predicted = predict(svm_model, newdata=svm_test, type="class"))
  prediction_success <- 100 * sum(diag(svm_test_table)/sum(svm_test_table))
  return(prediction_success)
}

svm_genotype_confusion <- function(dataset, datatable, kernel_type = "radial"){
  svm_wings_genotype <- strata_4var(dataset, "genotype", datatable, 2/3, 1, 73)
  svm_training <- svm_wings_genotype[[1]]
  svm_test <- svm_wings_genotype[[2]]
  svm_model <- svm(genotype ~., data = svm_training, kernel = kernel_type)
  svm_training_table <- table(actual = svm_training$genotype,
                              predicted = predict(svm_model, type="class"))
  svm_test_table <- table(actual = svm_test$genotype,
                          predicted = predict(svm_model, newdata=svm_test, type="class"))
  #prediction_success <- 100 * sum(diag(svm_test_table)/sum(svm_test_table))
  return(svm_test_table)
}

svm_genotype_info <- resample_default(svm_genotype_fun, average_wings, tabw_g, 5)
svm_genotype_info

wings_info <- svm_genotype_confusion( average_wings, tabw_g)
wings_info_perc <- 100*( wings_info / rowSums( wings_info ) )
myHeatMap( dim=4, names=Names, data=round(wings_info_perc, 0), Main="Random Forest" )



### Combined Confusion Matrix Heatmaps -----------
# this works FINE, but I don't want to lose the code... layout is much better
# split.screen( figs = c( 2, 1 ) )
# split.screen( figs = c( 1, 2 ), screen = 1 )
# 
# screen(2)
# image( 0:5, 0:1, matrix( c(100,80,60,40,20,0), byrow=T ), xaxt='n', yaxt='n', xlab='', ylab='', main="% correct" )
# text( 0:5, 0.5, c(0,20,40,60,80,100) )
# abline( v=0.5:5.5, col="grey")
# box( lwd=1 )
# 
# screen(3)
# par(mar = c(0,4,3,0))
# wings_info <- rf_genotype_confusion( average_wings, tabw_g)
# wings_info_perc <- 100*( wings_info / rowSums( wings_info ) )
# myHeatMap( dim=4, names=Names, data=round(wings_info_perc, 0), Main="Random Forest" )
# 
# screen(4)
# par(mar = c(0,4,3,0))
# wings_info <- nnet_genotype_confusion( average_wings, tabw_g)
# wings_info_perc <- 100*( wings_info / rowSums( wings_info ) )
# myHeatMap( dim=4, names=Names, data=round(wings_info_perc, 0), Main="Neural Net", yLab = "")
# 
# close.screen(all.screens = TRUE)

## OR ##
layout(matrix(c(1,2,3,4,4,5), 2, 3, byrow = TRUE), widths=c(2,2,2), heights=c(4,1.5))
par(mar = c(4,4,3,3))
wings_info <- rf_genotype_confusion( average_wings, tabw_g)
wings_info_perc <- 100*( wings_info / rowSums( wings_info ) )
myHeatMap( dim=4, names=Names, data=round(wings_info_perc, 0), Main="Random Forest", xLab = "Actual" )
wings_info <- nnet_genotype_confusion( average_wings, tabw_g)
wings_info_perc <- 100*( wings_info / rowSums( wings_info ) )
myHeatMap( dim=4, names=Names, data=round(wings_info_perc, 0), Main="Neural Net", yLab = "", xLab = "Actual")
wings_info <- svm_genotype_confusion( average_wings, tabw_g)
wings_info_perc <- 100*( wings_info / rowSums( wings_info ) )
myHeatMap( dim=4, names=Names, data=round(wings_info_perc, 0), Main="SVM" , yLab = "", xLab = "Actual")
image( 0:5, 0:1, matrix( c(100,80,60,40,20,0), byrow=T ), xaxt='n', yaxt='n', xlab='', ylab='', main="% correct" )
text( 0:5, 0.5, c(0,20,40,60,80,100) )
abline( v=0.5:5.5, col="grey")
box( lwd=1 )
