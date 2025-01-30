
library(rpart)
library(rpart.plot)
library(pROC)
library(neuralnet)
library(grpreg)
library(tidyverse)
library(cito)
library(sigmoid)
library(glmnet)
library(MLmetrics)
library(carat)


# which data to calculate statistics for
dat_source = "validation"# "test", "train" or "validation", "validation_up", "test_up"



# read in cutoffs for classification we calculated for training data
cutoffs <-   read.csv("data/derived/train_pred_cut.csv") %>% 
  dplyr::select(model, cutoff)


# read in data
diab_data <- read.csv(paste0("data/diabetes_clean_",dat_source,".csv"), 
                       stringsAsFactors =TRUE) %>% 
  mutate(charlson = exp(lcharlson) - 1) %>% 
  dplyr::select( -patient_mrn, -ID, -admit_dttm, -osp_full_name, - represent) %>% 
  relocate(readmit) %>% # put readmit in first column
  mutate(age = as.integer(age), # for tree model
         bgl_min=relevel(bgl_min, ">=4"),  
         med_surg = factor(med_surg, levels = 1:2, labels = c("Medical", "Surgical")),
         emergency = factor(emergency, levels = 1:2, labels = c("Emerg", "NotEmerg")),
         emergency = relevel(emergency, 2),
         charlson = as.integer(charlson), # for tree model
         disadvantage = as.integer(disadvantage), # for tree model
         length_of_stay_days = as.integer(length_of_stay_days),
         ethnicity = relevel(ethnicity, "aust"),
         age_s = (age - 70)/10, # scaling
         lcharlson_s = 2*(lcharlson - 0.5), # scaling
         disadvantage_s = (disadvantage-1000)/50,
         log_los_days_s = log(length_of_stay_days) - 1.5) # scaling

# create predictor matrix, needed for prediction
X <- model.matrix(~ age_s + sex + ethnicity + dmtypes + lcharlson + 
                    insulin + last90days + med_surg+
                    emergency + MaritalStatus + 
                    bgl_max+bgl_min + disadvantage_s + log_los_days_s,
                  data = diab_data)[,-1]


# load fitted models
load( file ="data/derived/diab_tree")
load( file ="data/derived/diab_lasso")
load( file ="data/derived/diab_elastic")

# #load model weights for nnet
layer_1_weights <- read.csv("data/derived/model/lin1.weight.csv")
layer_2_weights <- read.csv("data/derived/model/lin2.weight.csv")
layer_1_bias <- read.csv("data/derived/model/lin1.bias.csv")
layer_2_bias <- read.csv("data/derived/model/lin2.bias.csv")


###################### calculate model fit ####################

# calculating PPV and NPV
# https://www.ncbi.nlm.nih.gov/books/NBK430867/
prevalence = mean(diab_data$readmit)

# tree
diab_data$pred_tree = predict(diab_tree, newdata = diab_data)[,2]
diab_data$pred_tree_bin <- (diab_data$pred_tree > cutoffs$cutoff[cutoffs$model == "tree"])*1

tree_stats = caret::confusionMatrix(data = factor(diab_data$pred_tree_bin), 
                                    factor(diab_data$readmit),
                             positive = "1")
roc_all <- roc(diab_data$readmit,diab_data$pred_tree,ci=TRUE,plot=FALSE)
# write.csv(roc_all, file = paste0("data/derived/roc_tab_tree_",dat_source,".csv" ))



auc_tab = data.frame(score = roc_all$thresholds,
                     sensitivities = roc_all$sensitivities,
                     specificities = roc_all$specificities) %>% 
  mutate(PPV = (sensitivities * prevalence) / ( (sensitivities * prevalence) + ((1 - specificities) * (1 - prevalence))),
         NPV = (specificities * (1 - prevalence)) / ( (specificities * (1 - prevalence)) + ((1 - sensitivities) * prevalence) )) 



write.csv(auc_tab, file = paste0("data/derived/auc_tab_tree_",dat_source,".csv" ))

tree_roc <- roc_all$ci %>% 
  as.vector()

# lasso
diab_data$pred_lasso = predict(diab_lasso, X = X, 
                                lambda  = diab_lasso$lambda.min)
diab_data$pred_lasso_bin <- (diab_data$pred_lasso > cutoffs$cutoff[cutoffs$model == "gp_lasso"])*1

lasso_stats = caret::confusionMatrix(data = factor(diab_data$pred_lasso_bin), 
                                    factor(diab_data$readmit),
                                    positive = "1")

roc_all <- roc(diab_data$readmit,diab_data$pred_lasso,ci=TRUE,plot=FALSE)
# write.csv(roc_all, file = paste0("data/derived/roc_tab_lasso_",dat_source,".csv" ))

auc_tab = data.frame(score = roc_all$thresholds,
                     sensitivities = roc_all$sensitivities,
                     specificities = roc_all$specificities) %>% 
  mutate(PPV = (sensitivities * prevalence) / ( (sensitivities * prevalence) + ((1 - specificities) * (1 - prevalence))),
         NPV = (specificities * (1 - prevalence)) / ( (specificities * (1 - prevalence)) + ((1 - sensitivities) * prevalence) )) 

write.csv(auc_tab, file = paste0("data/derived/auc_tab_lasso_",dat_source,".csv" ))



lasso_roc <- roc_all$ci %>% 
  as.vector()




#elastic

diab_data$pred_elastic = predict(diab_elastic, newx = X, s = "lambda.min") %>% 
  as.vector()
diab_data$pred_elastic_bin <- (diab_data$pred_elastic > cutoffs$cutoff[cutoffs$model == "elastic"])*1

elastic_stats = caret::confusionMatrix(data = factor(diab_data$pred_elastic_bin), 
                                     factor(diab_data$readmit),
                                     positive = "1")

roc_all <- roc(diab_data$readmit,diab_data$pred_elastic,ci=TRUE,plot=FALSE)
# write.csv(roc_all, file = paste0("data/derived/roc_tab_elastic_",dat_source,".csv" ))
# 

auc_tab = data.frame(score = roc_all$thresholds,
                     sensitivities = roc_all$sensitivities,
                     specificities = roc_all$specificities) %>% 
  mutate(PPV = (sensitivities * prevalence) / ( (sensitivities * prevalence) + ((1 - specificities) * (1 - prevalence))),
         NPV = (specificities * (1 - prevalence)) / ( (specificities * (1 - prevalence)) + ((1 - sensitivities) * prevalence) )) 


write.csv(auc_tab, file = paste0("data/derived/auc_tab_elastic_",dat_source,".csv" ))

elastic_roc <- roc_all$ci %>% 
  as.vector()


#reformat bias matrices
layer_1_bias_mat = as.matrix(as.data.frame(t(layer_1_bias)) %>% slice(rep(1:n(), each = nrow(X))))
layer_2_bias_mat = as.matrix(as.data.frame(t(layer_2_bias)) %>% slice(rep(1:n(), each = nrow(X))))

#evaluate the neural network
layer_1 = X %*% t(as.matrix(layer_1_weights)) + layer_1_bias_mat
layer_1_relu <- relu(layer_1)
layer_2 = layer_1_relu %*% t(as.matrix(layer_2_weights)) + layer_2_bias_mat
diab_data$pred_nnet <- as.vector(sigmoid(layer_2))


# stats
# diab_data$pred_nnet_bin <- (diab_data$pred_nnet > cutoffs$cutoff[cutoffs$model == "nnet"])*1

# nnet_stats = caret::confusionMatrix(data = factor(diab_data$pred_nnet_bin), 
#                                        factor(diab_data$readmit),
#                                     positive = "1")

roc_all <- roc(diab_data$readmit,diab_data$pred_nnet,ci=TRUE,plot=FALSE) 

auc_tab = data.frame(score = roc_all$thresholds,
                     sensitivities = roc_all$sensitivities,
                     specificities = roc_all$specificities) %>% 
  mutate(PPV = (sensitivities * prevalence) / ( (sensitivities * prevalence) + ((1 - specificities) * (1 - prevalence))),
         NPV = (specificities * (1 - prevalence)) / ( (specificities * (1 - prevalence)) + ((1 - sensitivities) * prevalence) )) 


nnet_roc <- roc_all$ci %>% 
  as.vector()

write.csv(auc_tab, file = paste0("data/derived/auc_tab_nnet_",dat_source,".csv" ))


nnet_roc <- roc(diab_data$readmit,diab_data$pred_nnet,ci=TRUE,plot=FALSE)$ci %>% 
  as.vector()


## combine
AUC <- data.frame(rbind(tree_roc[c(2,1,3)], 
                  lasso_roc[c(2,1,3)], 
                  elastic_roc[c(2,1,3)], 
                  nnet_roc[c(2,1,3)]) )
colnames(AUC) <- c("AUC","LCL", "UCL")

AUC$model <- c("tree", "gp_lasso", "elastic",  "nnet")#


write.csv(AUC, file = paste0("data/derived/model_pred_",dat_source, ".csv" ),
          row.names = F)

write.csv(diab_data, file = paste0("data/derived/data_with_pred_",dat_source, ".csv" ),
          row.names = F)

