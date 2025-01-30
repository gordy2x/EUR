library(rpart)
library(rpart.plot)
library(pROC)
library(neuralnet)
library(grpreg)
library(tidyverse)
library(cito)
library(sigmoid)
library(glmnet)
library(OptimalCutpoints)


# TRUE to fit and plot fit diagnostics
# FALSE to summarise models and plot results

fit <- FALSE

#################################### data  #####################################

diab_train <- read.csv("data/diabetes_clean_train.csv", stringsAsFactors =TRUE) %>%
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
         disadvantage_s = (disadvantage-1000)/50,
         log_los_days_s = log(length_of_stay_days) - 1.5) # scaling




# create predictor and response matrices
X <- model.matrix(~ age_s + sex + ethnicity + dmtypes + lcharlson + 
                    insulin + last90days + med_surg+
                    emergency + MaritalStatus + 
                    bgl_max+bgl_min + disadvantage_s + log_los_days_s,
                  data = diab_train)[,-1]

y <- matrix(diab_train$readmit)

########################## fit and check models  ##########################

if(fit){
  
  # helper function
  foo <- function(x){max(x$cvm)} # helper function
  
  ## Tree model
  
  ## use weights to balance the data
  rel_01= table(diab_train$readmit == 1)
  weights=rep(1,nrow(diab_train))
  weights[diab_train$readmit == 0]=rel_01[2]/rel_01[1]
  
  # weighted sum is the same
  sum(weights*diab_train$readmit)
  sum(weights*!diab_train$readmit)
  
  # grow a big tree
  set.seed(10)
  diab_tree_unpruned = rpart(readmit ~ age + sex + ethnicity + dmtypes + charlson +
                               insulin + last90days + med_surg+ 
                               emergency + MaritalStatus + bgl_max+bgl_min +
                               disadvantage + length_of_stay_days, 
                             data = diab_train,
                             model=T, x=TRUE,
                             method="class",
                             weights = weights,
                             control = rpart.control(maxdepth=20, minbucket = 5))
  
  
  rpart.plot(diab_tree_unpruned)
  # prune the tree using cross validation
  diab_tree = prune(diab_tree_unpruned, cp = diab_tree_unpruned$cptable[which.min(diab_tree_unpruned$cptable[,"xerror"]),"CP"])
  
  # plot the tree, and save
  rpart.plot(diab_tree, type=4, extra=7,box.palette=0)
  save(diab_tree, file ="results/diab_tree")
  
  
  
  ## Group lasso
  
  # group categorical variables, check grouping
  gp_index = as.numeric(factor(substr(colnames(X),1,6),
                               levels = unique(substr(colnames(X),1,6))))
  cbind(gp_index,colnames(X))
  
  # run group lasso model
  set.seed(3)
  diab_lasso <- cv.grpreg(X, y, group = gp_index, family = "binomial")
  
  
  # check that the cv minimisation looks okay
  plot(diab_lasso)
  
  # save model
  save(diab_lasso, file ="results/diab_lasso")
  
  
  ## elastic net
  models = list()
  alphas <- c(1, 0.5, 0.25, 0.125)
  
  # need to keep cv folds to use the same ones for different alphas
  set.seed(5)
  foldid <- sample(1:10,size = nrow(X), replace = T)
  
  # fit elastic net with each alpha value
  for(imod in 1:length(alphas)){
    models[[imod]] <- cv.glmnet(X, y,foldid = foldid,
                                alpha = alphas[imod], family = "binomial", type.measure="auc")
  }
  
  # chose best alpha
  which_best <- which.max(unlist(lapply(models,foo)))
  diab_elastic <- models[[which_best]]
  
  # check that the cv maximisation looks okay
  plot(diab_elastic)
  
  # check model
  save(diab_elastic, file ="results/diab_elastic")
  
  

  
} else {
  
  load( file ="results/diab_tree")
  load( file ="results/diab_lasso")
  load( file ="results/diab_elastic")
}

  
  ########################## model output  ##########################
  
  ## helper functions
  calc_pred_stats <- function(dat, predname){
    
    out_stats <- optimal.cutpoints(as.formula(paste(predname, "~  readmit")), 
                                   data = dat,
                                   tag.healthy = 0,
                                   methods = c("Youden"),
                                   control = control.cutpoints(generalized.Youden = TRUE,
                                                               CFN = 5,
                                                               CFP = 1))

    data.frame(AUC = summary(out_stats)$p.table$Global$AUC_CI, 
               t(unlist(out_stats$Youden$Global$optimal.cutoff)[1:10]))
    
  }
  

  
  
  ########### Calculate prediction statistics
  # tree stats

  diab_train$pred_tree = predict(diab_tree, newdata = diab_train)[,2]
  tree_stats <- calc_pred_stats(dat = diab_train, predname = "pred_tree")
    

  # lasso
  
  diab_train$pred_lasso = predict(diab_lasso, X = X, 
                                  lambda  = diab_lasso$lambda.min)
  lasso_stats <- calc_pred_stats(dat = diab_train, predname = "pred_lasso")
  
  # elastic
  
  diab_train$pred_elastic = predict(diab_elastic, newx = X, s = "lambda.min") %>% 
    as.vector()
  
  elastic_stats <- calc_pred_stats(dat = diab_train, predname = "pred_elastic") 
  

  ## nnet


  #load model weights
  layer_1_weights <- read.csv("results/model/lin1.weight.csv")
  layer_2_weights <- read.csv("results/model/lin2.weight.csv")
  layer_1_bias <- read.csv("results/model/lin1.bias.csv")
  layer_2_bias <- read.csv("results/model/lin2.bias.csv")


  
  #reformat bias matrices
  layer_1_bias_mat = as.matrix(as.data.frame(t(layer_1_bias)) %>% slice(rep(1:n(), each = nrow(X))))
  layer_2_bias_mat = as.matrix(as.data.frame(t(layer_2_bias)) %>% slice(rep(1:n(), each = nrow(X))))
  

  #evaluate the neural network
  layer_1 = X %*% t(as.matrix(layer_1_weights)) + layer_1_bias_mat
  layer_1_relu <- relu(layer_1)
  layer_2 = layer_1_relu %*% t(as.matrix(layer_2_weights)) + layer_2_bias_mat
  layer_2_sigmoid <- sigmoid(layer_2)

  diab_train$pred_nnet = layer_2_sigmoid
  
  
  nnet_stats <- calc_pred_stats(dat = diab_train, predname = "pred_nnet")


  ## combine
  perform <- bind_rows(tree = tree_stats,
                       gp_lasso = lasso_stats,
                       elastic = elastic_stats)
                       # nnet = nnet_stats)
  
  
  perform$model <- c("tree", "gp_lasso", "elastic")#, "nnet")
  
  print(perform)
  
  write.csv(perform, "results/train_pred_cut.csv", row.names = F)
  
  

