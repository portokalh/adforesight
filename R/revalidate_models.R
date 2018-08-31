#properly evaluate my models


library( knitr )
library(ANTsR)
library(visreg)
library(robustbase)
library(groupdata2)
library(ggplot2)

#reuses fused eigenregions for best performing fold
#Alex Badea 24 august 2018
# Manganese Enhanced MRI   predicts cognitive  performnce
mypath<-'/Volumes/CivmUsers/omega/alex/GitHub/adforesight/'
# Load in Behavior and Imaging Data
mypath <- '/Users/alex/GitHub/adforesight/'
setwd(mypath)


mysp <- 0.05  # 0.01  # 0.05 #0.2 #0.05 # was 0.005 sparseness
mynvecs <- 50 #put back to 50 alex #tested with 10
myell1 <- 1 # make this smaller 0.5, Brian says it is not what i think just switch between l0 and l1
myits<-5 #put back to 5
mysmooth<-0.01  # was 0
myclus<-250
extension<-paste('JAC', 'sp', toString(mysp), 'vecs', toString(mynvecs), 's', toString(mysmooth),'clus', toString(myclus), sep='') # 'JACsp0p005s0'

output.path <- paste(mypath,'/mydata/outdata/',extension, '/', sep='') #sd2_projall_noscale/'
#setwd( '/Users/omega/alex/adforesight/' )
source(paste(mypath, '/R/myr2score.R',sep=''))

#output.path <- '/Users/omega/alex/adforesight/mydata/outdata/sd2_projall_noscale/'
#if (dir.exists(output.path)){ 1} else {dir.create(output.path, recursive=TRUE)}
#if (dir.exists(paste(output.path,'/',extension, '/', sep='')){ 1} else {dir.create(paste(output.path,'/',extension, '/'), recursive=TRUE)}

# Load in Behavior and Imaging Data
behavior <- read.csv('./mydata/All_Behavior.csv')
labled.set <-read.csv('./mydata/legendsCHASS2symmetric.csv')
labeled.brain.img <-  antsImageRead('./mydata/MDT_labels_chass_symmetric.nii.gz')
mask  <- antsImageRead('./mydata/MDT_mask_e3.nii')
mask <- thresholdImage( mask, 0.1, Inf )

# mang_files <- list.files(path = "./mydata/imdata/", pattern = "T2_to_MDT",full.names = T,recursive = T)
# mang_mat <- imagesToMatrix(mang_files,mask)

jac_files <- list.files(path = "./mydata/imdata/", pattern = "jac_to_MDT",full.names = T,recursive = T)
#jac_mat <- imagesToMatrix(jac_files,mask)

#suscept_files <- list.files(path = "./mydata/imdata/", pattern = "X_to_MDT",full.names = T,recursive = T)
#suscept_mat <- imagesToMatrix(suscept_files,mask)

#pick your contrast
mang_files <- jac_files
mang_mat <- imagesToMatrix(mang_files,mask)

rows.test <- as.integer(c(2,3,4,9,15,17))
rows.train <- as.integer(c(23,8,18,13,14,19,16,12,11,24,20,7,6,1,22,10,5,21))

set.seed(1)

#https://cran.r-project.org/web/packages/groupdata2/vignettes/cross-validation_with_groupdata2.html#creating-folds-for-cross-validation  
mygroup <- behavior$genotype[1:24]
myindex <- c(1:24)

mydfb <- data.frame("mysubject_index" = factor(as.integer(myindex)),"mygenotype"=mygroup)
kable(mydfb, align = 'c')

parts <- partition(mydfb, p = c(0.25), id_col = "mysubject_index", cat_col = 'mygenotype')
parts <- partition(mydfb, p = c(0.3), id_col = "mysubject_index", cat_col = 'mygenotype')

test_set <- parts[[1]]
train_set <- parts[[2]]

# %superseeding randomness
test_set$mysubject_index <- rows.test
test_set$mygenotype <- mygroup[rows.test]
train_set$mysubject_index <- rows.train
train_set$mygenotype <- mygroup[rows.train]

# Show test_set
test_set %>% kable()
train_set %>% kable()

train_set <- fold(train_set, k = 4, cat_col = 'mygenotype', id_col = 'mysubject_index')
train_set <- train_set[order(train_set$.folds),]
train_set %>% kable()

# Order by .folds
train_set <- train_set[order(train_set$.folds),]
train_set %>% kable()

set.seed(1)
k<-4



performances <- c()
myBICs  <- c()
myR2score<-c()
myps<-c()


for (fold in 1:k){
  # for (fold in 3){
  gc(verbose = TRUE, reset = FALSE, full = TRUE)
  print('fold:',fold)
  print(fold)
  load(file=paste(output.path , "model", toString(fold), ".Rdata", sep='')) 
  # Create training set for this iteration
  training_set <- train_set[train_set$.folds != fold, ]
  # training_set %>% kable()
  testing_set <- train_set[train_set$.folds == fold, ]
  dist4predicted<-mylm$fitted.values
  dist4observed<-behavior[training_set$mysubject_index, ]
  dist4observed<-dist4observed$d4
  RMSE2<-sqrt(mean((dist4predicted- dist4observed)^2)) 
  
  performances[fold]<-RMSE2
  myBICs[fold] <- BIC(mylm)
  myR2score[fold]<-myr2score(distpred,dist4.test)
  myps[fold]<-my.p
  
 
  dist4.test <- behav.test[,'d4'] 
  
  projs.train <- data.frame(cbind(dist4.train,imgmat_train_mang)) # column combined the behavior wth the projections
  colnames(projs.train) <- c('Dist_4', paste0('Proj', c(1:ncol(imgmat_train_mang)))) # insert column names
  
  projs.test <- data.frame(cbind(dist4.test,imgmat_test_mang)) # column combind the behavior wth the projections
  colnames(projs.test) <- c('Dist_4', paste0('Proj', c(1:ncol(imgmat_test_mang)))) # insert column names
  
  mylm2 <- lm('Dist_4 ~ .', data=projs.train) # behavior correlation with the number of projections
  
  distpred <- predict.lm(mylm, newdata=projs.test) # based on the linear model predict the distances for the same day
  modsum <-summary(lm(distpred~dist4.test))
  r2 <- modsum$adj.r.squared
  my.p <- modsum$coefficients[2,4]
  plot(dist4.test, distpred, xlab = 'Real Swim Distance on Day 4', ylab = 'Predicted Swim Distance on day 4', 
       main='Predicted vs. Real Swim Distance on Day 4', ylim = c(0,1200),xlim = c(0,1200)) # generate plot
  
  mymodel<-lm(distpred~dist4.test)
  
  
  
  RSS <- c(crossprod(mymodel$residuals))
  MSE <- RSS / length(mymodel$residuals)
  RMSE <- sqrt(MSE)
  
  RMSE2<-sqrt(mean((distpred - dist4.test)^2)) 
  
  performances[fold]<-RMSE2
  myBICs[fold] <- BIC(mymodel)
  myR2score[fold]<-myr2score(distpred,dist4.test)
  myps[fold]<-my.p
  
  
 
  
  
}
