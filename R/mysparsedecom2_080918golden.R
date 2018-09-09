#load libs
library( knitr )
library(ANTsR)
library(visreg)
library(robustbase)
library(groupdata2) 
library(ggplot2)

#Manganese Enhanced MRI   predicts cognitive  performnace Alex Badea and Natalie Delpratt
#reuses results of sparsedecom2 for best performing fold Alex Badea 7 dept 2018
#uses RMSE as prediction error, rather than goodness of fit error 30 August 2018

#legacy from Natalie to rememebr paths
#setwd( '/Users/omega/alex/adforesight/' )
#output.path <- '/Users/omega/alex/adforesight/mydata/outdata/sd2_projall_noscale/' #ND


mypath<-'/Volumes/CivmUsers/omega/alex/GitHub/adforesight/'
mypath <- '/Users/alex/GitHub/adforesight/' #flavors of serifos
mypath <- '/Users/alex/Documents/GitHub/adforesight/' #flavors of ithaka
setwd(mypath)
source(paste(mypath, '/R/myr2score.R',sep=''))

mysp <- 0.02 #0.05  # 0.01  # 0.05 #0.2 #0.05 # was 0.005 sparseness
mynvecs <- 2 # 10 shows ventral thalamic nuclei 10 # 50 #put back to 50 alex #tested with 10
myell1 <- 1 # make this smaller 0.5, Brian says it is not what i think just switch between l0 and l1
myits<-15 #5 #put back to 5
mysmooth<-0.1 #0.01 #0  # was 0.01
myclus<-250 #was 250

#build a place to save results
extension<-paste('CHI', 'sp', toString(mysp), 'vecs', toString(mynvecs), 's', toString(mysmooth),'clus', toString(myclus), sep='') # 'JACsp0p005s0'
output.path <- paste(mypath,'/mydata/outdata_sd2/',extension, '/', sep='') #sd2_projall_noscale/'
if (dir.exists(output.path)){ 1} else {dir.create(output.path, recursive=TRUE)}

# Load in Behavior and Imaging Data
behavior <- read.csv('./mydata/All_Behavior.csv')
labled.set <-read.csv('./mydata/legendsCHASS2symmetric.csv')
labeled.brain.img <-  antsImageRead('./mydata/MDT_labels_chass_symmetric.nii.gz')
mask  <- antsImageRead('./mydata/MDT_mask_e3.nii')
mask <- thresholdImage( mask, 0.1, Inf )

mang_files <- list.files(path = "./mydata/imdata/", pattern = "T2_to_MDT",full.names = T,recursive = T)
jac_files <- list.files(path = "./mydata/imdata/", pattern = "jac_to_MDT",full.names = T,recursive = T)
chi_files <- list.files(path = "./mydata/imdata/", pattern = "X_to_MDT",full.names = T,recursive = T)

#pick your contrast
mang_files <- chi_files
mang_mat <- imagesToMatrix(mang_files,mask)

rows.test <- as.integer(c(2,3,4,9,15,17))
rows.train <- as.integer(c(23,8,18,13,14,19,16,12,11,24,20,7,6,1,22,10,5,21))

set.seed(1)
res<-createFolds(behavior$genotype,4)
behavior$genotype[res$Fold1]

#https://cran.r-project.org/web/packages/groupdata2/vignettes/cross-validation_with_groupdata2.html#creating-folds-for-cross-validation  
mygroup <- behavior$genotype[1:24]
myindex <- c(1:24)

mydfb <- data.frame("mysubject_index" = factor(as.integer(myindex)),"mygenotype"=mygroup)
kable(mydfb, align = 'c')

#parts <- partition(mydfb, p = c(0.25), id_col = "mysubject_index", cat_col = 'mygenotype')
parts <- partition(mydfb, p = c(0.3), id_col = "mysubject_index", cat_col = 'mygenotype')

test_set <- parts[[1]]
train_set <- parts[[2]]

# this superseeds randomness for now
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
gfit<-c()

###build k models and retain the best performing one in terms of RMSE2
#LOOCV replaces folds
#k<-length(rows.train)-1

k<-4
res_train<-createFolds(behavior$genotype,k, list = TRUE, returnTrain = TRUE)
res_test<-createFolds(behavior$genotype,k)

for (fold in 1:k){
 # for (fold in 3){
  gc(verbose = TRUE, reset = FALSE)
  print('fold:',fold)
  print(fold)

 
  # training_set <- train_set[!(rows.train %in% rows.train[fold]), ] ### Create training set for this iteration
  # testing_set <- train_set[fold, ] ### training_set %>% kable() #showoffish
  # testing_set %>% kable() #showoffish
  
  rows.train<-as.integer(res_train[[fold]])
  rows.test<-as.integer(res_test[[fold]])
   
  mang.train <- mang_mat[rows.train, ]
  mang.test <- mang_mat[rows.test, ]
  behav.train <- behavior[rows.train, ]
  behav.test  <- behavior[rows.test, ]
  
  
  start_time <- Sys.time()
 
  myeig2<-sparseDecom2(inmatrix = list(mang.train,as.matrix(behav.train$d4)),its = 10, cthresh=c(myclus,0), smooth = mysmooth, mycoption = 0, sparseness = c(mysp,1), nvecs = mynvecs, verbose=1, statdir=paste(output.path2))
  end_time <- Sys.time()
  t1time<-end_time - start_time
  print(t1time)
  
  imgpredtrain_mang<-mang.train %*% (myeig2$eig1)
  imgmat_train_mang<-imgpredtrain_mang # just keep one of these
  imgpredtest_mang<-mang.test %*% (myeig2$eig1)
  imgmat_test_mang<-imgpredtest_mang
  
  dist4.train <- behav.train[,'d4']
  dist4.test <- behav.test[,'d4'] 
  
  projs.train <- data.frame(cbind(dist4.train,imgmat_train_mang)) # column combined the behavior with the projections
  projs.test <- data.frame(cbind(dist4.test,imgmat_test_mang)) # column combind the behavior with the projections
  colnames(projs.train) <- c('Dist_4', paste0('Proj', c(1:ncol( imgmat_train_mang))))  # insert column names
  colnames(projs.test)  <- c('Dist_4', paste0('Proj', c(1:ncol(imgmat_test_mang)))) # insert column names
  
  mylm <- lm('Dist_4 ~ .', data=projs.train) # behavior correlation with projections
 
  e2i<-matrixToImages(t((myeig2$eig1)),mask = mask)
  
  
  for (i in 1:mynvecs){
    antsImageWrite(e2i[[i]],paste(output.path,extension,'sd2eig' ,as.character(i), 'fold', toString(fold), '.nii.gz',sep=''))
    }
  
  distpred4 <- predict.lm(mylm, newdata=projs.test) # based on the linear model predict the distances for the same day
  #remove next lines for LOOCV
  mymodel<-lm(distpred4~dist4.test)
  modsum <-summary(mymodel)
  r2 <- modsum$r.squared #modsum$adj.r.squared
  my.p <- modsum$coefficients[2,4]
  RMSE2<-sqrt(mean((distpred4 - dist4.test)^2)) 
  performances[fold]<-RMSE2
  myR2score[fold]<-myr2score(distpred4,dist4.test)
  myps[fold]<-my.p
  myBICs[fold] <- BIC(mymodel)
  
  
  ###
  mytheme <- theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                   panel.background = element_rect(fill = "white"))
  myplot<- visreg(mymodel, gg=TRUE)
  myplot + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                 panel.background = element_rect(fill = "transparent", colour = NA),
                 axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) +
  ggtitle(paste("RMSE=",formatC(RMSE2,digits=2, format="f"), "p=",formatC(myps[fold],digits=2, format="f"), " BIC=", formatC(BIC(mymodel),digits=2, format="f")))

  ggsave(paste(output.path,extension,'Mnfold',toString(fold),'.pdf',sep=''), plot = last_plot(), device = 'pdf',
         scale = 1, width = 4, height = 4, units = c("in"),dpi = 300)

  save(mylm, file=paste(output.path , "model2", toString(fold), ".Rdata", sep=''))
  save(mymodel, file=paste(output.path , "behavmodelsd2", toString(fold), ".Rdata", sep=''))
  myperf<-data.frame(rbind(distpred4,dist4.test),row.names=c("d_predicted","d_valid"))
  write.csv(myperf, file = paste(output.path ,extension,'distances4_pv_fold' , toString(fold), '.csv',sep=''))
  myperf<-data.frame(c(RMSE2,myR2score[fold],myps[fold],myBICs[fold], r2),row.names=c("RMSE2","R2score","p","BIC", "R2"))
  write.csv(myperf, file = paste(output.path ,extension,'distances4_stats_fold' , toString(fold), '.csv',sep=''))
  #plot(dist4.test, distpred4, xlab = 'Real Swim Distance on Day 4', ylab = 'Predicted Swim Distance on day 4', 
  #      main='Predicted vs. Real Swim Distance on Day 4', ylim = c(0,1200),xlim = c(0,1200)) # generate plot
  # 
  
  for ( bestpred in 1:ncol(imgpredtrain_mang)) {
    cogpredtrain_mang <-behav.train$d4 %*% t(as.matrix(myeig2$eig2)[,bestpred])
    cogpredtest_mang <- behav.test$d4 %*% t(as.matrix(myeig2$eig2)[,bestpred])
    gvars<-paste("Proj",c(1:nrow(myeig2$eig2)),sep='',collapse='+')
    projs.train <- data.frame(cbind(cogpredtrain_mang,imgpredtrain_mang)) 
    colnames(projs.train) <- c('cognitive',paste0( 'Proj', c( 1:ncol(imgpredtrain_mang ) )))
    projs.test <- data.frame(cbind(cogpredtest_mang, imgpredtest_mang))
    colnames(projs.test) <- c('cognitive',paste0( 'Proj', c( 1:ncol(imgpredtest_mang ) )))
    myform<-as.formula( paste("cognitive~",gvars,sep='') )
    mylm <- lm(myform, data=projs.train) 
    distpred <- bigLMStats( mylm)
    cat(paste("Eig",bestpred,"is related to:\n"))
    #mycog<-colnames(behav.train.mat)[ abs(myeig2$eig2[,bestpred]) > 0 ]
    mycog<-colnames(behav.train)[5]
    cat( mycog )
    cat("\nwith weights\n")
    cat( abs(myeig2$eig2[,bestpred])[ abs(myeig2$eig2[,bestpred]) > 0 ])
    cat(paste("\nwith predictive correlation:", cor( cogpredtest_mang,predict(mylm,newdata=projs.test))))
      }
  
  #gc(verbose = TRUE, reset = FALSE, full = TRUE)
  gc(verbose = TRUE, reset = FALSE)
 
  }


###################################
#####    for validation now    ####
###################################

myminfold<-which(performances == min(performances), arr.ind = TRUE)
fold<-myminfold
load(file=paste(output.path , "model2", toString(fold), ".Rdata", sep='')) # loads mylm
#get model for fold 3 or whatever the minimum is

#besteigfiles<-list.files(path=paste(output.path, sep=''), pattern = paste('fold',toString(fold),'.nii.gz',sep=''),full.names = T,recursive = T)
#best_eig_jac_mat <- imagesToMatrix(besteigfiles,mask)
#eigsum<-antsAverageImages(besteigfiles,NORMALIZe=FALSE)
#antsImageWrite(eigsum,paste(output.path,extension,'EigSUM' ,as.character(i), 'fold', toString(fold), '.nii.gz',sep=''))
#myeig<-imagesToMatrix(eigsum,mask)

#read in validation/testing set
rows.valid <- as.integer(test_set$mysubject_index)
mang.valid <- mang_mat[rows.valid, ]

#read eigenregions for best fold
#paste(output.path,extension,'sd2eig' ,as.character(i), 'fold', toString(fold), '.nii.gz',sep='')
eig_files <- list.files(path = paste(output.path,sep=''), pattern=paste('*', 'fold', toString(fold), '.nii.gz', sep=''),full.names = T,recursive = T)
best_eig_mat <- imagesToMatrix(eig_files,mask)
imgmat_mang_valid <- mang.valid %*% t(best_eig_mat)

# i<-1
# imgmat_mang_valid <- (mang.valid %*% t(matrixToImages(eanat_mang$eigenanatomyimages,mask = mask)[[i]])) # load those for fold 2
# e1<-eanat_mang$eigenanatomyimages
# imgmat_mang_valid <- (mang.valid %*% (e1[i,]))

dist4.valid  <- behavior[rows.valid, 'd4']
projs.valid <- data.frame(cbind(dist4.valid,imgmat_mang_valid))
colnames(projs.valid) <- c('Dist_4', paste0('Proj', c(1:ncol(imgmat_mang_valid )))) # insert column names

distpred <- predict.lm(mylm, newdata=projs.valid) 

mymodel<-lm(distpred~dist4.valid)

RSS <- c(crossprod(mymodel$residuals))
MSE <- RSS / length(mymodel$residuals)
RMSE <- sqrt(MSE)
RMSE2<-sqrt(mean((distpred - dist4.valid)^2)) 
mysummary <-summary(mymodel)
r2pred <- mysummary$adj.r.squared
ppred <- mysummary$coefficients[2,4]


max(behavior$d4[1:24])
RMSE_valid<-RMSE2
BIC_valid <- BIC(mymodel)
R2score_valid<-myr2score(distpred,dist4.valid)

myplot<- visreg(mymodel, gg=TRUE, ylim=c(0,1200),xlim=c(0,1200)) 
myplot + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
               panel.background = element_rect(fill = "transparent", colour = NA), 
               axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
               axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')) +
    #xlim(0,1200)+ylim(0,1200)+coord_cartesian(xlim = c(1200, 1200),ylim = c(1200, 1200)) + coord_equal()+
  ggtitle(paste("RMSE=",formatC(RMSE_valid,digits=2, format="f"), 
               # "R2score=",formatC(R2score_valid,digits=2, format="f"), 
               # " R2=", formatC(r2pred,digits=2, format="f"), 
                " p= ", formatC(ppred,digits=4, format="f"),
                "  BIC=", formatC(BIC_valid,digits=2, format="f")))

ggsave(paste(output.path,extension,'MnValidationSet',toString(fold),'sd2.pdf',sep=''), plot = last_plot(), device = 'pdf', 
       scale = 1, width = 4, height = 4, units = c("in"),dpi = 300)

numcols<-dim(imgmat_mang_valid)[2]
rd4 <- t(t(dist4.valid)[rep(1,c(numcols)),])
pcor<-c(numcols)
corval<-c(numcols)

cor(rd4,imgmat_mang_valid)
for (i in 1:numcols) {
  mypcor<-cor.test(t(dist4.valid),t(imgmat_mang_valid[,i]))
  pcor[i]<-mypcor$p.value
  corval[i]<-mypcor$estimate
}

# for (i in 1:numcols){
#   res1eig<-matrixToImages(jeanat_mang$fusedlist,mask = mask)[[i]] #eig1
#   antsImageWrite(res1eig,paste(output.path,extension,'JoinEanat' ,as.character(i), '.nii.gz',sep=''))
# }

mycorsdf_eig2d4<-data.frame(rbind(pcor,corval),row.names=c("pcor","cor"))
write.csv(mycorsdf_eig2d4, file = paste(output.path ,extension,'fold', toString(fold), 'd4corsvalidsd2.csv',sep=''))

myperf<-data.frame(rbind(distpred,dist4.valid),row.names=c("d_predicted","d_valid"))
write.csv(myperf, file = paste(output.path ,extension,'fold', toString(fold), 'distances4_validsd2.csv',sep=''))


#redo fold min or save models

