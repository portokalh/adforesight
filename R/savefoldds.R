k<-4
set.seed(1)
res_train<-createFolds(behavior$genotype,k, list = TRUE, returnTrain = TRUE)
set.seed(1)
res_test<-createFolds(behavior$genotype,k)

f1train<-mydfb[res_train$Fold1,]
f1test<-mydfb[res_test$Fold1,]
write.csv(f1test,file = "f1test.csv")
write.csv(f1train,file = "f1train.csv")

f2train<-mydfb[res_train$Fold2,]
f2test<-mydfb[res_test$Fold2,]
write.csv(f2test,file = "f2test.csv")
write.csv(f2train,file = "f2train.csv")

f3train<-mydfb[res_train$Fold3,]
f3test<-mydfb[res_test$Fold3,]
write.csv(f3test,file = "f3test.csv")
write.csv(f3train,file = "f3train.csv")

f4train<-mydfb[res_train$Fold4,]
f4test<-mydfb[res_test$Fold4,]
write.csv(f4test,file = "f4test.csv")
write.csv(f4train,file = "f4train.csv")
