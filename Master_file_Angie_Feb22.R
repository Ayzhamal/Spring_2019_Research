
############################################################################
# Read in our training and testing data
library(data.table)
training <- read.csv('https://raw.githubusercontent.com/OliShawn/KmerResearch/master/Spring2019Kmer/12merFiles/train_sequences_12mer_stride1_tokens.csv.count.csv',header = TRUE)
#training <- read.csv('/Users/group14/Documents/Training_40.csv',header = TRUE)
training <- data.frame(training)
#training <- training[,names(training) != "DNA"] # use this line if you want to remove the 'name' of each row
training <- training[sample(nrow(training), nrow(training)), ] #randomizes the rows
training$class[training$class == "1"] <- "positive" # Assigns "Positive" to the 1 class
training$class[training$class == "0"] <- "negative" # Assigns 'Negative' to the 0 class
training$class <- factor(training$class)

#splitting the dataset 80% training, 20% testing and save them into csv files
set.seed(123)
indexes<-createDataPartition(training_dataset$class, times = 1, p=0.8, list = FALSE)
training<-training_dataset[indexes,]
testing<-training_dataset[-indexes,]
write.csv(training, "training_set.csv")
write.csv(testing, "testing_set.csv")


#Preparing testing data
testing = read.csv("https://raw.githubusercontent.com/OliShawn/KmerResearch/master/Spring2019Kmer/12merFiles/test_sequences_12mer_stride1_tokens.csv.count.csv", header = TRUE)
#testing = read.csv("/Users/group14/Documents/Testing_40.csv", header = TRUE)
testing <- data.frame(testing)
#testing <- testing[,names(testing) != "DNA"] # use this line if you want to remove the 'name' of each row
testing <- testing[sample(nrow(testing), nrow(testing)), ] #randomizes the rows
testing$class[testing$class == "1"] <- "positive" # Assigns "Positive" to the 1 class
testing$class[testing$class == "0"] <- "negative" # Assigns 'Negative' to the 0 class
testing$class <- factor(testing$class)
testing <- data.frame(testing)
#################################################################################

##############################################################################
# Libraries we need to use, majority use caret
suppressMessages(library(caret))
suppressMessages(library(e1071))


#############################################################################
#CARET Random Forest definition
do.RF <- function(training)
{  
  set.seed(313)
  n <- dim(training)[2]
  gridRF <- expand.grid(mtry = seq(from=0,by=as.integer(n/10),to=n)[-1]) #may need to change this depend on your data size
  ctrl.crossRF <- trainControl(method = "cv",number = 10,classProbs = TRUE,savePredictions = TRUE,allowParallel=TRUE)
  rf.Fit <- train(class ~ .,data = training,method = "rf",metric = "Accuracy",preProc = c("center", "scale"),
                  ntree = 200, tuneGrid = gridRF,trControl = ctrl.crossRF)
  rf.Fit
}

#training and testing
rf.Fit <- do.RF(training) #training done here
#print(rf.Fit)

Pred <-  predict(rf.Fit, testing)#prediction on the testing set
cm <- confusionMatrix(Pred, testing$class)
cm <- confusionMatrix(Pred, testing$class, mode = "prec_recall") 
result.roc <- plot(roc(testing$class, Pred$positive))
plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")

Pred <-  predict(rf.Fit, testing, type = "prob")
result.roc <- plot(roc(testing$class, Pred$positive))
plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")
auc(result.roc)

print("CM for RF:") 
print(cm)
saveRDS(rf.Fit, "RF.Rds") #saves the model to an rds file
#############################################################################


#############################################################################
#CARET boosted trees definition
do.Boost <- function(training)
{ 
  #trials = number of boosting iterations, or (simply number of trees)
  #winnow = remove unimportant predictors
  gridBoost <- expand.grid(model="tree",trials=seq(from=1,by=2,to=100),winnow=FALSE)
  set.seed(1)
  ctrl.crossBoost <- trainControl(method = "cv",number = 10,classProbs = TRUE,savePredictions = TRUE,allowParallel=TRUE)
  C5.0.Fit <- train(class ~ .,data = training,method = "C5.0",metric = "Accuracy",preProc = c("center", "scale"),
                    tuneGrid = gridBoost,trControl = ctrl.crossBoost)
  
  C5.0.Fit
}
#training
boost.Fit <- do.Boost(training)
#print(boost.Fit)

Pred <-  predict(boost.Fit,testing) #prediction on the testing set
cm <- confusionMatrix(Pred,testing$class)
print("CM for Boosted:")
print(cm)
saveRDS(boost.Fit, "Boost.rds") #saves the model to an rds file

#install.packages("pROC")
library(pROC)
Pred <-  predict(boost.Fit, testing, type = "prob")
result.roc <- plot(roc(testing$class, Pred$positive))
plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")

cm <- confusionMatrix(Pred, testing$class, mode = "prec_recall") 


RP.perf <- performance(Pred, "prec", "rec");
############################################################################


############################################################################
#CARET KNN 
#controls
#install.packages("knn", dependencies = TRUE)
grid = expand.grid(kmax=c(1:20),distance=10,kernel="optimal")
ctrl.cross <- trainControl(method="cv",number=10, classProbs=TRUE,savePredictions=TRUE)

#training
knnFit.cross <- train(class ~ .,
data = training, # training data
method ="kknn",  # model  
metric="Accuracy", #evaluation metric
preProc=c("center","scale"), # data to be scaled
tuneGrid = grid, # range of parameters to be tuned
trControl=ctrl.cross) # training controls
#print(knnFit.cross)
#plot(knnFit.cross)

#Testing
Pred <- predict(knnFit.cross,testing) #prediction on the testing set
cm<- confusionMatrix(Pred,testing$class)
cm <- confusionMatrix(Pred, testing$class, mode = "prec_recall") 
print("CM for KNN:")
print(cm)
saveRDS(knnFit.cross, "KNN.rds") #saves the model to an rds file

Pred <-  predict(knnFit.cross, testing, type = "prob")
result.roc <- plot(roc(testing$class, Pred$positive))
plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")
auc(result.roc)


#############################################################################


#############################################################################
#CARET Decision Tree:
#this is based on CARET, but sometimes doesn't run well, use the e1071 instead
do.DT <- function(training)
{
  set.seed(1)
  grid <- expand.grid(cp = 2^seq(from = -30 , to= 0, by = 2) )
  ctrl.cross <- trainControl(method = "cv", number = 5,classProbs = TRUE)
  dec_tree <-   train(class ~ ., data= training, perProc = c("center", "scale"),
                      method = 'rpart', #rpart for classif. dec tree
                      metric ='Accuracy',
                      tuneGrid= grid, trControl = ctrl.cross
  )
  dec_tree
}

dec.Fit <- do.DT(training)

Pred <- predict(dec.Fit,testing) #prediction on the testing set
cm<- confusionMatrix(Pred,testing$class) #
cm <- confusionMatrix(Pred, testing$class, mode = "prec_recall") # I need this cm
print("CM for DT:")
print(cm)
saveRDS(do.DT, "DT.rds") #saves the model to an rds file


Pred <-  predict(knnFit.cross, testing, type = "prob")
result.roc <- plot(roc(testing$class, Pred$positive))
plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")
auc(result.roc)
###############################################################################


###############################################################################
##### Our new decision tree. the e1071 package works better IMO################
suppressMessages(library(caret))
suppressMessages(library(e1071))

library(e1071)
do.e071.DT <- function(training)
{
  
  dec_tree <- tune.rpart(class ~ . , data = training , minsplit=c(5,10,15,20),cp = 2^seq(from = -20 , to= 0, by = 2) )
  dec_tree
}

DT.Fit <- do.e071.DT(training)
print(DT.Fit)
#predict using tuned DT.Fit
Pred <-  predict(DT.Fit$best.model,testing,type="class")
cm <- confusionMatrix(Pred,testing$class)
cm <- confusionMatrix(Pred, testing$class, mode = "prec_recall") 
print(cm)

Pred <-  predict(DT.Fit, testing, type = "prob")
result.roc <- plot(roc(testing$class, Pred[,2]))
plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")
auc(result.roc)


###############################################################################



################################################################################
#Regularization elastic-net logistic regression:
#install.packages("glmnet", repos = "http://cran.us.r-project.org", dependencies = TRUE)
library(glmnet)
traintest=rbind(training,testing)

X = sparse.model.matrix(as.formula(paste("class ~", paste(colnames(training[,-1]), sep = "", collapse=" +"))), data = traintest)

#Alpha = 0.1
logreg.fit01 = cv.glmnet(X[1:nrow(training),], training[,1], alpha = 0.1, family = "binomial",type.measure = "auc",nfolds = 10)
plot(logreg.fit01)
model$lambda.min
#predict on test set
pred01 = predict(model, s='lambda.min', newx=X[-(1:nrow(training)),], type="class")

#Alpha - 0.5
logreg.fit05 = cv.glmnet(X[1:nrow(training),], training[,1], alpha = 0.5, family = "binomial",type.measure = "auc",nfolds = 10)
plot(logreg.fit05)
model$lambda.min
#predict on test set
pred05 = predict(model, s='lambda.min', newx=X[-(1:nrow(training)),], type="class")

#Alpha = 0.8
logreg.fit08 = cv.glmnet(X[1:nrow(training),], training[,1], alpha = 0.8, family = "binomial",type.measure = "auc",nfolds = 10)
plot(logreg.fit08)
model$lambda.min

#predict on test set
pred08 = predict(model, s='lambda.min', newx=X[-(1:nrow(training)),], type = 'class')

#Confusion matrix:
cm01 <- confusionMatrix(factor(pred01),factor(testing$class), mode = "prec_recall")
cm05 <- confusionMatrix(factor(pred05),factor(testing$class))
cm08 <- confusionMatrix(factor(pred08),factor(testing$class))
print("CM for Log Reg, alpah = 0.1:")
print(cm01)
print("CM for Log Reg, alpah = 0.5:")
print(cm05)
print("CM for Log Reg, alpah = 0.8:")
print(cm08)
saveRDS(logreg.fit01, "/Users/group14/Documents/RDSfile/LogReg01_12mer.rds")
saveRDS(logreg.fit05, "/Users/group14/Documents/RDSfile//LogReg05_12mer.rds")
saveRDS(logreg.fit08, "/Users/group14/Documents/RDSfile//LogReg08_12mer.rds")

Pred <-  predict(DT.Fit$best.model,testing,type="class")
cm <- confusionMatrix(Pred,testing$class)
cm <- confusionMatrix(Pred, testing$class, mode = "prec_recall") 
print(cm)

# Pred <-  predict(DT.Fit$best.model, testing, type = "prob")
# result.roc <- plot(roc(testing$class, Pred[,2]))
# plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")
# auc(result.roc)
############################################################################


############################################################################





















