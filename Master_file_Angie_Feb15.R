
############################################################################
# Read in our training and testing data
library(data.table)
training = read.csv("/Users/group5/Documents/Ayzhamal/RStudio/training_set.csv", header = TRUE)
#training <- read.csv('/Users/group14/Documents/Training_40.csv',header = TRUE)
training <- data.frame(training)
#training <- training[,names(training) != "DNA"] # use this line if you want to remove the 'name' of each row
training <- training[sample(nrow(training), nrow(training)), ] #randomizes the rows
#training$class[training$class == "1"] <- "positive" # Assigns "Positive" to the 1 class
#training$class[training$class == "0"] <- "negative" # Assigns 'Negative' to the 0 class
training$class <- factor(training$class)



#Preparing testing data
testing = read.csv("/Users/group5/Documents/Ayzhamal/RStudio/testing_set.csv", header = TRUE)
#testing = read.csv("/Users/group14/Documents/Testing_40.csv", header = TRUE)
testing <- data.frame(testing)
#testing <- testing[,names(testing) != "DNA"] # use this line if you want to remove the 'name' of each row
testing <- testing[sample(nrow(testing), nrow(testing)), ] #randomizes the rows
#testing$class[testing$class == "1"] <- "positive" # Assigns "Positive" to the 1 class
#testing$class[testing$class == "0"] <- "negative" # Assigns 'Negative' to the 0 class
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
  gridRF <- expand.grid(mtry = seq(from=0,by=as.integer(n/40),to=n)[-1]) #may need to change this depend on your data size
  ctrl.crossRF <- trainControl(method = "cv",number = 5,classProbs = TRUE,savePredictions = TRUE,allowParallel=TRUE)
  rf.Fit <- train(class ~ .,data = training,method = "rf",metric = "Accuracy",preProc = c("center", "scale"),
                  ntree = 200, tuneGrid = gridRF,trControl = ctrl.crossRF)
  rf.Fit
}

#training and testing
rf.Fit <- do.RF(training) #training done here
#print(rf.Fit)

Pred <-  predict(rf.Fit, testing)#prediction on the testing set
cm <- confusionMatrix(Pred, testing$class)
print("CM for RF:") 
print(cm)
saveRDS(rf.Fit, "/Users/group13/Downloads/kmer_6_12/RF_6mer.rds") #saves the model to an rds file
#############################################################################


#############################################################################
#CARET boosted trees definition
do.Boost <- function(training)
{ 
  #trials = number of boosting iterations, or (simply number of trees)
  #winnow = remove unimportant predictors
  gridBoost <- expand.grid(model="tree",trials=seq(from=1,by=2,to=100),winnow=FALSE)
  set.seed(1)
  ctrl.crossBoost <- trainControl(method = "cv",number = 5,classProbs = TRUE,savePredictions = TRUE,allowParallel=TRUE)
  C5.0.Fit <- train(class ~ .,data = training,method = "C5.0",metric = "Accuracy",preProc = c("center", "scale"),
                    tuneGrid = gridBoost,trControl = ctrl.crossBoost)
  
  C5.0.Fit
}
#training
boost.Fit <- do.Boost(training)
#print(boost.Fit)

Pred <-  predict(boost.Fit,testing) #prediction on the testing set
cm <- confusionMatrix(Pred,testing$class)
cm <- confusionMatrix(Pred,testing$class, mode = "prec_recall")
print("CM for Boosted:")
print(cm)
saveRDS(boost.Fit, "/Users/group13/Downloads/kmer_6_12/Boost_6mer.rds") #saves the model to an rds file
############################################################################


############################################################################
#CARET KNN 
#controls
#install.packages("knn", dependencies = TRUE)
grid = expand.grid(kmax=c(1:20),distance=10,kernel="optimal")
ctrl.cross <- trainControl(method="cv",number=5, classProbs=TRUE,savePredictions=TRUE)

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
Pred <- predict(KNN_6mer,testing)
cm<- confusionMatrix(Pred,testing$class)
cm<- confusionMatrix(Pred,testing$class, mode = 'prec_recall')
print("CM for KNN:")
print(cm)
saveRDS(knnFit.cross, "/Users/group13/Downloads/kmer_6_12/KNN_6mer.rds") #saves the model to an rds file
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

DT.Fit <- do.DT(training)

Pred <- predict(DT_Fit,testing) #prediction on the testing set
Pred <- predict(RF_6mer,testing)
cm<- confusionMatrix(Pred,testing$class) 
cm<- confusionMatrix(Pred,testing$class, mode = 'prec_recall')
print("CM for DT:")
print(cm)
saveRDS(do.DT, "/Users/group13/Downloads/kmer_6_12/DT_6mer.rds") #saves the model to an rds file
###############################################################################



################################################################################
#Regularization elastic-net logistic regression:
#install.packages("glmnet", repos = "http://cran.us.r-project.org", dependencies = TRUE)
#install.packages('ROCR')
#install.packages('pROC')
library(glmnet)
library(pROC)
library(ROCR)
traintest=rbind(training,testing)

X = sparse.model.matrix(as.formula(paste("class ~", paste(colnames(training[,-1]), sep = "", collapse=" +"))), data = traintest)

#Alpha = 0.1
logreg.fit01 = cv.glmnet(X[1:nrow(training),], training[,1], alpha = 0.1, family = "binomial",type.measure = "auc",nfolds = 10)
plot(logreg.fit01)
logreg.fit01$lambda.min
#predict on test set
pred01 = predict(logreg.fit01, s='lambda.min', newx=X[-(1:nrow(training)),], type="response")
perf01 <- performance(prediction(pred01, testing$class), 'auc')
perf01@y.values

#Alpha - 0.5
logreg.fit05 = cv.glmnet(X[1:nrow(training),], training[,1], alpha = 0.5, family = "binomial",type.measure = "auc",nfolds = 10)
plot(logreg.fit05)
logreg.fit05$lambda.min
#predict on test set
pred05 = predict(logreg.fit05, s='lambda.min', newx=X[-(1:nrow(training)),], type="response")
perf05 <- performance(prediction(pred05, testing$class), 'auc')
perf05@y.values

#Alpha = 0.8
logreg.fit08 = cv.glmnet(X[1:nrow(training),], training[,1], alpha = 0.8, family = "binomial",type.measure = "auc",nfolds = 10)
plot(logreg.fit08)
logreg.fit08$lambda.min
#predict on test set
pred08 = predict(logreg.fit08, s='lambda.min', newx=X[-(1:nrow(training)),], type = 'response')
perf08 <- performance(prediction(pred08, testing$class), 'auc')
perf08@y.values

#Confusion matrix:
pred01 = predict(logreg.fit01, s='lambda.min', newx=X[-(1:nrow(training)),], type="class")
#cm01 <- confusionMatrix(factor(pred01),factor(testing$class))
cm01 <- confusionMatrix(factor(pred01),factor(testing$class), mode = 'prec_recall')
pred05 = predict(logreg.fit05, s='lambda.min', newx=X[-(1:nrow(training)),], type="class")
cm05 <- confusionMatrix(factor(pred05),factor(testing$class))
#cm05 <- confusionMatrix(factor(pred05),factor(testing$class), mode = 'prec_recall')
pred08 = predict(logreg.fit08, s='lambda.min', newx=X[-(1:nrow(training)),], type="class")
cm08 <- confusionMatrix(factor(pred08),factor(testing$class))
#cm08 <- confusionMatrix(factor(pred08),factor(testing$class), mode = 'prec_recall')

print("CM for Log Reg, alpah = 0.1:")
print(cm01)
print("CM for Log Reg, alpah = 0.5:")
print(cm05)
print("CM for Log Reg, alpah = 0.8:")
print(cm08)
saveRDS(logreg.fit01, "/Users/group5/Documents/Ayzhamal/RStudio/LogReg01.rds")
saveRDS(logreg.fit05, "/Users/group5/Documents/Ayzhamal/RStudio/LogReg05.rds")
saveRDS(logreg.fit08, "/Users/group5/Documents/Ayzhamal/RStudio/LogReg08.rds")

#area under the curve 
Pred05 <-  predict(logreg.fit05, testing, type = "class")
result.roc <- plot(roc(testing$class, sort(pred05[1])))
plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")
auc(logreg.fit05)
############################################################################


###########################################################################
#This section plots the ROC of a model and gives the AUC
library('pROC')
Pred <-  predict(RF_6mer, testing, type = "prob")
result.roc <- plot(roc(testing$class, Pred$positive))
plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")
auc(result.roc)
##########################################################################
#This decision tree model works better
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

#area under the curve 
Pred <-  predict(logreg.fit01, testing, type = "prob")
result.roc <- plot(roc(testing$class, sort(pred01[1])))
plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")
auc(logreg.fit01)

################################################################################



