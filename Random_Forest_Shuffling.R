
############################################################################
############################################################################
# Read in our testing data
library(data.table)

#Preparing testing data
testing = read.csv("https://raw.githubusercontent.com/Ayzhamal/Spring_2019_Research/master/datasets/rao_hmec_10kb%20dataset/testing_set.csv", header = TRUE)
testing <- data.frame(testing)
#testing <- testing[,names(testing) != "DNA"] # use this line if you want to remove the 'name' of each row
testing <- testing[sample(nrow(testing), nrow(testing)), ] #randomizes the rows
#testing$class[testing$class == "real"] <- "positive" # Assigns "Positive" to the real class
#testing$class[testing$class == "random"] <- "negative" # Assigns 'Negative' to the random class
testing$class <- factor(testing$class)
#################################################################################

##############################################################################
# Libraries we need to use, majority use caret
suppressMessages(library(caret))
suppressMessages(library(e1071))
library(doParallel) #parallel processing 
registerDoParallel(10)
#r=getOption("repos")
#r["CRAN"]="http://cran.us.r-project.org"
#options(repos=r)

setwd("/home/angiez1/output_files")
############################################################################


#############################################################################
# 1.CARET Random Forest definition
# do.RF <- function(training)
# {  
#   set.seed(313)
#   n <- dim(training)[2]
#   gridRF <- expand.grid(mtry = seq(from=0,by=as.integer(n/40),to=n)[-1]) #may need to change this depend on your data size
#   ctrl.crossRF <- trainControl(method = "cv",number = 5,classProbs = TRUE,savePredictions = TRUE,allowParallel=TRUE)
#   rf.Fit <- train(class ~ .,data = training,method = "rf",metric = "Accuracy",preProc = c("center", "scale"),
#                   ntree = 200, tuneGrid = gridRF,trControl = ctrl.crossRF)
#   rf.Fit
# }
# 
# #training and testing
# rf.Fit <- do.RF(training) #training done here
# rf.Fit

# Pred <-  predict(rf.Fit, testing) #prediction on the testing set
# cm <- confusionMatrix(Pred, testing$class) # will show accuracy, sensitivity, specifity
# print("CM for RF: Sensitivity...:") 
# print(cm)
# 
# cm <- confusionMatrix(Pred, testing$class, mode = "prec_recall") # will show Precision, Recall, F1 Score
# print("CM for RF: Recall...") 
# print(cm)

#install.packages("pROC")
# library('pROC')
# Pred<-predict(rf.Fit, testing, type="prob")
# result.roc <- plot(roc(testing$class, Pred$random))
# plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")
# auc(result.roc)

#install.packages('e1071', dependencies=TRUE)

#read in an Rds file
rf.Fit<-readRDS("/home/angiez1/rds_files/RF.Rds")
Prediction <-  predict(rf.Fit, testing) #prediction on the testing set
cm <- confusionMatrix(Prediction, testing$class) # will show accuracy, sensitivity, specifity
initialAccuracy<-cm$overall['Accuracy']
initialAccuracy

#Shuffle predictions for variable importance
AccShuffle<-NULL
shuffletimes<-100

outcomeName<-'class'
predictorNames<-setdiff(names(testing), outcomeName)
featureMeanAccs<-c()
for (feature in predictorNames){
  #print(feature)
  featureAccs<-c()
  shuffledData<-testing[,predictorNames]
  for (iter in 1:shuffletimes) {
    shuffledData[,feature]<-sample(shuffledData[,feature], length(shuffledData[,feature]))
    
    Pred<-predict(rf.Fit, shuffledData[, predictorNames])
    #print(dim(shuffledData[,predictorNames]))
    cm <- confusionMatrix(Pred, testing$class) # will show accuracy, sensitivity, specifity
    featureAccs<-c(featureAccs, cm$overall['Accuracy'])
  }
  meanAccs<-mean(featureAccs)
  difference<-initialAccuracy-meanAccs
  featureMeanAccs<-c(featureMeanAccs, difference)
}
#options(digits=5)
AccShuffle<-data.frame('feature'=predictorNames, 'accuracyDifference'=featureMeanAccs)
AccShuffle<-AccShuffle[order(AccShuffle$accuracyDifference, decreasing = TRUE),]
print(AccShuffle)

#save a csv file with results
write.csv(AccShuffle, file = "ShuffleResult_RF.csv")

