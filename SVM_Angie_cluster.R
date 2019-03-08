
############################################################################
############################################################################
# Read in our training and testing data
library(data.table)
training <- read.csv("https://raw.githubusercontent.com/Ayzhamal/Spring_2019_Research/master/datasets/rao_imr90%20dataset/training_set.csv",header = TRUE)
training <- data.frame(training)
#training <- training[,names(training) != "DNA"] # use this line if you want to remove the 'name' of each row
training <- training[sample(nrow(training), nrow(training)), ] #randomizes the rows
#training$class[training$class == "real"] <- "positive" # Assigns "Positive" to the real class
#training$class[training$class == "random"] <- "negative" # Assigns 'Negative' to the random class
training$class <- factor(training$class)


#Preparing testing data
testing = read.csv("https://raw.githubusercontent.com/Ayzhamal/Spring_2019_Research/master/datasets/rao_imr90%20dataset/testing_set.csv", header = TRUE)
testing <- data.frame(testing)
#testing <- testing[,names(testing) != "DNA"] # use this line if you want to remove the 'name' of each row
testing <- testing[sample(nrow(testing), nrow(testing)), ] #randomizes the rows
#testing$class[testing$class == "real"] <- "positive" # Assigns "Positive" to the real class
#testing$class[testing$class == "random"] <- "negative" # Assigns 'Negative' to the random class
testing$class <- factor(testing$class)
testing <- data.frame(testing)
#################################################################################

##############################################################################
# Libraries we need to use, majority use caret
suppressMessages(library(caret))
suppressMessages(library(e1071))
library(doParallel) #parallel processing 
registerDoParallel(10)
r=getOption("repos")
r["CRAN"]="http://cran.us.r-project.org"
options(repos=r)

setwd("/home/angiez1/rao_imr90")
############################################################################
#install.packages("kernlab")
#install.packages("pROC")
library('pROC')
suppressMessages(library(kernlab))
do.RadialKernelSVM <- function(training)
{
  set.seed(1)
  tmpTraining <- training
  tmpTraining$class <- NULL
  sigma=sigest(as.matrix(tmpTraining)) # sigest returns 3 values of sigma 
  grid <- expand.grid(sigma = sigma , C = 2^seq(from=-4,by = 1, to =8)) # set up sigma and cost parameters
  ctrl.cross <- trainControl(method = "cv", number = 5,classProbs = TRUE,savePredictions=TRUE)
  svm.Fit <- caret::train(class ~ ., data= training,perProc = c("center"),
                          method = 'svmRadial', 
                          metric ='Accuracy',
                          tuneGrid= grid,
                          trControl = ctrl.cross
  )
  svm.Fit
}
svm.radial.Fit<-do.RadialKernelSVM(training)
print(svm.radial.Fit)

Pred <- predict(svm.radial.Fit,testing)#prediction on the testing set
cm<- confusionMatrix(Pred,testing$class)
print("CM for SVM: Sensitivity...")
print(cm)

cm <- confusionMatrix(Pred, testing$class, mode = "prec_recall") # will show Precision, Recall, F1 Score
print("CM for SVM: Recall...")
print(cm)

saveRDS(svm.radial.Fit, "Svm.radial.rds") #saves the model to an rds file

Pred <-  predict(svm.radial.Fit, testing, type = "prob")
result.roc <- plot(roc(testing$class, Pred$random))
#plot(result.roc, print.thres="best", print.thres.best.method="closest.topleft")
auc(result.roc)
#############################################################################









