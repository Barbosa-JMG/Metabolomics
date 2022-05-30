##Code date: 2022 - 05 - 17
##Name: João Marcos Gonçalves Barbosa         
##ORCID address: https://orcid.org/0000-0003-3560-5838
##GitHub adress: https://github.com/Barbosa-JMG
##Obs:RScript applied to the paper "A cerumenolomic approach to bovine trypanosomosis diagnosis"
##In case of use, Please cite: 

                  ###############################################################################

                                          ####Model I####
##Data visualization of the Model I results 
##Use the Model_1.csv file contained the binary data for the 20 Healthy/Inoculated biomarkers (Previous selected by GAPLS)

                  ### HCA Ward Hamming 
df <- read.csv("Model_1.csv", sep=",",header = TRUE) 
Num <- df[,3:22]
Num.mat <- as.matrix(Num)
Step <- as.factor(df$Step)

library(e1071)
groupCodes <- as.character(df$Step)
dist.ham<-as.dist(hamming.distance(Num.mat))
rownames(Num.mat) <- make.unique(groupCodes)
colorCodes <- c(A="blue", PI="red")
hc <- hclust(dist.ham, method = "ward.D")
d <- as.dendrogram(hc)
require(dendextend)
labels_colors(d) <- colorCodes[groupCodes][order.dendrogram(d)]
par(cex=0.9)
plot(d, main = " Model 1; ward; dist = Ham")


###Validation for the 20 VOMs selected by GA-PLS (Model I, Chromosome 1)
library(caret)
library(pls)
df <- read.csv("Model_1.csv", sep=",",header = TRUE) 
Num <- df[,2:22]
Num$Step <- factor(Num$Step)
Num.mat <- as.matrix(Num)


set.seed(1)
ctrl <- trainControl(method = "repeatedcv",
                          number = 10,
                          repeats = 10,
                          classProbs = TRUE,
                          summaryFunction = twoClassSummary)

set.seed(7)
plsFit <- train(
          Step ~ .,
          data = Num.mat,
          method = "pls",
          trControl = ctrl,
          metric = "ROC",
          tuneLength = 5,
          tuneGrid = expand.grid(ncomp = 1:5))

plsFit
varImp(plsFit)
ggplot(plsFit)+ theme_bw()+theme(text=element_text(size = 16, family="serif"))

ggsave("Model1_20GA.tiff")



### GA-PLS (Model I, Chromosome 2)
#Use the BinaryData_Model1.csv file containing the 70VOMs detected in the study 
df <- read.csv("BinnaryData_Model1.csv", sep=",",header = TRUE) 
Step <- as.factor(df$Step)
Num <- df[,3:72]
Num.mat <- as.matrix(Num)
Mx <- df[,2:72]

library("plsVarSel")
library('caret')
set.seed(456)
X <- Num.mat
y <- Step
xga <- ga_pls(y, X, GA.threshold = 10, iters = 100, popSize = 100)
Xga <- Mx[,xga$ga.selection]  
Xga$Step <- Mx$Step

set.seed(8)
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 10,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

set.seed(16)
plsFit <- train(
  Step ~ .,
  data = Xga,
  method = "pls",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  tuneGrid = expand.grid(ncomp = 1:5))

plsFit
varImp(plsFit)
ggplot(plsFit)+ theme_bw()+theme(text=element_text(size = 16, family="serif"))

ggsave("Model1_SecondGA.tiff")

##GA-PLS (Model I, Chromosome 3)

set.seed(356)
X <- Num.mat
y <- Step
xga <- ga_pls(y, X, GA.threshold = 10, iters = 100, popSize = 100)
Xga <- Mx[,xga$ga.selection]  
Xga$Step <- Mx$Step

set.seed(32)
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 10,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

set.seed(46)
plsFit <- train(
  Step ~ .,
  data = Xga,
  method = "pls",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  tuneGrid = expand.grid(ncomp = 1:5))

plsFit
varImp(plsFit)
ggplot(plsFit)+ theme_bw()+theme(text=element_text(size = 16, family="serif"))

ggsave("Model1_ThirdGA.tiff")

##GA-PLS (Model I, Chromosome 4)
set.seed(24)
X <- Num.mat
y <- Step
xga <- ga_pls(y, X, GA.threshold = 10, iters = 100, popSize = 100)
Xga <- Mx[,xga$ga.selection]  
Xga$Step <- Mx$Step

set.seed(790)
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 10,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

set.seed(820)
plsFit <- train(
  Step ~ .,
  data = Xga,
  method = "pls",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  tuneGrid = expand.grid(ncomp = 1:5))

plsFit
varImp(plsFit)
ggplot(plsFit)+ theme_bw()+theme(text=element_text(size = 16, family="serif"))

ggsave("Model1_FourthGA.tiff")

##GA-PLS (Model I, Chromosome 5)

set.seed(1036)
X <- Num.mat
y <- Step
xga <- ga_pls(y, X, GA.threshold = 10, iters = 100, popSize = 100)
Xga <- Mx[,xga$ga.selection]  
Xga$Step <- Mx$Step

set.seed(2065)
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 10,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

set.seed(1820)
plsFit <- train(
  Step ~ .,
  data = Xga,
  method = "pls",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  tuneGrid = expand.grid(ncomp = 1:5))

plsFit
varImp(plsFit)
ggplot(plsFit)+ theme_bw()+theme(text=element_text(size = 16, family="serif"))

ggsave("Model1_FifthGA.tiff")


##GA-PLS (Model I, Chromosome 6)

set.seed(1236)
X <- Num.mat
y <- Step
xga <- ga_pls(y, X, GA.threshold = 10, iters = 100, popSize = 100)
Xga <- Mx[,xga$ga.selection]  
Xga$Step <- Mx$Step

set.seed(2265)
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 10,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

set.seed(1880)
plsFit <- train(
  Step ~ .,
  data = Xga,
  method = "pls",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  tuneGrid = expand.grid(ncomp = 1:5))

plsFit
varImp(plsFit)
ggplot(plsFit)+ theme_bw()+theme(text=element_text(size = 16, family="serif"))

ggsave("Model1_SixthGA.tiff")



                                        ####Model II####
#Data visualization of the Model II results 
##Use the Model_2.csv file contained the binary data for the 13 Chronic\acute phase biomarkers (Previous selected by GAPLS)

                              ##HCA Ward Hamming 
df <- read.csv("Model_2.csv", sep=",",header = TRUE) 
Num <- df[,3:15]
Num.mat <- as.matrix(Num)
Phase<- as.factor(df$Phase)

library(e1071)
groupCodes <- as.character(df$Phase)
dist.ham<-as.dist(hamming.distance(Num.mat))
rownames(Num.mat) <- make.unique(groupCodes)
colorCodes <- c(PIA="purple", PIC="green")
hc <- hclust(dist.ham, method = "ward.D")
d <- as.dendrogram(hc)
require(dendextend)
labels_colors(d) <- colorCodes[groupCodes][order.dendrogram(d)]
par(cex=0.9)
plot(d, main = " Model 2; ward; dist = Ham")

######Validation for the 13 VOMs selected by GA-PLS (Model II, Chromosome 1)
library(caret)
library(pls)
df <- read.csv("Model_2.csv", sep=",",header = TRUE) 
Num <- df[,2:15]
Num$Phase <- factor(Num$Phase)
Num.mat <- as.matrix(Num)


set.seed(59)
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 10,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

set.seed(771)
plsFit <- train(
  Phase ~ .,
  data = Num.mat,
  method = "pls",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  tuneGrid = expand.grid(ncomp = 1:5))

plsFit
varImp(plsFit)
ggplot(plsFit)+ theme_bw()+theme(text=element_text(size = 16, family="serif"))

ggsave("Model2_FirstGA.tiff")

##GA-PLS (Model II, Chromosome 2)
df <- read.csv("BinnaryData_Model2.csv", sep=",",header = TRUE) 
Phase <- as.factor(df$Phase)
Num <- df[,3:72]
Num.mat <- as.matrix(Num)
Mx <- df[,2:72]


library("plsVarSel")
library('caret')
set.seed(4184)
X <- Num.mat
y <- Phase
xga <- ga_pls(y, X, GA.threshold = 10, iters = 100, popSize = 100)
Xga <- Mx[,xga$ga.selection]  
Xga$Phase <- Mx$Phase

set.seed(8889)
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 10,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

set.seed(106)
plsFit <- train(
  Phase ~ .,
  data = Xga,
  method = "pls",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  tuneGrid = expand.grid(ncomp = 1:5))

plsFit
varImp(plsFit)
ggplot(plsFit)+ theme_bw()+theme(text=element_text(size = 16, family="serif"))

ggsave("Model2_SecondGA.tiff")

##GA-PLS (Model II, Chromosome 3)

set.seed(5555)
X <- Num.mat
y <- Phase
xga <- ga_pls(y, X, GA.threshold = 10, iters = 100, popSize = 100)
Xga <- Mx[,xga$ga.selection]  
Xga$Phase <- Mx$Phase

set.seed(66798)
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 10,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

set.seed(11106)
plsFit <- train(
  Phase ~ .,
  data = Xga,
  method = "pls",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  tuneGrid = expand.grid(ncomp = 1:5))

plsFit
varImp(plsFit)
ggplot(plsFit)+ theme_bw()+theme(text=element_text(size = 16, family="serif"))

ggsave("Model2_ThirdGA.tiff")

##GA-PLS (Model II, Chromosome 4)

set.seed(5552)
X <- Num.mat
y <- Phase
xga <- ga_pls(y, X, GA.threshold = 10, iters = 100, popSize = 100)
Xga <- Mx[,xga$ga.selection]  
Xga$Phase <- Mx$Phase

set.seed(600798)
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 10,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

set.seed(11966)
plsFit <- train(
  Phase ~ .,
  data = Xga,
  method = "pls",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  tuneGrid = expand.grid(ncomp = 1:5))

plsFit
varImp(plsFit)
ggplot(plsFit)+ theme_bw()+theme(text=element_text(size = 16, family="serif"))

ggsave("Model2_FourthGA.tiff")


##GA-PLS (Model II, Chromosome 5)

set.seed(24)
X <- Num.mat
y <- Phase
xga <- ga_pls(y, X, GA.threshold = 10, iters = 100, popSize = 100)
Xga <- Mx[,xga$ga.selection]  
Xga$Phase <- Mx$Phase

set.seed(777798)
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 10,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

set.seed(98966)
plsFit <- train(
  Phase ~ .,
  data = Xga,
  method = "pls",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  tuneGrid = expand.grid(ncomp = 1:5))

plsFit
varImp(plsFit)
ggplot(plsFit)+ theme_bw()+theme(text=element_text(size = 16, family="serif"))

ggsave("Model2_FifthGA.tiff")

##GA-PLS (Model II, Chromosome 6)

set.seed(42)
X <- Num.mat
y <- Phase
xga <- ga_pls(y, X, GA.threshold = 10, iters = 100, popSize = 100)
Xga <- Mx[,xga$ga.selection]  
Xga$Phase <- Mx$Phase

set.seed(70215)
ctrl <- trainControl(method = "repeatedcv",
                     number = 10,
                     repeats = 10,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

set.seed(10)
plsFit <- train(
  Phase ~ .,
  data = Xga,
  method = "pls",
  trControl = ctrl,
  metric = "ROC",
  tuneLength = 5,
  tuneGrid = expand.grid(ncomp = 1:5))

plsFit
varImp(plsFit)
ggplot(plsFit)+ theme_bw()+theme(text=element_text(size = 16, family="serif"))

ggsave("Model2_SixthGA.tiff")
