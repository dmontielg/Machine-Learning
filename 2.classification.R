######################################################
# Pattern Recognition Project (Week 2)
# Name(s)           Student Number
# Anne Tjallingii:  911024836020
# Diego Montiel:    880505580110
######################################################
# install.packages("tree")
# install.packages("randomForest")
# install.packages("gbm")
# install.packages("e1071")

library(gbm)
library(tree)
library(ISLR)
library(MASS)
library(class)
library(e1071)
library(randomForest)

data = read.table("expr4T.dat", header = TRUE)
d = data

data = d

# Last column contains tissue label
nn = ncol(data)-1
m = sapply(data[,1:nn],mean)
s = sapply(data[,1:nn],sd)
sm = s/m

#print the graphs of the m and sm
par(mfrow = c(2,2));plot(m);plot(sm);plot(m[m < 15]);plot(sm[sm < 1]);

minsm= 2.5
minm= 20

# In the following, we make sure that tissue label is kept;
# by giving it an artificial sufficiently high value for m and sm
m = c(m,minm+1)
sm = c(sm,minsm+1)

dim(data)
length(which(sm > minsm & m > minm))
data = data[,which(sm > minsm & m > minm)]
dim(data)

#Check if these two genes are in the new dataset
data[,"ENSG00000271043.1_MTRNR2L2"]
data[,"ENSG00000229344.1_RP5.857K21.7"]

#For getting the names of the tissues
unique(data[dim(data)[2]])

tissue1 = "brain_amygdala"
tissue2 = "brain_hippocampus"
data = data.frame(data)
mydat = data[which(data$tissue == tissue1|data$tissue==tissue2),]
mydat = droplevels(mydat)
mydat
#######################
# Classification tree
#######################

par(mfrow =c(1,1))

set.seed(2)
n = round(dim(mydat)[1]*.5)

train1 = sample(1:nrow(mydat), n)
mydat.test = mydat[-train1,]
tissue.test = mydat$tissue[-train1]
#Applying tree
tree.dataset1 = tree(tissue ~., mydat, subset = train1)
tree.pred = predict(tree.dataset1, mydat.test, type = "class")
tree.table = table(tree.pred, tissue.test)
tree.table

((tree.table[1] + tree.table[4]) / sum(tree.table))

plot(tree.dataset1)
text(tree.dataset1, pretty = 0)
tree.cv2 = cv.tree(tree.dataset1)
plot(tree.cv2$size, tree.cv2$dev, type = "b")
#optimal is 3 so prune back

prune.genes = prune.tree(tree.dataset1, best = 3)
plot(prune.genes)
text(prune.genes, pretty=0)

prediction.prune = predict(prune.genes, mydat.test, type = "class")
tree.table.prune = table(prediction.prune, tissue.test)
tree.table.prune
((tree.table[1] + tree.table[4]) / sum(tree.table))

#######################
# Regression tree
#######################
#fix(data)

n = round(dim(data)[1]*.5)

train2 = sample(1:nrow(data),n)
test2 = data[-train2,]

tree.dataset2 = tree(ENSG00000271043.1_MTRNR2L2 ~.-tissue, data = data, subset = train2)
par(mfrow = c(1,1))

summary(tree.dataset2)
plot(tree.dataset2); text(tree.dataset2, pretty = 0)
tree.cv = cv.tree(tree.dataset2)

plot(tree.cv$size, tree.cv$dev, type = "b")

#Prunning it doesnt make sense because the biggest value from the plot before
#is already taken in it = 7

yhat = predict(tree.dataset2, newdata = test2)
plot(yhat, test2$ENSG00000271043.1_MTRNR2L2); abline(0,1)
mean((yhat-test2$ENSG00000271043.1_MTRNR2L2)^2)

######################################
# RANDOM FOREST Regression
#####################################

install.packages("randomForest")
library(randomForest)
set.seed(1)
#The argument mtry=13 indicates that all 13 predictors 
#should be considered for each split of the tree-in other words,
#that bagging should be done
train3 = sample(1:nrow(data), nrow(data)/2)
test3 = data[-train3,]

bag = randomForest(ENSG00000271043.1_MTRNR2L2 ~.-tissue, data = data, subset = train3,
                   importance = TRUE, ntree=1000, mtry = 19)
importance(bag)
yhat.bag = predict(bag, newdata = test3)
mean((yhat.bag - test3$ENSG00000271043.1_MTRNR2L2)^2)
plot(yhat.bag, test3$ENSG00000271043.1_MTRNR2L2); abline(0,1)
plot(bag)

######################################
# RANDOM FOREST CLASSIFICATION
#####################################

set.seed(2)
n2 = round(dim(mydat)[1]*.5)
train.random.forest2 = sample(1:nrow(mydat), n2)
test.random.forest2 = mydat[-train.random.forest2,]
tissue.test2 = mydat$tissue[-train.random.forest2]
tissue.test2

bag.classification = randomForest(tissue~., data = mydat, subset = train.random.forest2,
                                  mtry = 19, importance = TRUE, ntree = 25)


yhat.bag2 = predict(bag.classification, test.random.forest2, type = "class")
yhat.bag2
tree.table2 = table(yhat.bag2, tissue.test2)
tree.table2
plot(yhat.bag2)
((tree.table2[1] + tree.table2[4]) / sum(tree.table2))

##############################################
# SVM RADIAL
##############################################

#install.packages("e1071")
library(e1071)
n = round(dim(data)[1]*.5)
svm.train = sample(1:nrow(data),n)
svm.test = data[-svm.train,]

svm.fit = svm(ENSG00000271043.1_MTRNR2L2 ~.-tissue, data = data[svm.train,], kernel ="radial")
summary(svm.fit)
svm.fit
plot(svm.fit, data[svm.train, ])

prediction.svm = predict(svm.fit, newdata = svm.test)

plot(prediction.svm, svm.test$ENSG00000271043.1_MTRNR2L2); abline(0,1)
mean((prediction.svm-svm.test$ENSG00000271043.1_MTRNR2L2)^2)

############################################
#SVM POLYNOMIAL
###########################################

svm.fit2 = svm(ENSG00000271043.1_MTRNR2L2 ~ . -tissue, data = mydat[svm.train,], kernel = "polynomial")

prediction.svm2 = predict(svm.fit2, newdata = svm.test)

plot(prediction.svm2, svm.test$ENSG00000271043.1_MTRNR2L2); abline(0,1)
mean((prediction.svm2-svm.test$ENSG00000271043.1_MTRNR2L2)^2)

#############################################
# Predicting using boosting (CLASSIFICATION)
#############################################
#install.packages("gbm")
library(gbm)

set.seed(1)
n = round(dim(mydat)[1]*.5)
train = sample(1:nrow(mydat), n)
test = mydat[-train,]
tissue.test = mydat$tissue[-train]
# tissue1 = "brain_amygdala"
# tissue2 = "brain_hippocampus"
mydat$tissue01 = ifelse(mydat$tissue == "brain_amygdala", 1, 0)

boost.best = gbm(tissue01 ~ .-tissue, data = mydat[train,], distribution = "bernoulli", 
                 n.trees = 500, interaction.depth = 4)

summary(boost.best)
mydat$tissue01 <- NULL
#plot(summary(boost.best))
par(mfrow =c(1,1))
#Plot the most relevant genes
sum2<-summary(boost.best); 
plot(sum2[1:100,]$rel.inf, xlab = "Most relevant genes")
yhat.boost = predict(boost.best, 
                     newdata = mydat[-train,], 
                     n.trees = 500)
yhat.boost
#Most relevant genes
# ENSG00000118271.5_TTR
# ENSG00000133048.8_CHI3L1
# ENSG00000183395.4_PMCH

#############################################
#RETRAINING CLASSIFICATION TREE WITH 1 GENES
#############################################

par(mfrow =c(1,1))

set.seed(2)
n = round(dim(mydat)[1]*.5)
train = sample(1:nrow(mydat), n)
test = mydat[-train,]
tissue.test = mydat$tissue[-train]
# fix(mydat)
# relevant genes
# ENSG00000118271.5_TTR
# ENSG00000133048.8_CHI3L1
# ENSG00000183395.4_PMCH
# ENSG00000187608.5_ISG15
tree.dataset = tree(tissue ~ 
                      ENSG00000133048.8_CHI3L1,
                      data = mydat, subset = train)
tree.pred = predict(tree.dataset, test, type = "class")
tree.table = table(tree.pred, tissue.test)
tree.table
((tree.table[1] + tree.table[4]) / sum(tree.table))

# 0.5421687 ENSG00000133048.8_CHI3L1
# 0.6144578 ENSG00000118271.5_TTR
# 0.7228916 ENSG00000183395.4_PMCH

plot(tree.dataset);text(tree.dataset, pretty = 0)
tree.cv = cv.tree(tree.dataset)
plot(tree.cv$size, tree.cv$dev, type = "b")
#optimal is 3 so prune back
prune.genes = prune.tree(tree.dataset, best = 2)
plot(prune.genes); text(prune.genes, pretty=0)
prediction.prune = predict(prune.genes, test, type = "class")
tree.table.prune = table(prediction.prune, tissue.test)
tree.table.prune
((tree.table.prune[1] + tree.table.prune[4]) / sum(tree.table.prune))

# 0.6506024 ENSG00000118271.5_TTR
# 0.686747 ENSG00000133048.8_CHI3L1
# 0.7228916 ENSG00000183395.4_PMCH

########################################
#CLASSIFICATION TREE WITHOUT THE 3 GENES
########################################

set.seed(2)
n = round(dim(mydat)[1]*.5)

train = sample(1:nrow(mydat), n)
test = mydat[-train,]
tissue.test = mydat$tissue[-train]

#relevant genes
# ENSG00000118271.5_TTR
# ENSG00000133048.8_CHI3L1
# ENSG00000183395.4_PMCH
# ENSG00000187608.5_ISG15
tree.dataset = tree(tissue ~.
                    -ENSG00000133048.8_CHI3L1,
                     data = mydat, subset = train)

tree.pred = predict(tree.dataset, test, type = "class")
tree.table = table(tree.pred, tissue.test)
tree.table
((tree.table[1] + tree.table[4]) / sum(tree.table))

# 0.626506 ENSG00000133048.8_CHI3L1
# 0.7228916 ENSG00000118271.5_TTR
# 0.746988 ENSG00000183395.4_PMCH

plot(tree.dataset);text(tree.dataset, pretty = 0)
tree.cv = cv.tree(tree.dataset)
plot(tree.cv$size, tree.cv$dev, type = "b")
#optimal is 3 so prune back
prune.genes = prune.tree(tree.dataset, best = 3)
plot(prune.genes);text(prune.genes, pretty=0)
prediction.prune = predict(prune.genes, test, type = "class")
tree.table.prune = table(prediction.prune, tissue.test)
tree.table.prune
((tree.table.prune[1] + tree.table.prune[4]) / sum(tree.table.prune)) * 100
# 65.06024 ENSG00000133048.8_CHI3L1
# 68.6747 ENSG00000183395.4_PMCH
# 0.7228916 ENSG00000118271.5_TTR
