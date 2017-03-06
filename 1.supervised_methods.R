######################################################
# Pattern Recognition Project (Week 1)
# Name(s)           Student Number
# Anne Tjallingii:  911024836020
# Diego Montiel:    880505580110
######################################################

#install.packages("ISLR")

library(ISLR)
library(class)
library(MASS)

data = read.table("expr4T.dat", header = TRUE)
tmp_data = data

#####################
#QUESTION 1
#####################
#Thresholds
minsm = 2.5
minm  = 20

# Last column contains tissue label
nn = ncol(data)-1
m = sapply(data[,1:nn],mean)
s = sapply(data[,1:nn],sd)
sm = s/m
par(mfrow = c(1,2))
plot(m[m < 20],main = "m")
plot(sm[sm < 2.5],main = "sm")

# In the following, we make sure that tissue label is kept;
# by giving it an artificial sufficiently high value for m and sm
m = c(m,minm+1)
sm = c(sm,minsm+1)

plot(sm)
plot(m)
#we choose threshold minimum 1 because we need the values to be evenly spread
#WHY??????????????????
dim(data)
length(which(sm > minsm & m > minm))
data = data[,which(sm > minsm & m > minm)]
dim(data)

#checking if the 2 genes are in the data
data[,"ENSG00000271043.1_MTRNR2L2"]
data[,"ENSG00000229344.1_RP5.857K21.7"]

#################
#QUESTION 2
#################

#1st regression with all tissues with 1st gene
set.seed(2)
n = round(dim(data)[1]*.5)
training.set = sample(1:nrow(data), n)
test.set = data[-training.set,]

lm.fit1 = lm(ENSG00000271043.1_MTRNR2L2 ~ .-tissue, subset = training.set, data = data)
par(mfrow = c(2,2))
plot(lm.fit1)
prediction = predict(lm.fit1, test.set)
prediction
par(mfrow = c(1,1))
plot(test.set$ENSG00000271043.1_MTRNR2L2,prediction)
cor.test(test.set$ENSG00000271043.1_MTRNR2L2,prediction)

mean((test.set$ENSG00000271043.1_MTRNR2L2-prediction)^2)
hist(data$ENSG00000271043.1_MTRNR2L2)

#doesnt look normally distributed
#graphs don't look so good, fitted values against standardized residuals is a pattern visible. 
#maybe a data transformation will make a better model. But for now it doesn't look like a good model for our data.
#linear regression doesn't seem like a good model for gene expression data. 
#in this case we use a gene as response >> to find a correlation between the other genes and the one gene which
#is your response

#2b
summary(lm.fit1)
lm.fit1.test = summary(lm.fit1)
x = lm.fit1.test$coefficients[,4]
sum(x < 0.05)
#2c
#What we can see is that a lot of parameter values are negative and positive. 
#Overall the beta's look really small.
#only 5 variables are significant

#1st regression with all tissues with 2nd gene
lm.fit2 = lm(ENSG00000229344.1_RP5.857K21.7 ~ .-tissue, data = data, subset= training.set)
prediction2 = predict(lm.fit2, test.set)
plot(test.set$ENSG00000229344.1_RP5.857K21.7,prediction2)
cor.test(test.set$ENSG00000229344.1_RP5.857K21.7,prediction2)
mean((test.set$ENSG00000229344.1_RP5.857K21.7-prediction2)^2)
#Same as in the first regression with the other gene.
#High mean square error, lots of negative values, doesnt look normally distributed. 
#correlation between test set en prediction looks really high. 
summary(lm.fit2)
par(mfrow = c(2,2))
plot(lm.fit2)

#For getting the names of the tissues
unique(data[dim(data)[2]])

######## Linear regression for 1st gene tissue = brain_amygdala

tissue1 = "brain_amygdala"
data1 = data.frame(data)
mydata1 = data1[which(data$tissue==tissue1),]
mydata1 = droplevels(mydata1)

set.seed(2)
n = round(dim(mydata1)[1]*.5)
training.amygdala = sample(1:nrow(mydata1), n)
test.amygdala = data[-training.amygdala,]

lm.fit3 = lm(ENSG00000229344.1_RP5.857K21.7 ~ ., data = mydata1[-dim(mydata1)[2]], subset = training.amygdala)
prediction.amygdala = predict(lm.fit3, test.amygdala)
plot(test.amygdala$ENSG00000229344.1_RP5.857K21.7,prediction.amygdala)
cor.test(test.amygdala$ENSG00000229344.1_RP5.857K21.7,prediction.amygdala)
mean((test.amygdala$ENSG00000229344.1_RP5.857K21.7-prediction.amygdala)^2)
#correlation en correlation graph between test data en prediction looks really bad. cor = 0.18
#mean square error really high 
summary(lm.fit3)

####### Linear regression for 1st gene tissue = brain_anteriorcortex

tissue2 = "brain_anteriorcortex"
data2 = data.frame(data)
mydata2 = data2[which(data$tissue==tissue2),]
mydata2 = droplevels(mydata2)

set.seed(2)
n = round(dim(mydata2)[1]*.5)
training.anteriorcortex = sample(1:nrow(mydata2), n)
test.anteriorcortex = data[-training.anteriorcortex,]

lm.fit4 = lm(ENSG00000229344.1_RP5.857K21.7 ~.,  data = mydata2[-dim(mydata2)[2]], subset = training.anteriorcortex)
prediction.anteriorcortex = predict(lm.fit3, test.anteriorcortex)
par(mfrow = c(2,2))
plot(test.anteriorcortex$ENSG00000229344.1_RP5.857K21.7,prediction.anteriorcortex)
cor.test(test.anteriorcortex$ENSG00000229344.1_RP5.857K21.7,prediction.anteriorcortex)
mean((test.anteriorcortex$ENSG00000229344.1_RP5.857K21.7-prediction.anteriorcortex)^2)
summary(lm.fit4)

#Add two more lm for the other relevant gene

########Linear regression for 2nd gene = brain_amygdala

tissue1 = "brain_amygdala"
data1 = data.frame(data)
mydata1 = data1[which(data$tissue==tissue1),]
mydata1 = droplevels(mydata1)

set.seed(2)
n = round(dim(mydata1)[1]*.5)
training.amygdala = sample(1:nrow(mydata1), n)
test.amygdala = data[-training.amygdala,]

lm.fit3 = lm(ENSG00000271043.1_MTRNR2L2 ~ ., data = mydata1[-dim(mydata1)[2]], subset = training.amygdala)
prediction.amygdala = predict(lm.fit3, test.amygdala)
plot(test.amygdala$ENSG00000271043.1_MTRNR2L2,prediction.amygdala)
cor.test(test.amygdala$ENSG00000271043.1_MTRNR2L2,prediction.amygdala)
mean((test.amygdala$ENSG00000271043.1_MTRNR2L2-prediction.amygdala)^2)
#correlation en correlation graph between test data en prediction looks really bad. cor = 0.10
#mean square error really high 
summary(lm.fit3)

#######Linear regression for 2nd gene = brain_anteriorcortex
tissue2 = "brain_anteriorcortex"
data2 = data.frame(data)
mydata2 = data2[which(data$tissue==tissue2),]
mydata2 = droplevels(mydata2)

set.seed(2)
n = round(dim(mydata1)[1]*.5)
training.anteriorcortex = sample(1:nrow(mydata2), n)
test.anteriorcortex = data[-training.anteriorcortex,]

lm.fit4 = lm(ENSG00000271043.1_MTRNR2L2 ~.,  data = mydata2[-dim(mydata2)[2]], subset = training.anteriorcortex)
prediction.anteriorcortex = predict(lm.fit4, test.anteriorcortex)
par(mfrow = c(2,2))
plot(lm.fit4)
cor.test(test.anteriorcortex$ENSG00000271043.1_MTRNR2L2,prediction.anteriorcortex)
mean((test.anteriorcortex$ENSG00000271043.1_MTRNR2L2-prediction.anteriorcortex)^2)
summary(lm.fit4)

#2d
log.data = log(data[1:20])
log.data
##################
#Question 3
##################

#Perform logistic regression, LDA and kNN to discriminate several 
#pairs of tissues and/or multiple tissues at once based on their gene expression patterns.
#Here it is using one of the two relevant gene expression levels

#To make a selection with two brain tissues and
#added in a new variable


#----------BRAIN AMYGDALA VS BRAIN HIPPOCAMPUS

#------------------brain_amygdala vs brain_hippocampus
tissue1 = "brain_amygdala"
tissue2 = "brain_hippocampus"
data = data.frame(data)
mydat = data[which(data$tissue == tissue1|data$tissue==tissue2),]
mydat = droplevels(mydat)

set.seed(3)
tissue01 = mydat$tissue
train = sample(1:nrow(mydat), round(dim(mydat)[1]*.5))
test = mydat[-train,]
tissue.test = tissue01[-train]
tissue.train = tissue01[train]

#Logistic regression

glm_amyg_hipp = glm(mydat$tissue ~., data = mydat, subset = train, family = binomial )
summary(glm_amyg_hipp)
#We cannot do a logistic regrression with too many variables

glm.probs = predict(glm_amyg_hipp, test, type = "response")
glm.pred = rep("brain_amygdala", length(glm.probs))
glm.pred[glm.probs > 0.5] = "brain_hippocampus"
table(glm.pred, tissue.test)
mean(glm.pred == tissue.test)
# result : 0.79518071 is correctly predicted

# LDA
lda.amyg.hipp = lda(mydat$tissue ~., data = mydat, subset = train)

lda.pred.train = predict(lda.amyg.hipp)
lda.class.train = lda.pred.train$class
table(lda.class.train, tissue.train)
mean(lda.class.train == tissue.train)
# result : 0.8192771 is correctly predicted

lda.pred = predict(lda.amyg.hipp, test)
table(lda.pred$class, tissue.test)
mean(lda.pred$class == tissue.test)
# result : 0.6987952 is correctly predicted

# KNN
set.seed(4)
#Scaling the data
std.data = scale(mydat[1:20])
training.data = std.data[train,]
testing.data = std.data[-train,]
training.tissue = tissue.train

knn.pred = knn(training.data, testing.data, training.tissue, k = 1)
table(knn.pred, tissue.test)
mean(knn.pred == tissue.test)
# result : 0.6626506 is correctly predicted

knn.pred = knn(training.data, testing.data, training.tissue, k = 3)
table(knn.pred, tissue.test)
mean(knn.pred == tissue.test)
# result : 0.686747 is correctly predicted

knn.pred = knn(training.data, testing.data, training.tissue, k = 10)
table(knn.pred, tissue.test)
mean(knn.pred == tissue.test)
# result : 0.7710843 is correctly predicted


#--------- brain_hippocampus vs brain_nucleusaccumbens
tissue1 = "brain_hippocampus"
tissue2 = "brain_nucleusaccumbens"
data = data.frame(data)
mydat = data[which(data$tissue == tissue1|data$tissue==tissue2),]
mydat = droplevels(mydat)

set.seed(3)
tissue01 = mydat$tissue
train = sample(1:nrow(mydat), round(dim(mydat)[1]*.5))
test = mydat[-train,]
tissue.test = tissue01[-train]
tissue.train = tissue01[train]

#Logistic regression
glm_amyg_hipp = glm(mydat$tissue ~., data = mydat, subset = train, family = binomial )
summary(glm_amyg_hipp)
#We cannot do a logistic regrression with too many variables
glm.probs = predict(glm_amyg_hipp, test, type = "response")
glm.pred = rep("brain_hippocampus", length(glm.probs))
glm.pred[glm.probs > 0.5] = "brain_nucleusaccumbens"
table(glm.pred, tissue.test)
mean(glm.pred == tissue.test)
# result : 0.7184466 is correctly predicted

# LDA
lda.amyg.hipp = lda(mydat$tissue ~., data = mydat, subset = train)

lda.pred.train = predict(lda.amyg.hipp)
lda.class.train = lda.pred.train$class
table(lda.class.train, tissue.train)
mean(lda.class.train == tissue.train)
# result : 0.7211538 is correctly predicted

lda.pred = predict(lda.amyg.hipp, test)
table(lda.pred$class, tissue.test)
mean(lda.pred$class == tissue.test)
# result : 0.6699029 is correctly predicted

# KNN
set.seed(4)
#Scaling the data
std.data = scale(mydat[1:20])
training.data = std.data[train,]
testing.data = std.data[-train,]
training.tissue = tissue.train

knn.pred = knn(training.data, testing.data, training.tissue, k = 1)
table(knn.pred, tissue.test)
mean(knn.pred == tissue.test)
# result : 0.6213592 is correctly predicted

knn.pred = knn(training.data, testing.data, training.tissue, k = 3)
table(knn.pred, tissue.test)
mean(knn.pred == tissue.test)
# result : 0.6699029 is correctly predicted

knn.pred = knn(training.data, testing.data, training.tissue, k = 10)
table(knn.pred, tissue.test)
mean(knn.pred == tissue.test)
# result : 0.6407767 is correctly predicted

#--------- brain_spinalcord vs brain_substantianigra
tissue1 = "brain_spinalcord"
tissue2 = "brain_substantianigra"
data = data.frame(data)
mydat = data[which(data$tissue == tissue1|data$tissue==tissue2),]
mydat = droplevels(mydat)

set.seed(3)
tissue01 = mydat$tissue
train = sample(1:nrow(mydat), round(dim(mydat)[1]*.5))
test = mydat[-train,]
tissue.test = tissue01[-train]
tissue.train = tissue01[train]

#Logistic regression
glm_amyg_hipp = glm(mydat$tissue ~., data = mydat, subset = train, family = binomial )
summary(glm_amyg_hipp)
#We cannot do a logistic regrression with too many variables
glm.probs = predict(glm_amyg_hipp, test, type = "response")
glm.pred = rep("brain_spinalcord", length(glm.probs))
glm.pred[glm.probs > 0.5] = "brain_substantianigra"
table(glm.pred, tissue.test)
mean(glm.pred == tissue.test)
# result 0.6865672 precision accuracy

# LDA
lda.amyg.hipp = lda(mydat$tissue ~., data = mydat, subset = train)

lda.pred.train = predict(lda.amyg.hipp)
lda.class.train = lda.pred.train$class
table(lda.class.train, tissue.train)
mean(lda.class.train == tissue.train)
# result 0.9402985 accuracy

lda.pred = predict(lda.amyg.hipp, test)
table(lda.pred$class, tissue.test)
mean(lda.pred$class == tissue.test)
# result 0.7313433 accuracy

# KNN
set.seed(4)
#Scaling the data
std.data = scale(mydat[1:20])
training.data = std.data[train,]
testing.data = std.data[-train,]
training.tissue = tissue.train

knn.pred = knn(training.data, testing.data, training.tissue, k = 1)
table(knn.pred, tissue.test)
mean(knn.pred == tissue.test)
# 0.6567164

knn.pred = knn(training.data, testing.data, training.tissue, k = 3)
table(knn.pred, tissue.test)
mean(knn.pred == tissue.test)
# 0.641791

knn.pred = knn(training.data, testing.data, training.tissue, k = 10)
table(knn.pred, tissue.test)
mean(knn.pred == tissue.test)
# 0.6567164

#--------- brain_cerebellum vs brain_amygdala
tissue1 = "brain_cerebellum"
tissue2 = "brain_amygdala"
data = data.frame(data)
mydat = data[which(data$tissue == tissue1|data$tissue==tissue2),]
mydat = droplevels(mydat)

set.seed(3)
tissue01 = mydat$tissue
train = sample(1:nrow(mydat), round(dim(mydat)[1]*.5))
test = mydat[-train,]
tissue.test = tissue01[-train]
tissue.train = tissue01[train]

#Logistic regression
glm_amyg_hipp = glm(mydat$tissue ~., data = mydat, subset = train, family = binomial )
summary(glm_amyg_hipp)
#We cannot do a logistic regrression with too many variables
glm.probs = predict(glm_amyg_hipp, test, type = "response")
glm.pred = rep("brain_cerebellum", length(glm.probs))
glm.pred[glm.probs > 0.5] = "brain_amygdala"
table(glm.pred, tissue.test)
mean(glm.pred == tissue.test)
# 0.02020202

# LDA
lda.amyg.hipp = lda(mydat$tissue ~., data = mydat, subset = train)

lda.pred.train = predict(lda.amyg.hipp)
lda.class.train = lda.pred.train$class
table(lda.class.train, tissue.train)
mean(lda.class.train == tissue.train)
# 0.9897959

lda.pred = predict(lda.amyg.hipp, test)
table(lda.pred$class, tissue.test)
mean(lda.pred$class == tissue.test)
# 0.979798

# KNN
set.seed(4)
#Scaling the data
std.data = scale(mydat[1:20])
training.data = std.data[train,]
testing.data = std.data[-train,]
training.tissue = tissue.train

knn.pred = knn(training.data, testing.data, training.tissue, k = 1)
table(knn.pred, tissue.test)
mean(knn.pred == tissue.test)
#0.9999
knn.pred = knn(training.data, testing.data, training.tissue, k = 3)
table(knn.pred, tissue.test)
mean(knn.pred == tissue.test)
# 0.959596

knn.pred = knn(training.data, testing.data, training.tissue, k = 10)
table(knn.pred, tissue.test)
mean(knn.pred == tissue.test)
# 0.9090909
View(data)
dim(mydat)
dim(data)
#--------- brain_cerebellarhemisphere vs brain_cerebellum
tissue1 = "brain_cerebellarhemisphere"
tissue2 = "brain_cerebellum"
data = data.frame(data)
mydat = data[which(data$tissue == tissue1|data$tissue==tissue2),]
mydat = droplevels(mydat)

set.seed(3)
tissue01 = mydat$tissue
train = sample(1:nrow(mydat), round(dim(mydat)[1]*.5))
test = mydat[-train,]
tissue.test = tissue01[-train]
tissue.train = tissue01[train]

#Logistic regression
glm_amyg_hipp = glm(mydat$tissue ~., data = mydat, subset = train, family = binomial )
summary(glm_amyg_hipp)
#We cannot do a logistic regrression with too many variables
glm.probs = predict(glm_amyg_hipp, test, type = "response")
glm.pred = rep("brain_cerebellarhemisphere", length(glm.probs))
glm.pred[glm.probs > 0.5] = "brain_cerebellum"
table(glm.pred, tissue.test)
mean(glm.pred == tissue.test)
# 0.6521739

# LDA
lda.amyg.hipp = lda(mydat$tissue ~., data = mydat, subset = train)

lda.pred.train = predict(lda.amyg.hipp)
lda.class.train = lda.pred.train$class
table(lda.class.train, tissue.train)
mean(lda.class.train == tissue.train)
# 0.7217391

lda.pred = predict(lda.amyg.hipp, test)
table(lda.pred$class, tissue.test)
mean(lda.pred$class == tissue.test)
# 0.6782609

# KNN
set.seed(4)
#Scaling the data
std.data = scale(mydat[1:20])
training.data = std.data[train,]
testing.data = std.data[-train,]
training.tissue = tissue.train

knn.pred = knn(training.data, testing.data, training.tissue, k = 1)
table(knn.pred, tissue.test)
mean(knn.pred == tissue.test)
# 0.5304348

knn.pred = knn(training.data, testing.data, training.tissue, k = 3)
table(knn.pred, tissue.test)
mean(knn.pred == tissue.test)
# 0.5913043

knn.pred = knn(training.data, testing.data, training.tissue, k = 10)
table(knn.pred, tissue.test)
mean(knn.pred == tissue.test)
# 0.6173913

#--------- brain_cortex vs brain_frontalcortex
tissue1 = "brain_cortex"
tissue2 = "brain_frontalcortex"
data = data.frame(data)
mydat = data[which(data$tissue == tissue1|data$tissue==tissue2),]
mydat = droplevels(mydat)

set.seed(3)
tissue01 = mydat$tissue
train = sample(1:nrow(mydat), round(dim(mydat)[1]*.5))
test = mydat[-train,]
tissue.test = tissue01[-train]
tissue.train = tissue01[train]

#Logistic regression
glm_amyg_hipp = glm(mydat$tissue ~., data = mydat, subset = train, family = binomial )
summary(glm_amyg_hipp)
#We cannot do a logistic regrression with too many variables
glm.probs = predict(glm_amyg_hipp, test, type = "response")
glm.pred = rep("brain_cortex", length(glm.probs))
glm.pred[glm.probs > 0.5] = "brain_frontalcortex"
table(glm.pred, tissue.test)
mean(glm.pred == tissue.test)
# 0.5495495

# LDA
lda.amyg.hipp = lda(mydat$tissue ~., data = mydat, subset = train)

lda.pred.train = predict(lda.amyg.hipp)
lda.class.train = lda.pred.train$class
table(lda.class.train, tissue.train)
mean(lda.class.train == tissue.train)
# 0.7657658

lda.pred = predict(lda.amyg.hipp, test)
table(lda.pred$class, tissue.test)
mean(lda.pred$class == tissue.test)
# 0.5495495

# KNN
set.seed(4)
#Scaling the data
std.data = scale(mydat[1:20])
training.data = std.data[train,]
testing.data = std.data[-train,]
training.tissue = tissue.train

knn.pred = knn(training.data, testing.data, training.tissue, k = 1)
table(knn.pred, tissue.test)
mean(knn.pred == tissue.test)
# 0.5135135

knn.pred = knn(training.data, testing.data, training.tissue, k = 3)
table(knn.pred, tissue.test)
mean(knn.pred == tissue.test)
# 0.5405405

knn.pred = knn(training.data, testing.data, training.tissue, k = 10)
table(knn.pred, tissue.test)
mean(knn.pred == tissue.test)
# 0.5585586

#--------- brain_putamen vs brain_cerebellum
tissue1 = "brain_putamen"
tissue2 = "brain_cerebellum"
data = data.frame(data)
mydat = data[which(data$tissue == tissue1|data$tissue==tissue2),]
mydat = droplevels(mydat)

set.seed(3)
tissue01 = mydat$tissue
train = sample(1:nrow(mydat), round(dim(mydat)[1]*.5))
test = mydat[-train,]
tissue.test = tissue01[-train]
tissue.train = tissue01[train]

#Logistic regression
glm_amyg_hipp = glm(mydat$tissue ~., data = mydat, subset = train, family = binomial )
summary(glm_amyg_hipp)
#We cannot do a logistic regrression with too many variables
glm.probs = predict(glm_amyg_hipp, test, type = "response")
glm.pred = rep("brain_putamen", length(glm.probs))
glm.pred[glm.probs > 0.5] = "brain_cerebellum"
table(glm.pred, tissue.test)
mean(glm.pred == tissue.test)
# 0.2072072

# LDA
lda.amyg.hipp = lda(mydat$tissue ~., data = mydat, subset = train)

lda.pred.train = predict(lda.amyg.hipp)
lda.class.train = lda.pred.train$class
table(lda.class.train, tissue.train)
mean(lda.class.train == tissue.train)
# 0.8918919

lda.pred = predict(lda.amyg.hipp, test)
table(lda.pred$class, tissue.test)
mean(lda.pred$class == tissue.test)
# 12

# KNN
set.seed(4)
#Scaling the data
std.data = scale(mydat[1:20])
training.data = std.data[train,]
testing.data = std.data[-train,]
training.tissue = tissue.train

knn.pred = knn(training.data, testing.data, training.tissue, k = 1)
table(knn.pred, tissue.test)
mean(knn.pred == tissue.test)
# 0.8108108

knn.pred = knn(training.data, testing.data, training.tissue, k = 3)
table(knn.pred, tissue.test)
mean(knn.pred == tissue.test)
# 0.8198198
knn.pred = knn(training.data, testing.data, training.tissue, k = 10)
table(knn.pred, tissue.test)
mean(knn.pred == tissue.test)
# 0.7657658

################
#QUESTION 4
###############

#Try including a non-linear term into (at least one of) 
#the classification models you developed in the previous step. 
#What do you observe?

# 0.7108434 = ENSG00000271043.1_MTRNR2L2
# 0.7108434 = ENSG00000229344.1_RP5.857K21.7
# 0.6987952 all genes

#check if this is the right way to do the non-linear term
non.linear = lda(tissue ~.+I(ENSG00000229344.1_RP5.857K21.7^3), 
                 data = data, subset = train)
summary(non.linear)
length(non.linear)
plot(non.linear)
lda.pred = predict(non.linear, test)
mean(lda.pred$class == tissue.test)
table(lda.pred$class, tissue.test)
par(mfrow = c(1,1))
plot(lda.pred$x)