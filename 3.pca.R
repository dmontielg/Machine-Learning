######################################################
# Pattern Recognition Project (Week 3)
# Name(s)           Student Number
# Anne Tjallingii:  911024836020
# Diego Montiel:    880505580110
######################################################

#install.packages("BiplotGUI")
#install.packages("cluster")
#install.packages("ISLR")

library(ISLR)
library(rgl)
library(ISLR)
library(MASS)
library(class)
library(cluster)
library(BiplotGUI)

data = read.table("expr4T.dat", header = TRUE)

#####################
# Principal settings
#####################

# Filter genes based on minimum mean and stdev/mean
#minsm: 2 = 10, 1.8 = 28, 1.9 = 25

minsm= 2.5
minm= 20

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
dim(data)
length(which(sm > minsm & m > minm))
data = data[,which(sm > minsm & m > minm)]
dim(data)
tissues = data$tissue

#checking if the 2 genes are in the data
data[,"ENSG00000271043.1_MTRNR2L2"]
data[,"ENSG00000229344.1_RP5.857K21.7"]

#########################################
# 1. PCA
#########################################

#This apply a function to a set of rows

#Calculate the mean and variance of the genes
apply(data[1:20],2, mean)
apply(data[1:20],2, var)

#Just give the complete matrix with a scale = true

# Scale
par(mfrow=c(1,2))
pr.out = prcomp(data[1:20], scale = TRUE)
plot(pr.out$x)

#biplot(pr.out, scale = 0)
# Biplot 
Biplots(data[1:20], PointLabels = data$tissue)

#Not scale
pr.out = prcomp(data[1:20], scale = FALSE)
plot(pr.out$x)
# pr.out$sdev
# biplot(pr.out, scale = 1)

# The center and scale components correspondd to the means and se
#of the variables that were used for scaling prior to implementing PCA

#A
summary(pr.out)
pve = 100*pr.out$sdev^2/sum(pr.out$sdev^2)
par(mfrow=c(1,2))
plot(pve, type="o", ylab="Proportion of Variance explained", xlab = "PC", col="blue")
plot(cumsum(pve), type="o", 
     ylab="Proportion of cumulative variance explained", xlab = "Acum PC", col="brown")

#B
set.seed(1)
labels = as.character(unique(tissues))
cols = rainbow(length(unique(tissues)))
colours = (cols[as.numeric(as.factor(tissues))])

data.tissues = data$tissue
par(mfrow =c(1,1))
# PCA1 vs PCA2: PCA2 explains more about green labels, 
# brain cerebellarhemisphere,cerebellum and cortex
plot(pr.out$x[,c(1,2)], col =colours, pch =19,xlab ="Z1",ylab="Z2");
legend(10, -12, c(labels), pch=19, col=c(cols), xjust = 0, cex = 0.7)
# PCA1 vs PCA3: PCA3 show most of the tissues clustered nicely together
plot(pr.out$x[,c(1,3) ], col =colours, pch =19,xlab ="Z1",ylab="Z3");
legend(10, -10, c(labels), pch=19, col=c(cols), xjust = 0, cex = 0.7)
# PCA1 vs PCA4: PCA4 explains more about blue labels.
# brain hypothalamus and brain nucleusaccumbens
plot(pr.out$x[,c(1,4)], col = colours, pch =19,xlab ="Z1",ylab="Z4");
legend(10, -3, c(labels), pch=19, col=c(cols), xjust = 0, cex = 0.6)
#PCA1 vs PCA5
plot(pr.out$x[,c(1,5) ], col =colours, pch =19,xlab ="Z1",ylab="Z5");
legend(10, 15, c(labels), pch=19, col=c(cols), xjust = 0, cex = 0.7)

# c which genes contribute substantially to first or second PC
relevant.genes <- pr.out$rotation[which(abs(pr.out$rotation[,1]) > 0.2 | abs(pr.out$rotation[,2]) >0.2), 1:2]
par(mfrow=c(1,1))
row.names(relevant.genes)

set.seed(1)
genes = colnames(data, do.NULL = TRUE, prefix = "col")[1:20]
labels = genes
cols = rainbow(length(genes))
colours = (cols[as.numeric(as.factor(genes))])

plot(relevant.genes, pch = 1, main="Most contributed genes to PC1 and PC2")
legend(0.2, 0, row.names(relevant.genes), pch=1, xjust = 0, cex = 0.7)
#identify(relevant.genes, n = 1, label = row.names(relevant.genes))

#########################################
# 2. Hierarchical Clustering
#########################################
#a)

#Run if this variable has not been create it from previous section
pr.out = prcomp(data[1:20], scale = TRUE)

sd.data = data[1:20]
par(mfrow = c(1,1))

# Euclidean complete
hc.euclidean.complete = hclust(dist(sd.data, method = "euclidean"),method = "complete")
plot(hc.euclidean.complete, main = "Complete Linkage (Euclidean)", cex = 0.5)
hc.euclidean.complete.cluster = cutree(hc.euclidean.complete, 13)
table(hc.euclidean.complete.cluster, tissues)
hist(hc.euclidean.complete.cluster)

# Euclidean average
hc.euclidean.average  = hclust(dist(sd.data, method = "euclidean"),method = "average")
plot(hc.euclidean.average, main = "Average Linkage (Euclidean)", cex = 0.5)
hc.euclidean.average.cluster = cutree(hc.euclidean.average,13)
table(hc.euclidean.average.cluster, tissues)
hist(hc.euclidean.average.cluster)

# Correlation pearson complete
dist_gene = as.dist(1-cor(t(data[,-21]), method = "pearson"))
hc.correlation.complete = hclust(dist_gene, method = "complete")
plot(hc.correlation.complete, main = "Complete Linkage (Pearson)", cex = 0.5)
hc.correlation.complete.cluster = cutree(hc.correlation.complete, 13)
table(hc.correlation.complete.cluster, tissues)
hist(hc.correlation.complete.cluster)

# Correlation pearson average
dist_gene = as.dist(1-cor(t(data[,-21]), method = "pearson"))
hc.correlation.average = hclust(dist_gene, method = "average")
plot(hc.correlation.average, main = "Average Linkage (Pearson)", cex = 0.5)
hc.correlation.average.cluster = cutree(hc.correlation.average, 13)
table(hc.correlation.average.cluster, tissues)
hist(hc.correlation.average.cluster)

# K means clustering vs Hierarchical clustering
set.seed(2)
km.out = kmeans(sd.data, 13, nstart = 10)
km.clusters = km.out$cluster
table(km.clusters, hc.euclidean.complete.cluster)
plot(km.clusters, hc.euclidean.complete.cluster, col = colours, pch =19,xlab ="K-means clusters",
     ylab="Hierarchical clusters" )
legend(3, 13, c(labels), pch=19, col=c(cols), xjust = 0, cex = 0.6)

# k means with 13 clusters
new.data = t(sd.data)
km.out = kmeans(new.data,13)
plot(new.data, col = km.out$cluster,
     main="k-means clustering with k=13",
     pch=20, cex=2)
legend("bottomright", legend = levels(data$tissue),fill=1:13, cex=0.5)

# applying pc vectors to k means clustering 
pcdata = data[1:20]
pca_km <- princomp(pcdata, cor = T)
pc.comp <- pca_km$scores
pc.comp1 <- -1*pc.comp[,1]
pc.comp2 <- -1*pc.comp[,2]
x <- cbind(pc.comp1,pc.comp2)
cl <- kmeans(x,13)
plot(pc.comp1,pc.comp2,col=data$tissue, main="PC Vector on k means clustering")
points(cl$centers, pch = 16)
legend("topright", legend = levels(data$tissue),fill=1:13, cex=0.45)

# clusplot(pcdata, cl$cluster, color=TRUE, shade = TRUE, labels = 2, lines = 0, cex=0.7)

###################################################
# 3.Retrain using PCA with the gene expression data 
#   by using LDA and KNN
###################################################

#hippocampus and nucleusaccumbens

# Train and test set
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

# Logistic regression

glm_amyg_hipp = glm(tissue ~ 
                      ENSG00000225972.1_MTND1P23+
                      ENSG00000237973.1_hsa.mir.6723+
                      ENSG00000229344.1_RP5.857K21.7+
                      ENSG00000271043.1_MTRNR2L2+
                      ENSG00000120738.7_EGR1+
                      ENSG00000204388.5_HSPA1B+
                      ENSG00000255633.3_MTRNR2L9+     
                      ENSG00000255823.1_MTRNR2L8 +    
                      ENSG00000132002.3_DNAJB1, 
                      data = mydat, subset = train, family = binomial )

summary(glm_amyg_hipp)
glm.probs = predict(glm_amyg_hipp, test, type = "response")
glm.pred = rep("brain_cerebellarhemisphere", length(glm.probs))
glm.pred[glm.probs > 0.5] = "brain_cerebellum"
table(glm.pred, tissue.test)
mean(glm.pred == tissue.test)

# LDA 
lda.amyg.hipp = lda(tissue ~ 
                      ENSG00000225972.1_MTND1P23+
                      ENSG00000237973.1_hsa.mir.6723+
                      ENSG00000229344.1_RP5.857K21.7+
                      ENSG00000271043.1_MTRNR2L2+
                      ENSG00000120738.7_EGR1+
                      ENSG00000204388.5_HSPA1B+
                      ENSG00000255633.3_MTRNR2L9+     
                      ENSG00000255823.1_MTRNR2L8 +    
                      ENSG00000132002.3_DNAJB1, 
                      data = mydat, subset = train)

lda.pred = predict(lda.amyg.hipp, test)
mean(lda.pred$class == tissue.test)
table(lda.pred$class, tissue.test)

# KNN

set.seed(4)
std.data = mydat[1:20]
training.data = std.data[train,]
testing.data = std.data[-train,]
training.tissue = tissue.train
knn.pred = knn(training.data, testing.data, training.tissue, k = 2)
mean(knn.pred == tissue.test)
table(knn.pred, tissue.test)
