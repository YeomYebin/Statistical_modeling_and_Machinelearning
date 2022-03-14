## 기말고사 MASTER CODE##

#################################################################
############### CH 2 ############################################
######## MISSING VALUE ##############

install.packages("tidyverse")
library(tidyverse)
library(mice)
library(rms)
library(finalfit)

# Data
head(nhanes)
nhanes$age = as.factor(nhanes$age)
nhanes$hyp = as.factor(nhanes$hyp)

# Description of data
describe(nhanes)

# Missing pattern
md.pattern(nhanes) ##미싱 패턴 확인

missing_pairs(nhanes,'chl',c('age','bmi','hyp') ) ##(data, y, x) ## Y변수의 miss에 대해 시각화

# Clustering for variables with missing values for the same obs.
# Missing pattern between variables 
# (rr: both observed, mm: both missing, 
#  mr: row variable is missing & column variable is observed)
md.pairs(nhanes) ##관측, NA 조합에 대해서 계산해줌

missing.clus = naclus(nhanes, method='average') ## NA로 코딩되어 있어야 함
missing.clus
plot(missing.clus)

# Check missing pattern of Y variable. ## Y가 NA인지 아닌지에 대해 Logistic Regression ## 유의한지 안한지..?
fit = glm(is.na(chl) ~ age + bmi + hyp, data=nhanes, family=binomial )
summary(fit)

########### MICE #################3
# Missing values should be coded as NA.
# m: The number of imputed datasets
# method: Imputation methods for each column.
# predictorMatrix: A matrix containing 0/1 data specifying 
#                 the set of predictors to be used for each target column.
# 

set.seed(1234) ##씨드고정!!

nhanes %>% head()
imp = mice(nhanes, m=10, method=c('','pmm','logreg','pmm'), print=F)

imp$predictorMatrix

pred = imp$predictorMatrix ## Predict Matrix에 Y를 0으로 넣어줌
pred[,'chl'] = 0
pred

imp1 = mice(nhanes, m=10, method=c('','pmm','logreg','pmm'), 
            predictorMatrix = pred, print=F)

imp1$predictorMatrix

# list the actual imputations for BMI
imp$imp$bmi

# The first imputed dataset
complete(imp,1)

# Checking the convergence of MICE.
plot(imp, c('bmi','hyp','chl')) ## 수렴 확인!!!!!!

# Compare the imputed data and observed data.
stripplot(imp, pch=20, cex=1.2)
# blue point: observed, red point: imputed

xyplot(imp,chl ~ bmi | .imp)
# The first plot: original complete set.

densityplot(imp, scales=list(relation='free'),layout=c(1, 2))
# blue line: density for observed data, red line: density for imputed data


# Prediction model with MICE
# Goal: predict chl based on age, hyp, bmi variables

set.seed(1234)

imp = mice(nhanes, m=10, method=c('','pmm','logreg','pmm'), print=F)
# In mice, all variables with missing values should be imputed, 
# even if it is Y variable.

# To delete obs with missing Y value, imputed Y is replaced with NA.
md = dim(imp$imp$chl)

iy = as.data.frame(matrix(NA,md[1],md[2]))
colnames(iy) = colnames(imp$imp$chl)
rownames(iy) = rownames(imp$imp$chl)

imp$imp$chl = iy


# Apply prediction model to each imputed dataset.
# E.g., prediction model => linear regression model

fit = with(imp, lm(chl ~ age + bmi + hyp))

# Model averaging.
summary(pool(fit))


# Checking imputation effect for significant X variables 
# E.g., to impute bmi variable, we used chl (Y) variable.

comp.dat = na.omit(nhanes)
fit1 = lm(chl ~ age + bmi + hyp, data=comp.dat)
summary(fit1)


############ Predict test obs.###################3
M = imp$m
imp.dat = vector(mode='list',length=M)
for (m in 1:M) imp.dat[[m]] = complete(imp,m)

p.model = function(dat) lm(chl ~ age + bmi + hyp, data=dat) ## 모델 넣기 ## GLM으로 바꿀 수 있을 듯...?

fit.imp = lapply(imp.dat, p.model)

test.obs = data.frame(age=c('2','1'),bmi = c(23.3,21.5),hyp=c('1','1'))

yhat = lapply(fit.imp, predict, newdata=test.obs)

yhat = matrix(unlist(yhat), nrow(test.obs), M)

yhat

apply(yhat,1,mean) ## 1:같은 행별, 2:같은 열별

####################################
##### Data transformation #############
####################################

# install.packages('rms')
library(rms)
library(e1071)

getHdata(titanic3)
dat = titanic3[,c('survived','pclass','age','sex','sibsp','parch')]

describe(dat)

is.na(dat) %>% apply(2, sum) ## 열별 NA개수 확인

md.pattern(dat)

imp = mice(dat,m=1,method='pmm')
imp.dat = complete(imp)

par(mfrow=c(1,3))
for (j in c('age','sibsp','parch')) {
  hist(imp.dat[,j], main=j,xlab = skewness(imp.dat[,j]))  
} ## skewness가 + : Positive skewness -> 정규분포보다 왼쪽으로 치우쳐 있음, Right-skewed되어 있다.

# Standardization or centering: Use 'scale()' function.
scale(imp.dat$age, cneter = TRUE, scale = TRUE) ##표준화

# Yeo-Johnson transformation
# install.packages('bestNormalize')
library(bestNormalize)

imp.dat1 = imp.dat
imp.dat1$sibsp = yeojohnson(imp.dat$sibsp)$x.t
imp.dat1$parch = yeojohnson(imp.dat$parch)$x.t

par(mfrow=c(1,2))
for (j in c('sibsp','parch')) {
  hist(imp.dat1[,j], main=j,xlab = skewness(imp.dat1[,j]))  
}

# Discretization of continuous variable

imp.dat2 = transform(imp.dat,
                     agec = ifelse(age < 21, 'child','adult'),
                     sibsp = ifelse(sibsp==0, 'no sibsp','sibsp'),
                     parch = ifelse(parch==0, 'no parch','parch')) ## 새로운 열 추가

head(imp.dat2)

##############################################################
########### CH 4. Dimension Reduction ########################
##############################################################

options(warn = -1)  # Turn off warning message

# Data: Breast Cancer Wisconsin Data

dat = read.csv("C:/Users/User/Desktop/통계적모델링과머신러닝/wdbc.csv", header = F)

x.name = c("radius", "texture", "perimeter", "area", "smoothness", 
           "compactness", "concavity", "concave_points", "symmetry", 
           "fractal_dimension")

names(dat) = c("id", "diagnosis", paste0(x.name,"_mean"), 
               paste0(x.name,"_se"), paste0(x.name,"_worst"))

dat = dat[,-1]

head(dat)

############### Principal Component Analysis (PCA) ###############

pr = prcomp(dat[,2:31], center = TRUE, scale = TRUE)
summary(pr)

# Scree plot
screeplot(pr, type = 'l', npcs = 15, main = 'Scree plot') ##일부 X에 Correlation 존재

pr$x[,1:3] ## PC 뽑기!!!!!!!!!!!!!!!!!!
predict(pr,dat[,2:31] )[,1:3] ## NEW data에 대해서!!!

##### Visulalization: Scatter plot matrix ########

library(lattice)

pc.dat = data.frame(type = dat$diagnosis, pr$x[,1:3])
pc.dat$type = as.factor(pc.dat$type)

splom(~pc.dat[,2:4], groups=type, data=pc.dat, panel=panel.superpose) ## 조건부 산점도 그래프

library(car)
scatterplotMatrix(~PC1+PC2+PC3|type, data=pc.dat)

library(GGally)
ggpairs(pc.dat, aes(colour = type, alpha = 0.4))


# Application to logistic regression:
library(boot)
dat$diagnosis = as.factor(dat$diagnosis)

# Comparison of CV statistics:
set.seed(1234)
fit = glm(diagnosis~., data=dat, family=binomial)
cv.glm(dat, fit, K=5)$delta

?cv.glm

fit1 = glm(type~., data=pc.dat, family=binomial)
cv.glm(pc.dat, fit1, K=5)$delta

############## Principal Curve ##############
library(princurve)

# Simple example
set.seed(1234)
n=100
x1 = runif(n,-1,1)
x2 = x1^2 + rnorm(n, sd = 0.1)
z = cbind(x1,x2)

cor(x1,x2)

fit = principal_curve(z)

plot(x1,x2)
lines(sort(x1),(sort(x1))^2,col='red',lty=3)
lines(fit ,col='blue')
whiskers(z, fit$s)

# WDBC data application:

fit = principal_curve(as.matrix(dat[,2:31]))

# Density of two groups in principal curve and PC1
par(mfrow=c(1,2))
plot(density(fit$lambda[dat$diagnosis == 'B']), col='red',
     xlim=c(500,5500),main='Principal Curve')
lines(density(fit$lambda[dat$diagnosis == 'M']),col='blue')

plot(density(pc.dat[dat$diagnosis == 'B',2]),col='red',
     xlim=c(-15,7),main='PC1')
lines(density(pc.dat[dat$diagnosis == 'M',2]),col='blue')

# Density for principal curve and PC1
par(mfrow=c(1,2))
plot(density(fit$lambda),col='red',main='Principal Curve')
plot(density(pc.dat[,2]),col='red',main='PC1')

dat1 = cbind(dat, pcurve=fit$lambda)
dat1$diagnosis = as.factor(dat1$diagnosis)

fit2 = glm(diagnosis~pcurve, data=dat1, family=binomial)
cv.glm(dat1, fit2, K=5)$delta

# Projection onto Principal curve
new.obs = as.matrix(dat[31:40,2:31])
project_to_curve(new.obs,fit$s)$lambda  
# arc-length along the curve.

################ Kernel PCA ################
install.packages('kernlab')
library(kernlab)

x = dat[,2:31]
fit = kpca(~., data=x, kernel='rbfdot', kpar=list(sigma=3), features=2)
# feature: # of PC's

# Kernel PC's
pc = pcv(fit)

B = pc[dat$diagnosis=='B',]
M = pc[dat$diagnosis=='M',]

par(mfrow=c(1,1))
plot(B, col='red', xlab='KPC1',ylab='KPC2')
points(M, col='blue')

# New observations
predict(fit, new.obs)


############ Non-negative matrix factorization ############ 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biobase")

install.packages('NMF')
library(NMF)

?esGolub

data(esGolub)

dim(esGolub)

res = nmf(esGolub, rank = 3, seed=123456)

# W matrix
W = basis(res)
dim(W)

# H matrix
H = coef(res)
dim(H)

if(requireNamespace("Biobase", quietly=TRUE))
{
  estim.r = nmf(esGolub, 2:6, nrun=10, seed=123456)
  plot(estim.r)
}

res = nmf(esGolub, rank = 5, seed=123456)

# Visualization
basismap(res, subsetRow=TRUE)

coefmap(res)

############ Independent component analysis ############

install.packages('fastICA')
library(fastICA)


# Ex1:

S = matrix(runif(10000), 5000, 2)
A = matrix(c(1, 1, -1, 3), 2, 2, byrow = TRUE)
X = S %*% A
a = fastICA(X, 2, alg.typ = "parallel", fun = "logcosh", alpha = 1,
            method = "C", row.norm = FALSE, maxit = 200,
            tol = 0.0001, verbose = TRUE)
par(mfrow = c(1, 3))
plot(a$X, main = "Pre-processed data")
plot(a$X %*% a$K, main = "PCA components")
plot(a$S, main = "ICA components")


# Ex2:

S = cbind(sin((1:1000)/20), rep((((1:200)-100)/100), 5))
A = matrix(c(0.291, 0.6557, -0.5439, 0.5572), 2, 2)
X = S %*% A
a = fastICA(X, 2, alg.typ = "parallel", fun = "logcosh", alpha = 1,
            method = "R", row.norm = FALSE, maxit = 200,
            tol = 0.0001, verbose = TRUE)

par(mfcol = c(2, 3))
plot(1:1000, S[,1], type = "l", main = "Original Signals",
     xlab = "", ylab = "")
plot(1:1000, S[,2], type = "l", xlab = "", ylab = "")
plot(1:1000, X[,1], type = "l", main = "Mixed Signals",
     xlab = "", ylab = "")
plot(1:1000, X[,2], type = "l", xlab = "", ylab = "")
plot(1:1000, a$S[,1], type = "l", main = "ICA source estimates",
     xlab = "", ylab = "")
plot(1:1000, a$S[, 2], type = "l", xlab = "", ylab = "")


##############################################################
########### CH 5. Variable Selection ########################
##############################################################

options(warn = -1)  # Turn off warning message

######### Variable Importance: Regression problem ##########

# Building data
# 90 economic variables and sales variable (output)
#install.packages("dplyr")
dat = read.table("C:/Users/User/Desktop/통계적모델링과머신러닝/building.csv", sep=',', header=T)

# Correlation coefficient ---------------------------
VI = cor(dat)[,'price']
SVI = sort(abs(VI), decreasing = T)[-1]
SVI

par(mfrow = c(1, 1))
plot(1:length(SVI),SVI, type='b', ylab='Size of correlation',
     xlab ='variables', main='Variable Importance', xaxt='n')
axis(side=1, at=1:length(SVI), labels=names(SVI), cex.axis=0.3,las=2)

# Spearman rank correlation coefficient -------------
SP = cor(dat, method='spearman')[,'price']
SPI = sort(abs(SP), decreasing = T)[-1]
SPI

plot(1:length(SPI),SPI, type='b', ylab='Size of correlation',
     xlab ='variables', main='Variable Importance', xaxt='n')
axis(side=1, at=1:length(SPI), labels=names(SPI), cex.axis=0.3,las=2)

# Pseudo R^2 ----------------------------------------
p = 90
PR2 = numeric(p)
names(PR2) = colnames(dat[,-91])
for (j in 1:p)
{
  fit = loess(price ~ dat[,j], data=dat)  # Local linear regression
  yhat = predict(fit, dat[,j])
  PR2[j] = 1-(sum((dat$price - yhat)^2)/sum((dat$price - mean(dat$price))^2))
}

SPR2 = sort(PR2, decreasing = T)
SPR2

plot(1:length(SPR2),SPR2, type='b', ylab='Pseudo R2',
     xlab ='variables', main='Variable Importance', xaxt='n')
axis(side=1, at=1:length(SPR2), labels=names(SPR2), cex.axis=0.3,las=2)


# Maximal information coefficient (MIC) -------------

#install.packages('minerva')
library(minerva)

MIC = mine(dat)
MIC = MIC$MIC[,'price']

SMIC = sort(MIC, decreasing = T)[-1]
SMIC

plot(1:length(SMIC),SMIC, type='b', ylab='MIC',
     xlab ='variables', main='Variable Importance', xaxt='n')
axis(side=1, at=1:length(SMIC), labels=names(SMIC), cex.axis=0.3,las=2)


######### Variable Importance: Classification problem ##########

# Data
#install.packages('mlbench')
library(mlbench)
data(BreastCancer)
dat = BreastCancer[,-1]

# Relief algorithm ------------------------------

#install.packages('CORElearn')
library(tidyverse)
library(CORElearn)

# Relief algorithm ## with classification!!!!!!
dat %>% head()
RE = attrEval(Class ~ ., data=dat, estimator='Relief',
              ReliefIterations=30)

SRE = sort(RE, decreasing = T)
SRE

plot(1:length(SRE),SRE, type='b', ylab='Separability',
     xlab ='variables', main='Variable Importance', xaxt='n')
axis(side=1, at=1:length(SRE), labels=names(SRE), cex.axis=0.8,las=2)


########### ReliefF algorithm #################
REF = attrEval(Class ~ ., data=dat, estimator='ReliefFequalK',
               ReliefIterations=30) ## estimator가 바뀜

SREF = sort(REF, decreasing = T)
SREF

plot(1:length(SREF),SREF, type='b', ylab='Separability',
     xlab ='variables', main='Variable Importance', xaxt='n')
axis(side=1, at=1:length(SREF), labels=names(SREF), cex.axis=0.8,las=2)


########## Variable Selection: Simulated Annealing ##########

#install.packages('mvtnorm')
library(mvtnorm)

# Data generation
set.seed(10)

n = 500
p = 20
S = matrix(0.3, nrow=p, ncol=p)
diag(S) = 1
X = rmvnorm(n, mean=rep(0,p), sigma=S)

XN = NULL
for (j in 1:p) XN = c(XN,paste('X',j,sep=''))
colnames(X) = XN

Y = 2 + 0.5*X[,1] - 0.3*X[,2] + 1.2*X[,3] + rnorm(n,sd=0.1)

# Simulated Annealing ---------------------------

#install.packages('caret')
library(caret)

ctrl = safsControl(functions=caretSA, method='cv', number=5)

obj = safs(x=X, y=Y, iters=20, safsControl=ctrl, method='lm')
obj$fit
obj

#################### ISIS ####################

#install.packages('SIS')
library(SIS)

?SIS

# Data generation
set.seed(0)
n = 400; p = 50; rho = 0.5
corrmat = diag(rep(1-rho, p)) + matrix(rho, p, p)
corrmat[,4] = sqrt(rho)
corrmat[4, ] = sqrt(rho)
corrmat[4,4] = 1
corrmat[,5] = 0
corrmat[5, ] = 0
corrmat[5,5] = 1
cholmat = chol(corrmat)
x = matrix(rnorm(n*p, mean=0, sd=1), n, p)
x = x%*%cholmat

# Linear regression
set.seed(1)
b = c(4,4,4,-6*sqrt(2),4/3)
y=x[, 1:5]%*%b + rnorm(n)


# ISIS with regularization
model11=SIS(x, y, family='gaussian', tune='bic', seed=1234)
model11$ix

model12=SIS(x, y, family='gaussian', tune='bic', varISIS='aggr', seed=1234) ## Vanilia SIS
model12$ix
model12$lambda

# logistic regression
set.seed(2)
feta = x[, 1:5]%*%b; fprob = exp(feta)/(1+exp(feta))
y = rbinom(n, 1, fprob)

# ISIS with regularization
model21=SIS(x, y, family='binomial', tune='bic', penalty='SCAD', perm=T, q=0.9)
model21$ix

model22=SIS(x, y, family='binomial', tune='bic', varISIS='aggr', seed=21)
model22$ix

######################################################
############# CH5. Imbalenaced Data Modeling##########
######################################################

options(warn = -1)  # Turn off warning message

# Data generation
#install.packages('mvtnorm')
library(mvtnorm)

set.seed(10)
n1 = 500; n2 = 50
mean1 = c(6,6)
mean2 = c(4,5)
sig1 = matrix(c(2,0,0,2),2,2)
sig2 = matrix(c(2,-0.8,-0.8,2),2,2)

x = rbind(rmvnorm(n1,mean1,sig1),rmvnorm(n2,mean2,sig2))
y = as.factor(c(rep(0,n1),rep(1,n2)))
train = data.frame(y,x)

x = rbind(rmvnorm(n1,mean1,sig1),rmvnorm(n2,mean2,sig2))
y = as.factor(c(rep(0,n1),rep(1,n2)))
test = data.frame(y,x)


# Scatter plot for the training data
plot(cbind(train$X1,train$X2),
     col=(3-as.numeric(train$y)),xlab='X1',ylab='X2')

# Classifier: Logistic regression 
fit = glm(y~., family=binomial, data=train)
phat.test = predict(fit, test, type='response')
yhat.test = ifelse(phat.test > 0.5, 1, 0)

cm  = table(true = test$y, predict=yhat.test)
cm


misclass = function(cm) 1 - sum(diag(cm))/sum(cm)
# cm: confusion matrix (higher value = positive class)
fmeasure = function(cm)
{
  TPR = cm[2,2]/sum(cm[2,])
  PPV = cm[2,2]/sum(cm[,2])
  return((2*TPR*PPV)/(TPR + PPV))
}

# Test misclassification rate
misclass(cm)

# Test F-measure
fmeasure(cm)

library(caret)
confusionMatrix(as.factor(yhat.test), as.factor(test$y), mode = "everything", positive="1")

gdata::drop.levels(test$y)


######## ROC curve ################
library(pROC)

phat.tr = predict(fit, train, type='response')

lr.roc = roc(train$y ~ phat.tr)
plot(lr.roc)
library(Epi)
ROC(test = phat.tr, stat = train$y, plot = "ROC", AUC = T, main = "Logistics Regression") ###멋진 ROC커브 그려줌

# AUC
auc(lr.roc)


########## Alternate cut-off (using training dataset) ##########
# Fird the closest point on ROC curve to (1,1).
th = coords(lr.roc,x='best', best.method = 'closest.topleft')
th

# Evaluation for the new cut-off
yhat.test1 = ifelse(phat.test > th$threshold, 1, 0)

cm1 = table(true = test$y, predict=yhat.test1)
cm1

# Test misclassification rate
misclass(cm1)

# Test F-measure
fmeasure(cm1)

################### Adjusting prior prob. ####################
library(MASS)

fit = lda(y~., data=train)
yhat.te = predict(fit,test)$class

cm = table(true = test$y, predict=yhat.te)
cm

# Test misclassification rate
misclass(cm)

# Test F-measure
fmeasure(cm)

# Adjust prior prob.
fit1 = lda(train$y, x = as.matrix(train[,-1]), prior=c(0.6,0.4))
yhat.te1 = predict(fit1 ,as.matrix(test[,-1]))$class

cm1 = table(true = test$y, predict=yhat.te1)
cm1

# Test misclassification rate
misclass(cm1)

# Test F-measure
fmeasure(cm1)

#################### Sampling methods ####################
#install.packages('caret')
library(caret)

# SVM model from the original imbalaned data
#install.packages('e1071')
library(e1071)

cv.fit = tune(svm, y~., data=train, kernel='radial', 
              ranges=list(cost=c(0.001,0.01,0.1,1,5,10,100),
                          gamma=c(0.01,0.1,0.5,1,2,3,4)))
summary(cv.fit)

# SVM with the best parameters
best.fit = cv.fit$best.model

# Prediction for test data
yhat.te = predict(best.fit, test)
cm = table(true = test$y, predict=yhat.te)

misclass(cm)

fmeasure(cm)


################## Upsampling #################
# Sampling with replacement from the small class.
set.seed(10)
uptrain = upSample(x = train[,-1], y=train$y, yname='y')

dim(uptrain)

table(uptrain$y)

# Scatter plot for the upsampled training data
plot(cbind(uptrain$X1,uptrain$X2),
     col=(3-as.numeric(uptrain$y)),xlab='X1',ylab='X2')


# SVM for upsampled data.
set.seed(10)
cv.fit1 = tune(svm, y~., data=uptrain, kernel='radial', 
               ranges=list(cost=c(0.001,0.01,0.1,1,5,10,100),
                           gamma=c(0.01,0.1,0.5,1,2,3,4)))
summary(cv.fit1)

# SVM with the best parameters
best.fit1 = cv.fit1$best.model

# Prediction for test data
yhat.te1 = predict(best.fit1, test)
cm1 = table(true = test$y, predict=yhat.te1)

misclass(cm1)

fmeasure(cm1)


################ Dwonsampling ################
# Randomly remove obs. in the large class.
set.seed(1)
dntrain = downSample(x = train[,-1], y=train$y, yname='y')

dim(dntrain)

table(dntrain$y)

# Scatter plot for the downsampled training data
plot(cbind(dntrain$X1,dntrain$X2),
     col=(3-as.numeric(dntrain$y)),xlab='X1',ylab='X2')


# SVM for downsampled data.
set.seed(10)
cv.fit2 = tune(svm, y~., data=dntrain, kernel='radial', 
               ranges=list(cost=c(0.001,0.01,0.1,1,5,10,100),
                           gamma=c(0.01,0.1,0.5,1,2,3,4)))
summary(cv.fit2)

# SVM with the best parameters
best.fit2 = cv.fit2$best.model

# Prediction for test data
yhat.te2 = predict(best.fit2, test)
cm2 = table(true = test$y, predict=yhat.te2)

misclass(cm2)

fmeasure(cm2)

########################### SMOTE ###########################
install.packages('DMwR')
library(DMwR)

set.seed(1)
smtrain = SMOTE(y~., data=train, perc.over=200, k = 5, perc.under=200)

dim(smtrain)
table(smtrain$y)


# Scatter plot for the training data from SMOTE
plot(cbind(smtrain$X1,smtrain$X2),
     col=(3-as.numeric(smtrain$y)),xlab='X1',ylab='X2')



# SVM for data from SMOTE.
set.seed(10)
cv.fit3 = tune(svm, y~., data=smtrain, kernel='radial', 
               ranges=list(cost=c(0.001,0.01,0.1,1,5,10,100),
                           gamma=c(0.01,0.1,0.5,1,2,3,4)))
summary(cv.fit3)

# SVM with the best parameters
best.fit3 = cv.fit3$best.model

# Prediction for test data
yhat.te3 = predict(best.fit3, test)
cm3 = table(true = test$y, predict=yhat.te3)

misclass(cm3)

fmeasure(cm3)

################# One-class learning ##################
# Support vector data description (SVDD) -----------------

# Training data for the large class.
train.x0 = train[train$y == 0,-1]

result = NULL
for (nu in c(0.001,0.01,0.1,0.3,0.5,0.7,0.9))
{
  for (gamma in c(0.01,0.1,0.5,1,2,3,5))
  {
    svddfit = svm(x=train.x0, type='one-classification', kernel='radial',
                  nu=nu, gamma=gamma)
    
    minor = predict(svddfit,test[,-1])
    # predict: TRUE: small class(outlier), FALSE: large class
    yhat.te = numeric(nrow(test))
    yhat.te[minor==TRUE] = 1
    cm = table(true = test$y, predict=yhat.te)
    result = rbind(result, c(nu,gamma,fmeasure(cm)))
  }
}
# Best result for F-measure.
names(result) <- c('nu','gamma','F-measure')
result[which.max(result[,3]),]


############## Cost-sensitive Learning ################

# Class weighted SVM -----------------------------
# svm function using 'class.weight' option

wts = 500 / table(train$y)

set.seed(10)
cv.fit = tune(svm, y~., data=train, kernel='radial', 
              class.weights=wts,
              ranges=list(cost=c(0.001,0.01,0.1,1,5,10,100),
                          gamma=c(0.01,0.1,0.5,1,2,3,4)))

?tune
summary(cv.fit)

# SVM with the best parameters
best.fit = cv.fit$best.model

# Prediction for test data
yhat.te = predict(best.fit, test)
cm = table(true = test$y, predict=yhat.te)

misclass(cm)

fmeasure(cm)

################# Ensemble-based methods ###################
install.packages('ebmc')
library(ebmc)

################# SMOTE Boost#################
set.seed(10)

fit1 = sbo(y~., data=train, size=200, alg='cart', over=300)
# y should be encoded by (0,1); 0 large class, 1 small class
# size: # of boosting iterations
# alg: weak learner
# over: oversampling rate (multiple of 100 is only acceptible)

yhat.te = predict(fit1, test, type='class')
cm = table(true = test$y, predict=yhat.te)

misclass(cm)

fmeasure(cm)


################# SMOTE Bagging #################
set.seed(10)

fit2 = sbag(y~., data=train, size=300, alg='cart')
# y should be encoded by (0,1); 0 large class, 1 small class

yhat.te = predict(fit2, test, type='class')
cm = table(true = test$y, predict=yhat.te)

misclass(cm)

fmeasure(cm)






