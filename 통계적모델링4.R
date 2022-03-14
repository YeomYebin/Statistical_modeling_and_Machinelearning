#####################################################
# 통계적 모델링과 머신러닝 실습 중간고사 MasterCode #
######################################################

## LSE = MLE : 선형회귀분석
## 조건 : 1) Error's independent, 2) Variance of Yi doesn't depend on its mean value 3) Normal조건..?

library(tidyverse)
library(datasets)
?attenu
dim(attenu)

fit1 = lm(accel ~ mag + dist, data=attenu)
## par(mforw = c(2, 2)) ## Plot 그리드 설정
plot(fit1, 1) ## Residuals vs Fitted : 선형성, 오차의 등분산성 확인 가능
plot(fit1, 2) ## Nomaal QQ plot : 정규성 확인 가능, 대표적인 비모수적 방법이자 시각적 방법.
plot(fit1, 3) ## Scale - Location : 첫번째와 비슷하게 확인 가능, 선형성과 오차의 등분산성 확인 가능, 잔차에 절댓값이 씌워진 형태
plot(fit1, 4) ## Residuals vs Leverage : Outliar 확인 가능. Cook's Distance 0.5의 경계

fit1$coefficients ## 계수 : Beta
fit1$residuals ## Residuals
fit1$fitted.values ## Yhat : 예측값

## NonLinear 선형회귀분석 ##

## GAM -> ANOVA -> Function세우기 -> Object Function세우기 : RSS -> Gradient vector of object function -> Optimization하기

## GAM : Generalized addictive model
library(gam)
fit2 = gam(accel ~ s(mag,5) + s(dist,5), data=attenu)
par(mfrow=c(1,2))
plot(fit2)

fit2_1 = gam(accel ~ mag + s(dist,5),data=attenu)
plot(fit2_1)

anova(fit2,fit2_1) ##유의한지 안한지 확인

## Nonlinear Model with Constant variance
# Y = accel, X1 = mag, X2 = dist.
# Nonlinear model: Y = beta1 + beta2*X1 + beta3*X1^2 + beta4*exp(-beta5*X2)

## (1) Function세우기
f = function(beta,X)
{
  X1 = X[,1]; X2 = X[,2] ## X변수 지정
  beta[1] + beta[2]*X1 + beta[3]*X1^2 + beta[4]*exp(-beta[5]*X2)
}

## (2) Objective function: RSS 세우기
RSS = function(beta,Y,X) sum((Y-f(beta,X))^2)

## (3) Gradient vector of the objective function
grv = function(beta,Y,X)
{
  X1 = X[,1]; X2 = X[,2]
  R = Y - f(beta,X)
  c(-2*sum(R), -2*sum(R*X1), -2*sum(R*X1^2), -2*sum(R*exp(-beta[5]*X2)), 2*beta[4]*sum(R*X2*exp(-beta[5]*X2))) ## Beta(i)에 대한 편미분
}

# Optimization
X = cbind(attenu$mag,attenu$dist) ## X변수 지정
colnames(X) = c('mag', 'dist')
Y = attenu$accel
ml1 = optim(rep(0.1,5), RSS, gr=grv, method='BFGS', X=X, Y=Y) ## optim(초기 Beta값, Object Function, Gradient, Lower, Upper)
ml1

ml1$convergence ## 0이면 수렴함

beta.hat = ml1$par
beta.hat

# Fitted value : Yhat
Yhat = f(beta.hat,X)

# Residual plot
r = Y - Yhat
par(mfrow=c(1,1))
plot(Yhat,r,ylim=c(-0.5,0.5)) ##Y의 경계 설정 바꾸기
lines(c(-10,10),c(0,0),col='red') ## 빨간선 그리기 ## 첫번째 플랏과 비교
# Linearly increasing variance pattern.


## 1. Non Constant Variance ##
## Weighted Least squares

### E ~ N(0, sigma_i^2)
### Y ~ MVN(F, SIGMA)
### => WLSE

### 첫번째 플랏을 보고 등분산성 확인 후 Variance Function 적용

library(matrixcalc) # To check whether a matrix is singular or not
X = cbind(attenu$mag,attenu$dist) ## X변수 지정
colnames(X) = c('mag', 'dist')
Y = attenu$accel

# function 설정
f = function(beta,X)
{
  X1 = X[,1]; X2 = X[,2]  
  beta[1] + beta[2]*X1 + beta[3]*X1^2 + beta[4]*exp(-beta[5]*X2)
}

# Objective function for mean function: Genearalized least square method.
obj.mean = function(beta,Y,X,S) t(Y-f(beta,X)) %*% solve(S) %*% (Y-f(beta,X))
# S: Covariance matrix

# Gradient vector of the objective function
gr.mean = function(beta,Y,X,S)
{
  sigma2 = diag(S) ##sigma^2 설정!!! ##
  X1 = X[,1]; X2 = X[,2]
  R = Y - f(beta,X)
  c(-2*sum(R/sigma2), -2*sum(R*X1/sigma2), -2*sum(R*X1^2/sigma2), 
    -2*sum(R*exp(-beta[5]*X2)/sigma2), 
    2*beta[4]*sum(R*X2*exp(-beta[5]*X2)/sigma2))
}

# Linear variance function: |r| = gam1 + gam2*Yhat. ## Variance Fucntion 지정
# For linear variance function, we can consider absolute residuals,
# instead of squared residuals.
# gam.hat = (Z^T W Z)^(-1) Z^T W |r|. ##WLSE

beta.new = ml1$par      # initial parameter.
W = diag(rep(1,length(Y))) ##초기값 설정 ##SIGMA^(-1)
mdif = 100000

while(mdif > 0.000001)
{
  Yhat = f(beta.new,X)
  r = Y - Yhat
  Z = cbind(1,Yhat) ## Variance Function 식 여기서 설정
  gam.hat = solve(t(Z) %*% W %*% Z) %*% t(Z) %*% W %*% abs(r) ##WLSE : Variance의 gamma 먼저 구하기
  sigma = Z %*% gam.hat ## Sigma구하기!!!!
  S = diag(as.vector(sigma^2)) ## SIGMA
  
  ## Inverse matrix구하기!
  if (is.non.singular.matrix(S)) W = solve(S)
  else W = solve(S + 0.000000001*diag(rep(1,nrow(S))))
  
  ml2 = optim(beta.new, obj.mean, gr=gr.mean, method='BFGS', Y=Y, X=X, S=S) ## S : SIGMA
  beta.old = beta.new 
  beta.new = ml2$par ## Beta 업데이트
  mdif = max(abs(beta.new - beta.old)) ##수렴확인
}

beta.new ##최종 Beta

Yhat = f(beta.new, X)
sigma = Z %*% gam.hat ## 최종 SIGMA
r = (Y - Yhat)/sigma

# Residual plot
plot(Yhat,r, ylim=c(-4,4))
lines(c(0,10), c(0,0),col='red') ## 등분산이 나아졌는지 확인!!

gam.hat


## Method2 : Variance Function ; lmvar 패키지 이용 ##
##### Linear regression with nonconstant variance #####
# lmvar: 
# Linear mean function
# Linear variance function: log(sigma) = X*beta
#install.packages('lmvar')
library(lmvar)

X = cbind(attenu$mag,attenu$dist)
colnames(X) = c('mag', 'dist')

X_s = cbind(attenu$mag, attenu$dist)
colnames(X_s) = c('mag', 'dist')

fit3 = lmvar(attenu$accel, X, X_s)
summary(fit3)

?lmvar

ms = predict(fit3, X_mu=X, X_sigma=X_s)
r1 = (Y - ms[,1])/ms[,2]
plot(ms[,1], r1)
lines(c(-10,10),c(0,0),col='red')

################################################################################################
## 2. Error's are not independant
## time correlation

tsdat = read.table('C:/Users/User/Desktop/통계적모델링과머신러닝/tsdat.txt',header=T)

fit = lm(y ~ x, data=tsdat)
summary(fit)

par(mfrow=c(2,2))
plot(fit)

# Durbin-Watson test : 독립성 검정
## H0 : 귀무가설(H0)는 잔차들 사이에 자기상관관계가 없다, 즉 독립적이다.
# install.packages('lmtest')
library(lmtest)

dwtest(fit)

# Check ACF & PACF
# install.packages('astsa')
library(astsa)

# AR(p): ACF: Exponentially decreasing; PACF: Non-zero values at first p lags.
# MA(q): ACF: Non-zero values at first q lags; PACF: Exponentially decreasing.
# ARMA(p,q): ACF: Similar to ACF of AR(p); PACF: Similar to PACF of MA(q).

acf2(residuals(fit))
acf(residuals(fit))
pacf(residuals(fit))

# ar1 = arima(residuals(fit), c(1,0,0))  #AR(1)
#ar1$coef

library(forecast)
#auto.arima(residuals(fit))

ar1 = sarima(residuals(fit), 1,0,0, no.constant=T)  
ar1$fit$coef
ar1$fit$sigma2
ar1$fit

X = cbind(1,tsdat$x)
Y = tsdat$y
n = length(Y)
S = diag(rep(1,n))# initial covariance matrix

mdif = 1000000
beta.old = rep(100000,2)

while(mdif > 0.0000001)
{
  beta.new = as.vector(solve(t(X) %*% solve(S) %*% X) %*% t(X) %*% solve(S) %*% Y)
  r = as.vector(Y - (X %*% beta.new))
  ar1 = sarima (r, 1,0,0, no.constant=T, details=F)
  alpha = ar1$fit$coef
  sigma2 = ar1$fit$sigma2
  
  mdif = max(abs(beta.new - beta.old))
  beta.old = beta.new
  
  # Construct covariance matrix
  S = matrix(nrow=n,ncol=n)
  for (i in 1:n)
  {
    for (j in 1:n)
    {
      if (i == j) S[i,j] = 1
      if (i != j) S[i,j] = alpha^(abs(i-j))
    }
  }
  S = (sigma2 / (1-alpha^2)) * S
}

round(beta.new,4)

## Method2 : conditional pdf 사용 ##
# MLE: Product of conditional distribution (Approximation)
# Y_t | Y_t-1 ~ N(X_t*beta + alpha*epsilon_t-1, sigma^2)

fit = lm(y ~ x, data=tsdat)

Yt = tsdat$y[2:n]
Xt = tsdat$x[2:n]
et = residuals(fit)[1:(n-1)]
mdif = 10000
b.old = rep(0,3)

while(mdif > 0.0000001)
{
  fit.temp = lm(Yt ~ Xt + et)  ## Y = B0 + B1X_t + alpha*E(t-1) + eta_t , eta_t ~ N(0, sigma2)
  b.new = fit.temp$coefficient
  mdif = max(abs(b.new[1:2] - b.old[1:2]))
  
  et = (Y - X %*% b.new[1:2])[1:(n-1)]
  b.old = b.new
}

round(b.new,4)

# Built-in function 
# cochrane.orcutt => f: linear model, error: AR(p) process.
# install.packages("orcutt")
library(orcutt)

fit = lm(y ~ x, data=tsdat)
cochrane.orcutt(fit)

## Regression with MA1
X = cbind(1,tsdat$x)
Y = tsdat$y
n = length(Y)
S = diag(rep(1,n))# initial covariance matrix

mdif = 1000000
beta.old = rep(100000,2)

while(mdif > 0.0000001)
{
  beta.new = as.vector(solve(t(X) %*% solve(S) %*% X) %*% t(X) %*% solve(S) %*% Y)
  r = as.vector(Y - (X %*% beta.new))
  ma1 = sarima (r, 0,0,1, no.constant=T, details=F)
  theta = ma1$fit$coef
  sigma2 = ma1$fit$sigma2
  
  mdif = max(abs(beta.new - beta.old))
  beta.old = beta.new
  
  # Construct covariance matrix
  S = matrix(nrow=n,ncol=n)
  for (i in 1:n)
  {
    for (j in 1:n)
    {
      if (i == j) S[i,j] = (theta + 1)^2
      if (i != j) S[i,j] = 0
    }
  }
  S = (sigma2) * S
}

round(beta.new,4)


##################################################################################################
## 2. Error's are not independant
## Spatial correlation

library(spdep)
library(spatialreg)

data(oldcol)

?COL.OLD
# 'COL.nb' has the neighbors list.

# 2-D Coordinates of observations
crds = cbind(COL.OLD$X, COL.OLD$Y)  

# Compute the maximum distance
mdist = sqrt(sum(diff(apply(crds, 2, range))^2))   ## 1 : 행, 2 : 열 ## range(X) : c(min(X), max(X))

# All obs. between 0 and mdist are identified as neighborhoods.
dnb = dnearneigh(crds, 0, mdist)
?dnearneigh

# Compute Euclidean distance between obs.
dists = nbdists(dnb, crds)

# Compute Power distance weight d^(-2)
glst = lapply(dists, function(d) d^(-2))

## K-Nearest Neightbor Weights :
'''K_near = function(d, k){
  w =c()
  for ( i in 1:length(d) ){
    if (d[i] <= sort(d)[k]){
      w[i] = 1
    } else w[i] = 0
  }
  return(w)
}

glst = lapply(dists, K_near, k)'''


## Radial distance weights :
'''Radial_d = function(d, d1){
  w = c()
  for (i in 1:length(d)){
    if(d[i]<= d1){
      w[i] = 1
    }else w[i] = 0
  }
  return(w)
}
glst = lapply(dists, Radial_d, d1)'''

## Exponential distance weights :
'''alpha = 0.5
glst = lapply(dists, function(d) exp(-alpha*d))'''

## Double-Power distance weights :
'''k = 2
d1 = mdist/2
glst = lapply(dists, function(d) {(1-((d/d1)^k))^k} )'''


# Construct weight matrix with normalization 
# style='C': global normalization; 'W': row normalization
lw = nb2listw(dnb, glist=glst, style='C')

# Spatial Autoregressive Model
fit = lagsarlm(CRIME ~ HOVAL + INC, data = COL.OLD, listw=lw )
?lagsarlm
summary(fit)

# Fitted values
predict(fit)

##########################################################################3
## 3. GLM
library(ordinal)

### (1) : Logistic Regression  - binomial
library(survival)
data(colon)
fit = glm(status ~ sex + age + obstruct + perfor, data = colon, family =binomial)
summary(fit)
prob = predict(fit, newdata = colon, type = 'response')
Yhat = ifelse(prob > 0.5 , 1 , 0)
summary(fit) ## 오즈비 이용해서 해석 -> 범주교안 참고

## (2) 다범주 로짓 모형
data(iris)
iris$Species
train_data = iris
train_data$Species <- relevel(train_data$Species, "virginica")
library(nnet)
mlogit <- multinom(Species ~. , data=train_data)
summary(mlogit)
fit <- fitted(mlogit) ## 각 클래스에 속할 확률
pred_Species <- predict(mlogit, data=test_data, type='probs') ## 각 클래스에 속할 확률 ## type = 'class'로도 가능

## (3) Cumulative Regresstion
library(ordinal)
fit = clm(rating ~ temp + contact, data=wine, link = 'logit')
summary(fit)

## (4) Poisson regression model ##

# For data
# install.packages('lme4')
library(lme4)
data(grouseticks)
head(grouseticks)
par(mfrow = c(1, 1))
hist(grouseticks$TICKS,breaks=0:90)

fit = glm(TICKS ~ HEIGHT*YEAR, data = grouseticks, family=poisson)
summary(fit)

## 과대산포 검정
library(AER)
dispersiontest(fit) ##평균과 분산이 같다는 귀무가설을 기각합니다.


## (5) Negative binomial regression model ###
library(MASS)

fit1 = glm.nb(TICKS ~ HEIGHT*YEAR, data = grouseticks, link=log)
summary(fit1)

####################################################################################
########## Proportional hazard model ##########
library(survival)

# For data
# install.packages('carData')
library(carData)

?Rossi

Rossi$week
Rossi$arrest

Surv(Rossi$week, Rossi$arrest)

fit = coxph(Surv(week,arrest) ~ fin + age+ race + wexp + mar + prio, 
            data=Rossi)
summary(fit)

# Estimated survival function
plot(survfit(fit),ylim=c(0.6,1),xlab="Weeks", ylab="Prop.of Not Rearrested")


# Estimated survival functions for financial aid
Rossi.fin = with(Rossi, data.frame(fin=c(0, 1), age=rep(mean(age), 2), 
                                   race=rep(mean(race=='other'),2), 
                                   wexp=rep(mean(wexp=="yes"),2), 
                                   mar=rep(mean(mar=="not married"),2),
                                   prio=rep(mean(prio),2)))

plot(survfit(fit,newdata=Rossi.fin), conf.int=TRUE,
     lty=c(1, 2), ylim=c(0.6, 1), col=c('red','blue'), 
     xlab="Weeks", ylab="Prop. of Not Rearrested")

legend("bottomleft", legend=c("fin = no","fin = yes"), 
       lty=c(1 ,2),col=c('red','blue'))
