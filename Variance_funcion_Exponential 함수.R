library(tidyverse)
library(datasets)
?attenu
dim(attenu)

## Variance Funcrtion이 지수함수 : Optim 2번 돌리기 ##

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
obj.mean = function(beta,Y,X,S){
  if (is.non.singular.matrix(S)) W = solve(S)
  else W = solve(S + 0.000000001*diag(rep(1,nrow(S))))
  t(Y-f(beta,X)) %*% W %*% (Y-f(beta,X))}
# S: Covariance matrix

# Gradient vector of the objective function
gr.mean = function(beta,Y,X,S)
{
  sigma2 = diag(S) ##sigma^2 설정!!! ##
  X1 = X[,1]; X2 = X[,2]
  R = Y - f(beta,X)
  c(-2*sum(R/sigma2), -2*sum(R*X1/sigma2), -2*sum(R*X1^2/sigma2), 
    -2*sum(R*exp(-beta[5]*X2)/sigma2), 2*beta[4]*sum(R*X2*exp(-beta[5]*X2)/sigma2))
}

f.var = function(gamma,Z)
{
  #Z = Z[,1]
  gamma[1] + gamma[2]*exp(gamma[3]*Z)
}

obj.var = function(gamma , ri, Z, S){
  if (is.non.singular.matrix(S)) W = solve(S)
  else W = solve(S + 0.000000001*diag(rep(1,nrow(S))))
  t(abs(ri) - f.var(gamma,Z)) %*% W %*% (abs(ri) - f.var(gamma,Z))
} 

## (3) Gradient vector of the objective function
grv.var = function(gamma, ri , Z, S)
{
  #Z = Z[,1]
  sigma2 = diag(S)
  R.var = abs(ri) - f.var(gamma, Z)
  c(-2*sum(R.var/sigma2), -2*sum(R.var*exp(gamma[3]*Z)/sigma2), -2*gamma[2]*sum(R.var*Z*exp(gamma[3]*Z)/sigma2) ) ## gamma(i)에 대한 편미분
}


# Linear variance function: |ri| = gam1*exp(gam2*Yhat). ## Variance Fucntion 지정

beta.new = rep(0.1, 5)
gam.new = rep(0.1, 3) # initial parameter.
W = diag(rep(1,length(Y)))##초기값 설정 ##SIGMA^(-1)
S = diag(rep(1, length(Y)))
mdif = 100000
iter = 0

while(mdif > 0.00001 | iter <= 1000){
  Yhat = f(beta.new, X)
  r = Y - Yhat
  Z = Yhat
  ri = r
  ####
  ## (2) Objective function: RSS 세우기
  #RSS.var = function(gamma, ri_2 , Z) sum((ri_2 - f.var(gamma, Z))^2)
  
  
  # Optimization
   #cbind(Yhat) ## Z변수 지정
  #colnames(Z) = c('Yhat')
  ml.var = optim(gam.new , obj.var, gr=grv.var , method='BFGS', Z = Z , ri = ri , S = S) ## optim(초기 Beta값, Object Function, Gradient, Lower, Upper)
  gam.old = gam.new
  gam.new = ml.var$par
  ###
  
  sigma = f.var(gam.new, Z) ## sigma
  S = diag(as.vector(sigma)^2) ## SIGMA
  
  ## Inverse matrix구하기!
  #if (is.non.singular.matrix(S)) W = solve(S)
  #else W = solve(S + 0.000000001*diag(rep(1,nrow(S))))
  
  ml2 = optim(beta.new, obj.mean, gr=gr.mean, method='BFGS', Y=Y, X=X, S=S) ## S : SIGMA
  beta.old = beta.new 
  beta.new = ml2$par ## Beta 업데이트
  
  mdif = max(abs(beta.new - beta.old)) ##수렴확인
  iter  = iter + 1}

iter
beta.new ##최종 Beta

Yhat.final = f(beta.new, X)
sigma.final = f.var(gam.new, Yhat.final) ## 최종 SIGMA
r.final = (Y - Yhat.final)/sigma.final

# Residual plot
plot(Yhat.final,r.final)
lines(c(0,10), c(0,0),col='red') ## 등분산이 나아졌는지 확인!!

gam.new




######################################
## Variance Function ri^2 사용
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
obj.mean = function(beta,Y,X,S){
  if (is.non.singular.matrix(S)) W = solve(S)
  else W = solve(S + 0.000000001*diag(rep(1,nrow(S))))
  t(Y-f(beta,X)) %*% W %*% (Y-f(beta,X))
}
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

f.var = function(gamma,Z)
{
  #Z = Z[,1]
  gamma[1] + gamma[2]*exp(gamma[3]*Z)
}

obj.var = function(gamma , ri, Z, S_2){  
  if (is.non.singular.matrix(S_2)) W_2 = solve(S_2)
  else W_2 = solve(S_2 + 0.000000001*diag(rep(1,nrow(S_2))))
  t(ri^2 - f.var(gamma,Z)^2) %*% W_2 %*% (ri^2 - f.var(gamma,Z)^2)
  }

## (3) Gradient vector of the objective function
grv.var = function(gamma, ri , Z, S_2)
{
  #Z = Z[,1]
  sigma4 = diag(S_2)
  R.var = ri^2 - f.var(gamma, Z)^2
  c(-2*sum(R.var/sigma4), -2*sum(R.var*exp(gamma[3]*Z)/sigma4), -2*gamma[2]*sum(R.var*Z*exp(gamma[3]*Z)/sigma4) ) ## gamma(i)에 대한 편미분
}


beta.new = rep(0.1, 5)
gam.new = rep(0.01, 3) # initial parameter.
S_2 = diag(rep(1,length(Y)))##초기값 설정 ##SIGMA^(-1)
S = diag(rep(1, length(Y)))
mdif = 100000
iter = 0

while(mdif > 0.00001 | iter <= 1000){
  Yhat = f(beta.new, X)
  r = Y - Yhat
  Z = Yhat
  ri = r
  ####
  ## (2) Objective function: RSS 세우기
  #RSS.var = function(gamma, ri_2 , Z) sum((ri_2 - f.var(gamma, Z))^2)
  
  
  # Optimization
  #cbind(Yhat) ## Z변수 지정
  #colnames(Z) = c('Yhat')
  ml.var = optim(gam.new, obj.var, gr= grv.var , method='BFGS', Z = Z , ri = ri , S_2 = S_2) ## optim(초기 Beta값, Object Function, Gradient, Lower, Upper)
  gam.old = gam.new
  gam.new = ml.var$par
  ###
  
  sigma_2 = f.var(gam.new, Z) ## sigma
  S = diag(as.vector(sigma_2)) ## SIGMA
  S_2 = diag(as.vector(sigma_2^2))
  
  ## Inverse matrix구하기!
  #if (is.non.singular.matrix(S)) W = solve(S)
  #else W = solve(S + 0.000000001*diag(rep(1,nrow(S))))
  
  ml2 = optim(beta.new, obj.mean, gr=gr.mean, method='BFGS', Y=Y, X=X, S=S) ## S : SIGMA
  beta.old = beta.new 
  beta.new = ml2$par ## Beta 업데이트
  
  mdif = max(abs(beta.new - beta.old)) ##수렴확인
  iter  = iter + 1}

iter
beta.new ##최종 Beta

Yhat.final = f(beta.new, X)
sigma.final = f.var(gam.new, Yhat.final) ## 최종 SIGMA^2
r.final = (Y - Yhat.final)/sqrt(sigma.final)

# Residual plot
plot(Yhat.final, r.final )
lines(c(0,10), c(0,0),col='red') ## 등분산이 나아졌는지 확인!!

gam.new


