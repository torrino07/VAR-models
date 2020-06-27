filepath='md_eu.csv'
df <- read.csv(filepath, header=TRUE, sep=",")
df$m <- log(df$M)
df$p <- log(df$P)
df$y <- log(df$Y)
df$whh <- log(df$WHH)
df$`(m-p)` <- df$m - df$p
df$diffwhh <- NA
df$diffwhh[2:nrow(df)] <- diff(df$whh)

library(vars)
library(urca)
library(tsDyn)
library(plotly)
library(xtable)

plot_ly(x=as.Date(df$X,"%Y-%M") ,y=df$`(m-p)` ,name="Real money stock(m-p)" ,type="scatter" ,mode="lines") %>% layout(showlegend=TRUE)
plot_ly(x=as.Date(df$X,"%Y-%M") ,y=df$y ,name="Natural log of RGDP" ,type="scatter" ,mode="lines") %>% layout(showlegend=TRUE)
plot_ly(x=as.Date(df$X,"%Y-%M") ,y=df$RS ,name="Nominal interest rate" ,type="scatter" ,mode="lines") %>% layout(showlegend=TRUE)
plot_ly(x=as.Date(df$X[2:nrow(df)],"%Y-%M") ,y=df$diffwhh[2:nrow(df)] ,name="Growth rate housing wealth" ,type="scatter" ,mode="lines") %>% layout(showlegend=TRUE)

var_estimation <- function(data, order , constant){
  Y <-  t(na.omit(data[order+1:nrow(data),]))
  HelpMatrix <- t(data)
  n_col<- ncol(data)
  n_row<- nrow(data)
  if(constant==1){
    Z<- data.frame(matrix(1, nrow = order*n_col + 1, ncol = n_row -order))
    for(i in 1:ncol(Z)){
      for(j in 1:order){
        Z[(2+ n_col * (j-1)):(1 +(n_col*j) ) , i] <- HelpMatrix[,i+(order-j)]
      }
      
    }
  }
  if(constant==0){
    Z<- data.frame(matrix(0, nrow = order*n_col, ncol = n_row -order))
    for(i in 1:ncol(Z)){
      for(j in 1:order*n_col){
        Z[(1+n_col * (j-1)): (n_col*j)  , i] <- HelpMatrix[,i+(order-j)]
      }
    }
  }
  Z <-as.matrix(Z)
  Y <- as.matrix(Y)
  B_hat <- Y %*% t(Z) %*% solve(Z %*% t(Z))
  U_hat <- Y - B_hat %*% Z
  return(list(B_hat, U_hat,Z))
}

calc_sigma_lse <- function(U_hat, order){

  sigma_u_hat <- 1/(ncol(U_hat) - order) * crossprod(t(U_hat))
 
}

AIC <- function(sigma ,order , sample_size ){
  return( log(det(sigma) ) + 2/sample_size * (order* ncol(sigma)^2 + ncol(sigma)))
 
}

HQ <- function(sigma,order, sample_size){
  return ( log(det(sigma)) + (2*log(log(sample_size)))/sample_size * (ncol(sigma)^2*order+ncol(sigma)))
  
}

SIC <- function(sigma,order, sample_size){
  return(log(det(sigma))  + log(sample_size)/sample_size * (ncol(sigma)^2*order+ncol(sigma)))
}
##########################################################
##########################################################
## estimate lag of order
# can be done by estimating var(p) 
df_crit <- data.frame(matrix(0,nrow=5, ncol=3 ))
colnames(df_crit) <- c("AIC", "HC", "SIC")
rownames(df_crit) <- c("p=2", "p=3","p=4", "p=5", "p=6")
##df$X[2:111]

data <- df[2:111, c("(m-p)","y" , "RS")]

for( i in 2:6){
  c_i_res <- var_estimation(data[(6-i):nrow(data),], order =i , 1)
  U_hat <- as.data.frame(c_i_res[2])
  sigma <- calc_sigma_lse(U_hat, order = i)
  df_crit[i-1, 1] <- AIC(sigma,i,nrow(U_hat))
  df_crit[i-1, 2] <- HQ(sigma, i, nrow(U_hat))
  df_crit[i-1, 3] <- SIC(sigma, i, nrow(U_hat))
}

VECM_model = VECM(data,lag=1,r=1 ,include= c("const") , LRinclude = c("trend") , estim ="ML")
summary(VECM_model)
rank.test(VECM_model, type="trace")
#Rank selected: 1 (first trace test with pval above 5 %: 37.8 %)


##########################################################
##########################################################
df_crit <- data.frame(matrix(0,nrow=5, ncol=3 ))
colnames(df_crit) <- c("AIC", "HC", "SIC")
rownames(df_crit) <- c("p=2", "p=3","p=4", "p=5", "p=6")
data_c <- df[2:nrow(df) ,c("(m-p)","y" , "RS")]
for( i in 2:6){
  c_i_res <- var_estimation(data_c[(6-i):nrow(data_c),], order =i , 1)
  U_hat <- as.data.frame(c_i_res[2])
  sigma <- calc_sigma_lse(U_hat, order = i)
  df_crit[i-1, 1] <- AIC(sigma,i,nrow(U_hat))
  df_crit[i-1, 2] <- HQ(sigma, i, nrow(U_hat))
  df_crit[i-1, 3] <- SIC(sigma, i, nrow(U_hat))
}


VECM_model = VECM(data_c,lag=1,r=1 ,include= c("const") , LRinclude = c("trend") , estim ="ML")
rank.test(VECM_model, type="trace")
#Rank selected: 1 (first trace test with pval above 5 %: 41.6 %)

##########################################################
##########################################################

df_crit <- data.frame(matrix(0,nrow=5, ncol=3 ))
colnames(df_crit) <- c("AIC", "HC", "SIC")
rownames(df_crit) <- c("p=2", "p=3","p=4", "p=5", "p=6")
data <- df[2:111, c("(m-p)","y" , "RS", "diffwhh")]
for( i in 2:6){
  c_i_res <- var_estimation(data[(6-i):nrow(data),], order =i , 1)
  U_hat <- as.data.frame(c_i_res[2])
  sigma <- calc_sigma_lse(U_hat, order = i)
  df_crit[i-1, 1] <- AIC(sigma,i,nrow(U_hat))
  df_crit[i-1, 2] <- HQ(sigma, i, nrow(U_hat))
  df_crit[i-1, 3] <- SIC(sigma, i, nrow(U_hat))
}


VECM_model = VECM(data,lag=1,r=1 ,include= c("const") , LRinclude = c("trend") , estim ="ML")
rank.test(VECM_model, type="trace")
#Rank selected: 1 (first trace test with pval above 5 %: 14.6 %)


##########################################################
##########################################################
df_crit <- data.frame(matrix(0,nrow=5, ncol=3 ))
colnames(df_crit) <- c("AIC", "HC", "SIC")
rownames(df_crit) <- c("p=2", "p=3","p=4", "p=5", "p=6")
data <- df[2:nrow(df) , c("(m-p)","y" , "RS", "diffwhh")]
for( i in 2:6){
  c_i_res <- var_estimation(data[(6-i):nrow(data),], order =i , 1)
  U_hat <- as.data.frame(c_i_res[2])
  sigma <- calc_sigma_lse(U_hat, order = i)
  df_crit[i-1, 1] <- AIC(sigma,i,nrow(U_hat))
  df_crit[i-1, 2] <- HQ(sigma, i, nrow(U_hat))
  df_crit[i-1, 3] <- SIC(sigma, i, nrow(U_hat))
}


VECM_model = VECM(data,lag=1,r=1 ,include= c("const") , LRinclude = c("trend") , estim ="ML")
rank.test(VECM_model, type="trace")
#Rank selected: 1 (first trace test with pval above 5 %: 32.9 %)