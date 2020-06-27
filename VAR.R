############################
##########Libraries#########
############################
library(plotly)
library(tidyr)
library(dplyr)
library(ggplot2)

############################
########Functions###########
############################

readTheData <- function(filepath) {
  data <- read.csv(filepath, dec=',',header=TRUE, sep=";")
  return(data)
}

dates <- function(data){
  dates <- as.Date(data$DATE,"%d.%m.%Y")
  return(dates)
}

vec <- function(my.matrix){
	N <- nrow(my.matrix)
	M <- ncol(my.matrix)
	my.vec <- matrix(c(my.matrix), nrow=N*M)
	my.vec
}

transform_log_diff <- function(data, nDiff){

  n_col=ncol(data)
  n_row=nrow(data)
  log_diff =data.frame(matrix(0, nrow = n_row-nDiff, ncol = n_col))
  names<-colnames(data)
  for(i in 2:n_col){
    log_diff[,i] <- diff(log(data[,i]))
    colnames(log_diff)[i] <- paste("log_diff_", nDiff,"_", colnames(data)[i], sep="")
  }
  return(log_diff)
}

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

stability_check <- function(B_hat, constant,order){ 
  n_row <- nrow(B_hat)
  if(constant==1){
  stab <- (B_hat[,2:ncol(B_hat)])
  colnames(stab)<-seq(1,ncol(B_hat) -1 )
    for(i in 1:(order-1)){
      help<- data.frame(matrix(0, nrow = n_row, ncol = n_row *order))
      help[,(1+(i-1)*(n_row)):(n_row*i) ] <- diag(n_row)
      colnames(help)<-colnames(stab)
      stab <- rbind(stab,help)
    }
  }
  
  if(constant==0){  
  }
  return(sort(abs(eigen(stab)$values), decreasing = TRUE))

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

calc_wald_stat <- function(beta_hat, Z, sigma, C){
  ## function assumes test for non casuality
  inter_1 <- C %*% vec(beta_hat)
  inter_2 <- solve((Z %*% t(Z)))
  inter_2 <- kronecker(inter_2, sigma)
  inter_2 <- C %*% inter_2 %*% t(C)
  inter_2 <- solve(inter_2)
  return(t(inter_1) %*% inter_2 %*% inter_1)
}


calc_f_stat <- function(beta_hat , Z , sigma, C, order){
  ## function assumes test for non casuality
  ## and that number rows of C is rank(C)
  stat <- calc_wald_stat(beta_hat, Z, sigma, C) / nrow(C)
  pval <- 1-pf(stat, df1=nrow(C) , df2= nrow(beta_hat)*ncol(Z)-nrow(beta_hat)^2*order -nrow(beta_hat))
  return(list(stat,pval))
}

calc_phi_ma <- function(x, order, K,i){
  phi_test <- diag(K)
  phi_test <-  cbind(x[,(1+K*(1-1)):(K*1)], phi_test)
  if(i==1){
    return(phi_test)
  }
  if(i>1){
    for(i in 2: i){
      phi_interim <- matrix(0,K,K)
      for(j in 1:i){
        if(j>order) next
        phi_interim  <- phi_interim + as.matrix(x[,(1+K*(j-1)) : (K*j)] )%*% as.matrix(phi_test[,(1+K*(j-1)) : (K*j)])
      }  
      phi_test <- cbind(phi_interim, phi_test)
      
    }
    
  }
  return(phi_test)
  
}


IR <- function(h,mB){
  A.hat.array <- array(0, dim=c(2,2,h))
  A.hat.array[,,1:4] = array(mB4[,-1], dim=c(2,2,4))
  IR.array <- array(NA, dim=c(2,2,(h+1)))
  IR.array[,,1] = diag(2)

  for(i in 2:(h+1)){
    IR.here <- matrix(0, nrow=2, ncol=2)
  for(j in 1:(i-1)){
    IR.here <- IR.here + (IR.array[,,(i-j)] %*% A.hat.array[,,j])
  }
  IR.array[,,i]= IR.here
  }
  return(IR.array)
}

accIR <- function(h,IR.array){
  acc.IR.array <- array(NA, dim=c(2,2,(h+1)))
  acc.IR.array[,,1] = IR.array[,,1]

  for(i in 2:(h+1)){
    acc.IR.array[,,i] = apply(IR.array[,,1:i], c(1,2), sum)
  }
  return(acc.IR.array)
}

###################################
###################################
filepath = 'us_data.csv'

data = readTheData(filepath)
returns = transform_log_diff(data,1)
dates = dates(data)

###################################
#############Exercises#############
###################################

###################################
############### a) ################
###################################
x1<- plot_ly(data, x =~dates, y = ~IP, name = 'IP',type='scatter',mode='lines')%>%
   layout(showlegend=TRUE)

x2<- plot_ly(data, x =~dates, y = ~CPI, name = 'CPI', type='scatter',mode='lines')%>%
   layout(showlegend=TRUE)

x3<- plot_ly(data, x =~dates, y = ~FEDFUNDS, name='FED FUNDS', type='scatter', mode='lines')%>%
   layout(showlegend=TRUE)

x <- subplot(x1,x2,x3)

y4<- plot_ly(returns, x =~dates[2:300], y = ~log_diff_1_IP, name = "IP",type='scatter',mode='lines')%>%
   layout(showlegend=TRUE)

y5<- plot_ly(returns,x =~dates[2:300], y = ~log_diff_1_CPI, name = "CPI", type='scatter',mode='lines')%>%
   layout(showlegend=TRUE)

y <- subplot(y4,y5)


###################################
############### b) ################
###################################
dataSet1 = as.data.frame(cbind(returns[[2]],returns[[3]]))
dataSet2 = as.data.frame(cbind(returns[[2]],returns[[3]],returns[[4]]))

mB6 = var_estimation(dataSet1, 6, 1)[[1]]
vU6 = var_estimation(dataSet1, 6, 1)[[2]]
lambda = stability_check(mB6, 1, 6)

###################################
############### c) ################
###################################

df_crit <- data.frame(matrix(0,nrow=3, ncol=3 ))
for( i in 4:6){
  c_i_res <- var_estimation(dataSet1[(7-i):nrow(dataSet1),1:2], order =i , 1)
  U_hat <- as.data.frame(c_i_res[2])
  sigma <- calc_sigma_lse(U_hat, order = i)
  df_crit[i-3, 1] <- AIC(sigma,i,nrow(U_hat))
  df_crit[, 2] <- replicate(3,'&')
  df_crit[i-3, 3] <- HQ(sigma, i, nrow(U_hat))
  df_crit[, 4] <- replicate(3,'&')
  df_crit[i-3, 5] <- SIC(sigma, i, nrow(U_hat))
}

###################################
############### d) ################
###################################

##chosen model order is 4
mB4 = var_estimation(dataSet1, 4, 1)[[1]]
IR.array = IR(10,mB4)

IR.table = cbind(IR.array[1,1,], IR.array[1,2,], IR.array[2,1,], IR.array[2,2,])
colnames(IR.table) = c("X1", "X2","X3","X4")

g1<- plot_ly(as.data.frame(IR.table),y=~X1, name = 'IR: IP -> IP',type='scatter',mode='lines')%>%
   layout(showlegend=TRUE)

g2<- plot_ly(as.data.frame(IR.table),y=~X2, name = 'IR: CPI -> IP', type='scatter',mode='lines')%>%
   layout(showlegend=TRUE)

g3<- plot_ly(as.data.frame(IR.table),y=~X3, name='IR: IP -> CPI', type='scatter', mode='lines')%>%
   layout(showlegend=TRUE)

g4<- plot_ly(as.data.frame(IR.table),y=~X4, name='IR: CPI -> CPI', type='scatter', mode='lines')%>%
   layout(showlegend=TRUE)

g <- subplot(g1,g2,g3,g4,nrows=(2))
#htmlwidgets::saveWidget(g, "IR.html")



acc.IR.array = accIR(10,IR.array)
acc.IR.table = cbind(acc.IR.array[1,1,], acc.IR.array[1,2,], acc.IR.array[2,1,], acc.IR.array[2,2,])
colnames(acc.IR.table) = c("X1", "X2", "X3", "X4")


k1<- plot_ly(as.data.frame(acc.IR.table),y=~X1, name = 'acc IR: IP -> IP',type='scatter',mode='lines')%>%
   layout(showlegend=TRUE)

k2<- plot_ly(as.data.frame(acc.IR.table),y=~X2, name = 'acc IR: CPI -> IP', type='scatter',mode='lines')%>%
   layout(showlegend=TRUE)

k3<- plot_ly(as.data.frame(acc.IR.table),y=~X3, name='acc IR: IP -> CPI', type='scatter', mode='lines')%>%
   layout(showlegend=TRUE)

k4<- plot_ly(as.data.frame(acc.IR.table),y=~X4, name='acc IR: CPI -> CPI', type='scatter', mode='lines')%>%
   layout(showlegend=TRUE)

k <- subplot(k1,k2,k3,k4,nrows=(2))
htmlwidgets::saveWidget(k, "accIR.html")


###################################
############### e) ################
###################################

mB2 = var_estimation(as.data.frame(dataSet2), 2, 1)[[1]]
vU2 = var_estimation(as.data.frame(dataSet2), 2, 1)[[2]]
mZ = var_estimation(as.data.frame(dataSet2), 2, 1)[[3]]
SigmaU = calc_sigma_lse(vU2, 2)

###################################
############### f) ################
###################################

SigmaU = calc_sigma_lse(vU2, 2)
C_d = rbind(c(rep(0,9), 1, rep(0,11)),c(rep(0,10), 1, rep(0,10)),c(rep(0,18), 1, rep(0,2)),c(rep(0,19), 1, rep(0,1)))
Fstat = calc_f_stat(mB2 , mZ , SigmaU, C_d, 2)[[1]]
pval = calc_f_stat(mB2 , mZ , SigmaU, C_d, 2)[[2]]



