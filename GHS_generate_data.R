######
## Generate data for precision matrix estimation
######

rm(list=ls())

if (!require(mvtnorm)) install.packages('mvtnorm')
if (!require(moments)) install.packages('moments')
if (!require(Matrix)) install.packages('Matrix')
library(mvtnorm)
library(moments)
library(Matrix)

#Set number of dimension and number of observations

p=100; n=120
# p=200; n=120

for(precision_str in c(1,2,3,4))
{
  if(precision_str == 1)
  {
    ############################################
    #Case One: random structure 
    set.seed(2016)
    #generate off-diagonal elements in \Omega
    a = rep(NA,p*(p-1)/2)
    for (i in 1:(p*(p-1)/2)) {
      if(p == 100)
      {
        if (runif(1)>0.01)
        {
          a[i] = 0 #elements nonzero with prob=0.01, for p=100
        }else{
          a[i] = -runif(1,0.2,1)
        }    
      }else{
        if (runif(1)>0.002)
        {
          a[i] = 0 #elements nonzero with prob=0.002, for p=200
        }else{
          a[i] = -runif(1,0.2,1)
        }    
      }
    }
    A = matrix(NA,nrow=p,ncol=p)
    A[upper.tri(A)] = a; diag(A) = 1
    sigma_inv = forceSymmetric(A); sigma_inv = as.matrix(sigma_inv)
    #if \Omega is not positive-definite, generate again till it is
    while (min(eigen(sigma_inv)$values) < 0.01) {
      for (i in 1:(p*(p-1)/2)) {
        if(p == 100)
        {
          if (runif(1)>0.01)
          {
            a[i] = 0 #elements nonzero with prob=0.01, for p=100
          }else{
            a[i] = -runif(1,0.2,1)
          }    
        }else{
          if (runif(1)>0.002)
          {
            a[i] = 0 #elements nonzero with prob=0.002, for p=200
          }else{
            a[i] = -runif(1,0.2,1)
          }    
        }
      }
      A = matrix(NA,nrow=p,ncol=p)
      A[upper.tri(A)] = a; diag(A) = 1
      sigma_inv = forceSymmetric(A); sigma_inv = as.matrix(sigma_inv)
    }
    eigen(sigma_inv)$values
    sigma = solve(sigma_inv); sigma = as.matrix(sigma)
    file_name = paste0("GHS_sim_p",p,"randomn_sigmainv.csv")
    write.table(sigma_inv,file=file_name,sep=",",row.names=FALSE,col.names=FALSE)
    
  }else if(precision_str == 2){
    ############################################
    #Case Two: hubs structure
    sigma_inv = diag(1,nrow=p,ncol=p)
    alpha = 0.25      #for magnitude of nonzero partial covariance
    num_hub = p/10    #number of hubs
    for (i in 1:num_hub) {
      sigma_inv[(1+p/num_hub*(i-1)),(2+p/num_hub*(i-1)):(p/num_hub+p/num_hub*(i-1))] = alpha
      sigma_inv[(2+p/num_hub*(i-1)):(p/num_hub+p/num_hub*(i-1)),(1+p/num_hub*(i-1))] = alpha
    }
    eigen(sigma_inv)$values
    sigma = solve(sigma_inv); sigma = as.matrix(sigma)
    file_name = paste0("GHS_sim_p",p,"hubs_sigmainv.csv")
    write.table(sigma_inv,file=file_name,sep=",",row.names=FALSE,col.names=FALSE)
    
  }else if(precision_str == 3){
    ############################################
    #Case Three: cliques positive structure
    k=3    #number of features in a group; ten features make a group
    alpha = rep(0.75,k*(k-1)/2)    # for cliques positive
    set.seed(2016)
    sigma_inv = matrix(0,nrow=p,ncol=p)
    for (i in 1:(p/10)) {
      rind = sample((1+10*(i-1)):(10+10*(i-1)),size=k)
      sigma_inv[rind,rind][upper.tri(sigma_inv[rind,rind],diag=FALSE)] = alpha
      sigma_inv[rind,rind] = sigma_inv[rind,rind]+t(sigma_inv[rind,rind])
    }
    diag(sigma_inv) = 1
    eigen(sigma_inv)$values
    sigma = solve(sigma_inv); sigma = as.matrix(sigma)

    file_name = paste0("GHS_sim_p",p,"cliquespos75_sigmainv.csv")
    write.table(sigma_inv,file=file_name,sep=",",row.names=FALSE,col.names=FALSE)
    
  }else{
    ############################################
    #Case Four: cliques negative structure
    k=3    #number of features in a group; ten features make a group
    alpha = rep(-0.45,k*(k-1)/2)   # for cliques negative
    set.seed(2016)
    sigma_inv = matrix(0,nrow=p,ncol=p)
    for (i in 1:(p/10)) {
      rind = sample((1+10*(i-1)):(10+10*(i-1)),size=k)
      sigma_inv[rind,rind][upper.tri(sigma_inv[rind,rind],diag=FALSE)] = alpha
      sigma_inv[rind,rind] = sigma_inv[rind,rind]+t(sigma_inv[rind,rind])
    }
    diag(sigma_inv) = 1
    eigen(sigma_inv)$values
    sigma = solve(sigma_inv); sigma = as.matrix(sigma)
    file_name = paste0("GHS_sim_p",p,"cliquesneg45_sigmainv.csv")
    write.table(sigma_inv,file=file_name,sep=",",row.names=FALSE,col.names=FALSE)
    
  }
  
  ######################################################
  #Generate m data sets
  m = 1 # m = 50 (in paper)
  omega_elements = t(sigma_inv)[lower.tri(sigma_inv,diag=FALSE)]
  sum(omega_elements!=0)
  
  if(precision_str == 1)
  {
    data_file_name = paste0("GHS_sim_p",p,"randomn",n,"_data")
    RDATA_file_name = paste0("GHS_sim_p",p,"randomn",n,"_data.RData")
  }else if(precision_str == 2)
  {
    data_file_name = paste0("GHS_sim_p",p,"hubs",n,"_data")
    RDATA_file_name = paste0("GHS_sim_p",p,"hubs",n,"_data.RData")
  }else if(precision_str == 3)
  {
    data_file_name = paste0("GHS_sim_p",p,"cliquespos75",n,"_data")
    RDATA_file_name = paste0("GHS_sim_p",p,"cliquespos75",n,"_data.RData")
  }else{
    data_file_name = paste0("GHS_sim_p",p,"cliquesneg45",n,"_data")
    RDATA_file_name = paste0("GHS_sim_p",p,"cliquesneg45",n,"_data.RData")
  }
  
  xx_list = list(NA)
  for (i in 1:m) {
    set.seed(2050+i)
    xx_list[[i]] = rmvnorm(n=n,mean=rep(0,p),sigma=sigma)
    
    write.table(xx_list[[i]],file=paste0(data_file_name, i, ".csv"),sep=",",row.names=FALSE,col.names=FALSE)
   }
  
  save(xx_list,sigma_inv,sigma,omega_elements,n,p,file = RDATA_file_name)
}





