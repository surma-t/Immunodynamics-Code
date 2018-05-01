n=5
v <- 2^n

source('~/Dropbox/R Surma/getode_functions.R')

print("input y in this order")
ny=ynames(n)
print(ynames(n))

p <- c(k=0.01,df=0.5,db=0.5,s=1,phi=100,a=0.1,da=0.1,alpha=0,delta=0)  

betamat <- matrix(0,n,n)
betamat[row(betamat)!=col(betamat)]=0.0008
betamat[row(betamat)==3]=0.1
betamat[col(betamat)==3]=0.2
betamat[row(betamat)==col(betamat)]=1
betamat



#basic model
#p <- c(k=0.01,df=0.5,db=0.5,s=1,phi=100,a=0.1,da=0.1,alpha=0,delta=1)


