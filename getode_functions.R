
#function to get names for y state variables : example
#[1] "y0"   "y1"   "y2"   "y3"   "y12"  "y13"  "y23"  "y123" "b1"   "b1"   "b1"   "a1"   "a1"   "a1"  
ynames <- function(n){
  v <- 2^n # number of variables
  z <- lapply(0:(n), function(x) combn(n,x)) # matrix of states for each level
  z[[1]] <- 1 
  # the vector of states : statevec
  index=1
  vec1 <- c()
  for(i in 2:length(z)){
    for(j in 1:ncol(z[[i]])){
      vec1[index] <- as.numeric(paste(z[[i]][,j],collapse=""))
      index=index+1
    }
  }
  statevec <- c(0,vec1)
  # the names of states in the y vector 
  namevec <- c(paste0("y",statevec,collpase=""),rep("b1",n),rep("a1",n))
  return(namevec)
}


#data frame for current state level and concentration of y variables
st.data <- function(y,ny){
  lv <- sapply(0:n, function(x) choose(n,x)) 
  gp <- rep(1:(n+1),c(lv)) 
  stname <- as.numeric(substring(ny, 2)) 
  stdata <- data.frame(y=y[1:v],level=gp,state=stname[1:v]) 
  return(stdata)
}

#vector of A concentration, to ve updated in every iteration
a.vec <- function(y){
  #avec <- y[names(y)=="a1"] #A values
  avec <- y[(v+n+1):(v+2*n)]
  }
#vector of B concentration, to ve updated in every iteration
b.vec <- function(y){
   # bvec <- y[names(y)=="b1"] #B values
    bvec <- y[(v+1):(v+n)] 
}
#function to calculate concentration of free antigen: level 1 
l1 <- function(n,y,avec,p,dframe){
  with(as.list(c(p)),{
    dy1 <- 0
    stdata = dframe
    #avec <- y[names(y)=="a1"]
    #free antigen equation
    dy1 = dy1 -k*stdata[which(stdata$state==0),"y"]*sum(avec) - df*stdata[which(stdata$state==0),"y"]
    return(dy1)
  })
}


#function to calculate concentration of antigens at level 2 
l2 <- function(n,y,p,betamat,dframe){
  with(as.list(c(p)),{
    l=2
    avec <- y[names(y)=="a1"]
    stdata <- dframe
    dy2 <- c(rep(0,n))
    for(i in 1:n){
      for(nvec in 1:n){
        if(nvec!=i){ #for i=1, nvec = 2,3, dy2 = dy2-k*a2*(1-beta12)*y[1]-k*a3*(1-beta13)*y[1]
          dy2[i] = dy2[i] - k*avec[nvec]*(1-betamat[i,nvec])*stdata[which(stdata$state==i),"y"]
        }
      } #dy2 = dy2 + k*a1*y[0] - db*y[1] for i=1
      dy2[i] = dy2[i] + k*avec[i]*stdata[which(stdata$level==(l-1)),"y"] - db*stdata[which(stdata$state==i),"y"]
    }
    dy2
  })
}

#function to calculate concentration of antigens at level 3
l3 <- function(n,y,p,betamat,dframe){
  with(as.list(c(p)),{
    l=3
    nvstate=0
    avec <- y[names(y)=="a1"]
    stdata <- dframe
    dy3 <- c(rep(0,choose(n,(l-1)))) 
    lmat=combn(n,(l-1))
    #loop for all possible state variables in this level
    for(i in 1:ncol(lmat)){
      nmat = lmat[,i] # current state in number : c(1,2) or c(3,4)
      nvstate = as.numeric(paste0(nmat,collapse="")) # current state in paste form "12" to match with data frame
      nvec = 1:n
      newvec <- nvec[!(nvec%in%nmat)] #newvec has epitope numbers to be added to the current state : 12 -> 123 : newvec=3
      bterm <- rep(1,(length(nmat))) #bterm gets the (1-beta12) terms
      
      #get the in terms: for i=1, state=12, dy3=k*avec[2](1-beta12)y[1]+k*avec[1](1-beta21)y[2]
      for(j in 1:length(nmat)){ 
        nvec = 1:n
        n2 = nmat[j]
        n1 = nvec[(nvec%in%nmat) & (nvec!=nmat[j])]
        dy3[i] = dy3[i] + k*avec[n1]*(1-betamat[n2,n1])*stdata[which(stdata$state==n2),"y"]
      }
      
      #get the out terms : for i=1,state=12->123 dy3 =-k*a3*(1-beta13)*(1-beta23)*y[12]-k*a4*(1-beta14)*(1-beta24)*y[12] 
      for(nv in 1:length(newvec)){  # newvec=3,4
        bterm <- rep(1,length(nmat)) # initiate bterm
        for(ih in 1:length(nmat)){  #nmat=1,2
          bterm[ih] = bterm[ih]*(1-betamat[nmat[ih],newvec[nv]]) # (1-beta13) (1-beta23)
        }
        dy3[i] = dy3[i] - k*avec[newvec[nv]]*prod(bterm)*stdata[which(stdata$state==nvstate),"y"]
      }
      
      #get the dy3 = df*y[12] term
      dy3[i] = dy3[i] - df*stdata[which(stdata$state==nvstate),"y"]
    }
    dy3
  })
}

#for l>=4 for any n : cross checked with hand-written calculations
l4n <- function(n,y,p,betamat,dframe,l){
  with(as.list(c(p)),{
    avec <- y[names(y)=="a1"]
    stdata <- dframe
    nvstate=0; ystate=0; 
    lmat=combn(n,(l-1))
    dyn <- c(rep(0,ncol(lmat))) 
    
    #loop through all possible states in a level
    for(i in 1:ncol(lmat)){
      nmat = lmat[,i] # i=1, nmat = 1 2 3
      ystate=as.numeric(paste0(nmat,collapse="")) #nvstate = 123
      nvec = 1:n
      
      #get the in terms: for i=1,l=4,  state=123, 
      #dyn=k*avec[1](1-beta21)(1-beta31)y[23]+k*avec[2](1-beta12)(1-beta32)y[13]+k*avec[3](1-beta13)(1-beta23)y[12]
      
      for(j in 1:length(nmat)){ #j = 1,2,3
        nv = nvec[(nvec%in%nmat) & (nvec != nmat[j])] #nv=23 for j=1
        nvstate = as.numeric(paste0(nv,collapse="")) 
        bterm = rep(1,length(nv)) # bterm to collect all the (1-beta) terms
        
        for(ih in 1:length(nv)){
          bterm[ih] = bterm[ih]*(1-betamat[nv[ih],nmat[j]]) #(1-beta21) (1-beta31) terms
        }
        dyn[i]=dyn[i] + k*avec[nmat[j]]*prod(bterm)*stdata[which(stdata$state==nvstate),"y"]
      }
      #get the out terms : for i=1,l=4, state=123->1234 
      #dyn =-k*avec[4](1-beta14)(1-beta24)(1-beta34)y[123]-k*avec[5](1-beta15)(1-beta25)(1-beta35)y[123] 
      
      newvec = nvec[!(nvec%in%nmat)] # i=1, nmat=123, newvec=45, n=5
      bterm2 = rep(1,length(nmat))
      for(nw in 1:length(newvec)){
        for(n2 in 1:length(nmat)){
          bterm2[n2] = bterm2[n2]*(1-betamat[nmat[n2],newvec[nw]]) # collect (1-beta) terms
        }
        dyn[i] = dyn[i] - k*avec[newvec[nw]]*prod(bterm2)*stdata[which(stdata$state==ystate),"y"] #add the outgoing terms to dy
      }
      dyn[i] = dyn[i] - db*stdata[which(stdata$state==ystate),"y"] #add the db*y[state] term
    }
    dyn
  })
}


#for l=n+1 for any n: cross checked with hand-written calculations
ln1 <- function(n,y,p,betamat,dframe){
  with(as.list(c(p)),{
    l=n+1
    avec <- y[names(y)=="a1"]
    stdata <- dframe
    nvstate=0; ystate=0; 
    lmat <- combn(n,(l-1))
    dyn1 <- c(rep(0,ncol(lmat))) 
    nmat = lmat[,1] # n=4, nmat = 1 2 3 4
    ystate=as.numeric(paste0(nmat,collapse="")) #nvstate = 1234
    nvec = 1:n # nvec = 1,2,3,4
    
    #for n=4, l=5,  state=1234, 
    #dyn=k*avec[1](1-beta21)(1-beta31)(1-beta41)y[234]+k*avec[2](1-beta12)(1-beta32)(1-beta42)y[134]
    #   +k*avec[3](1-beta13)(1-beta23)(1-beta43)y[124]+k*avec[4](1-beta14)(1-beta24)(1-beta34)y[123]
    #loop through all possible states in a level
    for(j in 1:length(nmat)){ #j = 1,2,3,4
      nv = nvec[(nvec%in%nmat) & (nvec != nmat[j])] #nv=234 for j=1
      nvstate = as.numeric(paste0(nv,collapse="")) 
      bterm = rep(1,length(nv)) # bterm to collect all the (1-beta) terms
      
      for(ih in 1:length(nv)){
        bterm[ih] = bterm[ih]*(1-betamat[nv[ih],nmat[j]]) #(1-beta21) (1-beta31) (1-beta41) terms
      }
      dyn1[1]=dyn1[1] + k*avec[nmat[j]]*prod(bterm)*stdata[which(stdata$state==nvstate),"y"] # collect only the incoming terms
    }
    dyn1[1] = dyn1[1] - db*stdata[which(stdata$state==ystate),"y"] #add the db*y[state] term
    dyn1
  })
}


#for da1/dt terms: cross checked with hand-written calculations
yan <- function(n,y,bvec,avec,p,betamat,dframe){
  with(as.list(c(p)),{
    stdata <- dframe
    #avec <- y[names(y)=="a1"]
    #bvec <- y[names(y)=="b1"]
    nvstate=0
    daterm <- c(rep(0,n))

    # kAy[0] term
    for(i in 1:n){
      for(lvl in 1:n){
       # print(c(1,i,lvl))
        if(lvl == 1){
          daterm[i] <- daterm[i] - k*avec[i]*y[1]
         # print(c(2,i,lvl,daterm))
        }
        else if(lvl==2){
          
          # loop for k*A*(1-beta)*y terms
          for(nvec in 1:n){
            #print(c(3,i,lvl,nvec,daterm))
            if(nvec != i){
              daterm[i] = daterm[i] - k*avec[i]*stdata[which(stdata$state==nvec),"y"]*(1-betamat[nvec,i])
              #print(c(4,i,lvl,nvec,daterm))
            }
          }
        }
        # loop for terms: k*A*y*(1-beta)(1-beta)... 
        else if(lvl>2){
          cmat <- combn(n,(lvl-1))
          for(j in 1:ncol(cmat)){
            nmat <- cmat[,j]
            if(!(i %in% nmat)){
              nvstate=as.numeric(paste0(nmat,collapse=""))
              bterm <- rep(1,length(nmat))
              for(ih in 1:length(nmat)){
                bterm[ih] = bterm[ih]*(1-betamat[nmat[ih],i]) # collect the beta term
              }
              daterm[i] = daterm[i] - k*avec[i]*stdata[which(stdata$state==nvstate),"y"]*prod(bterm)
            }
          }
        }
      }
     # print(c(4,i,lvl,nvec,daterm))
      daterm[i] =  daterm[i] + a*bvec[i] - da*avec[i] 
     # print(c(5,i,lvl,nvec,daterm))
    }
    daterm
  })
}



#for db1/dt terms: cross checked with hand-written calculations
ybn <- function(n,y,bvec,avec,p,betamat,dframe){
  with(as.list(c(p)),{
    stdata <- dframe
    #avec <- y[names(y)=="a1"]
    #bvec <- y[names(y)=="b1"]
    nvstate=0
    dbnume <- c(rep(0,n))
    dbterm <- c(rep(0,n))
    for(i in 1:n){
      for(lvl in 1:(n+1)){
        #print(c(1,i,lvl,dbnume))
        if(lvl == 1){
          #dbnume[i] <- dbnume[i] + y[1]
        }
        else if(lvl==2){
          for(nvec in 1:n){
           # print(c(2,i,lvl,nvec,dbnume))
            if(nvec == i){
              dbnume[i] = dbnume[i] + delta*stdata[which(stdata$state==nvec),"y"]
            }
            else if(nvec != i){
             # print(c(3,i,lvl,nvec,dbnume))
              dbnume[i] = dbnume[i] + stdata[which(stdata$state==nvec),"y"]*(1-betamat[nvec,i])
            }
          }
        }
        else if(lvl>2){
          cmat <- combn(n,(lvl-1))
          for(j in 1:ncol(cmat)){
            nmat <- cmat[,j]
            if(!(i %in% nmat)){
              nvstate=as.numeric(paste0(nmat,collapse=""))
              bterm <- rep(1,length(nmat))
              for(ih in 1:length(nmat)){
                bterm[ih] = bterm[ih]*(1-betamat[nmat[ih],i])
              }
              dbnume[i] = dbnume[i] + stdata[which(stdata$state==nvstate),"y"]*prod(bterm)
            }
            else if(i %in% nmat){
              nvstate=as.numeric(paste0(nmat,collapse=""))
              dbnume[i] = dbnume[i] + delta*stdata[which(stdata$state==nvstate),"y"]
            }
          }
        }
      }
      dbnume[i] = dbnume[i] + y[1]
      #print(c(4,i,dbnume,dbterm))
      dbterm[i] = dbterm[i]+(s*bvec[i]*dbnume[i])/(phi+dbnume[i])
      #print(c(5,i,dbnume,dbterm))
    }
    dbterm
  })
}

#function to create root function for any n
rHf <- function(n,t,y,parms){
  v=2^n
  rootvec <- c()
  
  myfunc <- function(t,y,parms){
    for(i in 1:v){
      rootvec <- c(rootvec,(y[i]+abs(y[i])-exp(log(1e-300))))    
    }
    return(rootvec)
  }
  myfunc
}

#function to create event function for any n 
eHf <- function(n,t,y,parms){
  v=2^n
  eventvec <- c()
  
  myfunc <- function(t, y, parms) {
    
    for(i in 1:v){
      y[i][y[i] + abs(y[1]) <= 0 ] <- 0
    }
    
    for(j in 1:(v+2*n)){
      eventvec <- c(eventvec,y[j])
    }
    return (eventvec)
  }
  myfunc
}

