
source('~/Dropbox/R Surma/getode_parameters.R')

#getode calls the dataframe, A vector, B vector for every iteration
#myfunc calls the li,l2,l3,l4,l4n,ln1 functions to calculate the differential equation formula in ode
#getode returns myfunc : myfunc is in exact form in which ode takes function input func(t,y,parms)
#getode can input n and betamat and produce a myfunc which can go to ode

getode <- function(n,t,y,p,betamat){
 
  v <- 2^n # number of state variables

  names(y) = ynames(n) #update names of y vector
  stdata <- st.data(y,ny) #update dataframe
  avec = a.vec(y) #update a vector
  bvec = b.vec(y) #update b vector
  
  myfunc <- function(t,y,p){
    with(as.list(c(y,p)),{
      
      v <- 2^n 
      
      names(y) = ynames(n) #update names of y vector
      stdata <- st.data(y,ny) #update dataframe
      avec = a.vec(y) #update a vector
      bvec = b.vec(y) #update b vector
    
      dyode <- c() #initiate dyode vector to list dy for all state variables
      
      #if n=1, have a simple loop 
      #collect dyode for level 1 (l1), level 2(l2), da/dt (yan), db/dt (ybn)
      if(n==1){
        #if(t>1.6){
       # print("stop here")
       #  }
        dyode = c(dyode,l1(n,y,avec,p,dframe=stdata)) 
        dyode = c(dyode,l2(n,y,p,betamat,dframe=stdata))
        dyode = c(dyode,ybn(n,y,bvec,avec,p,betamat,dframe=stdata))
        dyode = c(dyode,yan(n,y,bvec,avec,p,betamat,dframe=stdata))
      }
      else if(n!=1){
        lv <- sapply(0:(n-1), function(x) choose(n,x)) #lv gives me number of levels and number of state variables in each level
        
        for(j in 1:length(lv)){
        if(j==1){
          dyode = c(dyode,l1(n,y,avec,p,dframe=stdata)) #dyode updtaed for level 1 
         }
        if(j==2){
          dyode = c(dyode,l2(n,y,p,betamat,dframe = stdata)) #dyode updtaed for level 2
        }
        if(j==3){
           dyode = c(dyode,l3(n,y,p,betamat,dframe=stdata)) #dyode updtaed for level 3
        }
        if(j>3){
          dyode = c(dyode,l4n(n,y,p,betamat,dframe=stdata,l=j)) #dyode updtaed for levels>3
        }
       }
        #(n+1) terms in the oder : b,a
      dyode = c(dyode,ln1(n,y,p,betamat,dframe=stdata)) #the last level
      #A and B terms in the oder : b,a
      #if(t>2.4){
      #print("stop here")
     # }
      dyode = c(dyode,ybn(n,y,bvec,avec,p,betamat,dframe=stdata)) #collect b terms
      dyode = c(dyode,yan(n,y,bvec,avec,p,betamat,dframe=stdata)) #collect a terms
    }
    
     return(list(dyode))
      
    })
  }
  
 myfunc
 
}

#getode(n,t,y,p,betamat)

#func1 <- getode(n,t,y,p,betamat)

