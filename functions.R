library(pracma)
library(igraph)


m_s=function(y,s){
  a=mean(y^s)
  a=a^(1/s)
  return(a)
}
obj_power.k.means=function(X,theta,s){
  s1=0
  n=dim(X)[1]
  k=dim(theta)[1]
  y=numeric(k)
  for(i in 1:n){
    for(j in 1:k){
      y[j]=sum((X[i,]-theta[j,])^2)
    }
    s1=s1+m_s(y,s)
  }
  return(s1)
}

grad_theta=function(X,theta,s){
  #naieve code
  epsa=1e-5
  n=dim(X)[1]
  p=dim(X)[2]
  k=dim(theta)[1]
  m=matrix(0,n,k)
  coef=numeric(n)
  gr.theta=matrix(0,k,p)
  for(i in 1:n){
    for(l in 1:k){
      m[i,l]=sum((X[i,]-theta[l,])^2)+1e-4
    }
    coef[i] = (sum(m[i,]^s))^((1/s)-1)
  }
  for(j in 1:k){
    for(l in 1:p){
      for(i in 1:n){
        gr.theta[j,l]=gr.theta[j,l]+((m[i,j])^(s-1))*coef[i]*(X[i,l]-theta[j,l])
      }
    }
  }
  return(-gr.theta)
}


MOMPKM=function(X,k,L,s=-1,eta=1.02,alpha=0.1,verbose=FALSE,tmax=100){
  n=dim(X)[1]
  B=floor(n/L)
  X=X[1:(L*B),]
  n=dim(X)[1]
  p=dim(X)[2]
  permutation=sample(n,n)
  mu=median(X)
  f_val=numeric(tmax)
  val=numeric(L)
  ind=numeric(tmax)
  theta=X[sample(n,k),]+rand(k,p)*0.01
  b=rep(0,k)
  for(t in 1:tmax){
    for(l in 1:L){
      I=permutation[((l-1)*B+1):(l*B)]
      val[l]=obj_power.k.means(X[I,],theta,s)
    }
    val_prime=sort(val)
    m=val_prime[floor((L+1)/2)]
    #m=median(val)
    f_val[t]=m
    if((t>5)&&(abs((f_val[t]-f_val[t-1])/f_val[t-1])<1e-4)){
      break
    }
    index=which(val==m)
    ind[t]=index
    I=permutation[((index-1)*B+1):(index*B)]
    gr.theta=grad_theta(X[I,],theta,s)
    b=b+rowSums(gr.theta^2)
    theta=theta-alpha/sqrt(b)*gr.theta
    if(t%%2==0){
      s=s*eta
    }
   # permutation=sample(n,n)
    if(verbose==TRUE){
      if(t %%50 ==0){
        cat('Iteration no. ')
        cat(t)
        cat('\n')
      }
    }

  }
  label=numeric(n)
  dist=numeric(k)
  for(i in 1:n){
    for(j in 1:k){
      dist[j]=sum((X[i,]-theta[j,])^2)
    }
    label[i]=which.min(dist)
  }
  list1=list(label,theta)
  names(list1)=c('label','theta')
  #list(label,theta,s,f_val,b)
  return(list1)
}

data_generate=function(n,M,prob1,sigma,sigma2){
  p=dim(M)[2]
  X=matrix(0,n,p)
  k=dim(M)[1]
  label=numeric(n)
  for(i in 1:n){
    s=sample(1:k,size=1,prob=prob1)
    for(l in 1:p){
      if(M[s,l]==0){
        X[i,l]=rnorm(1,M[s,l],sigma2)
      }else{
        X[i,l]=rnorm(1,M[s,l],sigma)
      }

    }
    label[i]=s
  }
  ls=list(X,label)
  names(ls)=c('data','label')
  return(ls)
}

power.k.means=function(X,s=-1,k,eta=1.04,
                       tmax=200,tol=10^(-3)){
  n=dim(X)[1]
  p=dim(X)[2]
  W=matrix(0,nrow=n,ncol=k)
  m=matrix(0,nrow=n,ncol=k)
  label=numeric(n)
  dd=numeric(k)
  #initialization
  sam=sample(n,k)
  theta=X[sam,]
  #theta=rand(k,p)
  #theta=rand(k,p)
  for(i in 1:p){
    theta[,i]=theta[,i]+rnorm(k,0,0.1)
  }
  val1=theta
  coef = numeric(n)
  for(t in 1:tmax){
    # if(t==tmax){
    #   cat(coef)
    #   cat("\n")
    #   cat(W[1,])
    #   cat("\n")
    # }
    #update weights
    for(i in 1:n){
      for(l in 1:k){
        m[i,l]=sum((X[i,]-theta[l,])^2)
      }
      coef[i] = (sum(m[i,]^s))^((1/s)-1)
    }
    for(i in 1:n){
      for(l in 1:k){
        W[i,l]=m[i,l]^(s-1)*coef[i]
      }
    }
    #update centroids
    for(l in 1:k){
      for(d in 1:p){
        theta[l,d]=(sum(W[,l]*X[,d]))/sum(W[,l])
      }
    }
    #update s
    if(t%%2==0){s=eta*s}
    # cat(t)
    # cat('\n')
    t=t+1
    val2=theta
    if(norm(val1-val2)<tol){
      break
    }else{
      val1=val2
    }
  }
  for(i in 1:n){
    for(j in 1:k){
      dd[j]=sum((X[i,]-theta[j,])^2)
    }
    label[i]=which.min(dd)
  }
  list1=list(theta,label)
  names(list1)=c('theta','label')
  return(list1)
}