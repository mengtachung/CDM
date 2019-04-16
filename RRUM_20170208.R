for ( zz in 1:20){



rm(list=ls())
##############
# data #######
##############

alpha <- function (Q, N, J, K, R){
  m = matrix(R, K, K)
  diag(m) = 1	#m <- ifelse(row(m)==col(m), 1, R)
  ch = chol(m)  # Choleski decomposition
  #u = matrix(runif(N*K), ncol=K)  	# a random matrix unif(0,1)
  u = matrix(rnorm(N*K), ncol=K)  	# a random matrix normal(0,1)
  uc = u %*% ch
  #cr = uc/max(uc)
  cr = pnorm(uc) 		# theta
  # alpha <- matrix(qbinom(cr,1,0.5), nrow=N)	# cutoff 0.5
  alpha = matrix(0, N, K)	
  for (i in 1:K) alpha[,i] = ifelse(cr[,i]>i/(K+1), 1, 0)   
  
  out=list(cr=cr, alpha=alpha)
  out
  
}

m = 1000 ############################# How many students

qq <- matrix(c(1,0,0,0,0,1,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,
               0,1,0,0,0,0,1,0,0,0,1,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0,
               0,0,1,0,0,0,0,1,0,0,0,1,0,0,1,0,0,1,1,0,1,0,0,1,1,0,1,1,0,1,
               0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,1,0,1,0,1,0,1,0,1,0,1,1,0,1,1,
               0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,1,0,1,1,0,0,1,0,1,1,0,1,1,1), nrow=30)


out=alpha(Q=qq, N=m, J=nrow(qq), K=ncol(qq), R=0.1)
#cor(out$cr)
#cor(out$alpha)
a = out$alpha ### simulated alpha

n <- matrix(1, nrow=m,ncol=30)
for(i in 1:m){
  for(j in 1:30){
    s <- 0.2
    g <- 0.2
    r <- g/(1-s) # r = 0.2/0.8 = 0.25
    cs <- (1-s)
    n[i,j]<-(cs^(qq[j,1]+qq[j,2]+qq[j,3]+qq[j,4]+qq[j,5]))*(r^((1-a[i,1])*qq[j,1]+(1-a[i,2])*qq[j,2]+(1-a[i,3])*qq[j,3]+(1-a[i,4])*qq[j,4]+(1-a[i,5])*qq[j,5]))
  }
}
pro <- matrix(runif(m*nrow(qq)), nrow=m, ncol=nrow(qq))
yy = n - pro
for(i in 1:m){
  for(j in 1:nrow(qq)){
    yy[i,j] <- ifelse(yy[i,j]>0, 1, 0)  ### yy is Simulated DATA
  }
}


#######################################################################################
#######################################################################################
Q <- matrix(c(1,0,0,0,0,1,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,
              0,1,0,0,0,0,1,0,0,0,1,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0,
              0,0,1,0,0,0,0,1,0,0,0,1,0,0,1,0,0,1,1,0,1,0,0,1,1,0,1,1,0,1,
              0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,1,0,1,0,1,0,1,0,1,0,1,1,0,1,1,
              0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,1,0,1,1,0,0,1,0,1,1,0,1,1,1), nrow=30)
#######################################################################################

as.binary <- function(x){
  ans <- NULL
  while(any(x!=0)){
    ans <- cbind(x%%2,ans)
    x <- floor(x/2)
  }
  ans
}


EstQMCMC = function(Y, K, r.start = NULL, pai.start=NULL, pi.start=NULL, q, niter){
  
  N = nrow(Y)
  J = ncol(Y)
  Y = as.matrix(Y)
  
  if(is.null(r.start))
    #r = matrix(rep(.5, J*K), nrow=J)
    r = matrix(runif(J*K, 0, 1), nrow=J)
  else
    r = r.start
  
  if(is.null(pai.start))
    #pai = rep(.5, J)
    pai = runif(J, 0, 1)
  else
    pai = pai.start
  
  if(is.null(pi.start))
    pi = runif(2^K, 0, 1)
  else
    pi = pi.start
  
  alpha.out = matrix(0, N, K)
  all.a = as.binary(0:(2^K-1))
  arrayall=array(rep((1-all.a), J),c(2^K,K,J))
  natt = apply(Q, 1, sum)
  
  
  pi.out = pi
  r.out = r
  pai.out = pai
  #rat=rep(0, J)    # accept rate
  
  for(ii in 1:niter){
    
    ############################
    # UPDATE individual ALPHAS #
    ############################
    ta=r^Q
    arrayta=array(rep(ta, 2^K),c(J,K,2^K))
    arraytap=aperm(arrayta, c(3,2,1))
    tare=arraytap^arrayall #ta^re   
    arrayprod=tare[,1,]*tare[,2,]*tare[,3,]*tare[,4,]*tare[,5,] #arrayprod=apply(tare,c(1,3),prod)
    repai=t(matrix(rep(pai,2^K),nrow=J))
    predall=t(arrayprod*repai)
    
    ll=Y%*%log(predall)+(1-Y)%*%log(1-predall)
    ll = sweep(ll,2,log(pi),'+')
    pp = exp(ll)
    pp = apply(pp,1,cumsum)
    pp = sweep(pp,2,pp[2^K,],'/')
    u = runif(N)
    alpha = apply(sweep(pp,2,u,'<'),2,sum)
    alpha=as.binary(c(2^K-1,alpha))[-1,]
    alpha.out = rbind(alpha.out,alpha)
    #############
    # Update pi #
    #############
    cc = as.vector(alpha%*%(2^((K-1):0)))
    cc = apply(outer(cc, 0:(2^K-1), '=='), 2, sum)
    
    pi = rgamma(2^K, 1+cc)
    pi = pi/sum(pi)
    #pi.out = rbind(pi.out, pi)
    
    #############
    # UPDATE r  #
    #############
    arrayalpha=array(rep((1-alpha), J),c(N,K,J))
    arrayta=array(rep(ta, N),c(J,K,N))
    arraytap=aperm(arrayta, c(3,2,1))
    tare=arraytap^arrayalpha #ta^re
    
    arrayprod=tare[,1,]*tare[,2,]*tare[,3,]*tare[,4,]*tare[,5,] #arrayprod=apply(tare,c(1,3),prod)
    repai=t(matrix(rep(pai,N),nrow=J))
    preda=t(arrayprod*repai)
    
    tmp=Y*log(t(preda))+(1-Y)*log(t(1-preda)) #1000x30
    
    d0 <- colSums(tmp)+ rowSums(dbeta(r, 1, 1, log=TRUE))
    
    r1 = matrix(runif(J*K, r-0.052, r+0.052),J,K)
    r1=ifelse(r1<0, 1e-8, r1)
    r1=ifelse(r1>1, 1-1e-8, r1)
    
    ta=r1^Q   
    arrayta=array(rep(ta, N),c(J,K,N))
    arraytap=aperm(arrayta, c(3,2,1))
    tare=arraytap^arrayalpha	#ta^re
    
    arrayprod=tare[,1,]*tare[,2,]*tare[,3,]*tare[,4,]*tare[,5,]	#arrayprod=apply(tare,c(1,3),prod)
    preda=t(arrayprod*repai)
    
    tmp=Y*log(t(preda))+(1-Y)*log(t(1-preda))	#1000x30
    
    d1 <- colSums(tmp)+ rowSums(dbeta(r1, 1, 1, log=TRUE))
    
    accept <- d1-d0
    accept <-ifelse(accept>0,0, accept)
    accept <- ifelse(runif(J)<exp(accept),1,0) 
    
    accept <- which(accept==1)
    r[accept,] <- r1[accept,]
    r.out = rbind(r.out,r)
    
    ##############
    # update pai #
    ##############
    ta = r^Q
    arrayta=array(rep(ta, N),c(J,K,N))
    arraytap=aperm(arrayta, c(3,2,1))
    tare=arraytap^arrayalpha #ta^re
    
    arrayprod=tare[,1,]*tare[,2,]*tare[,3,]*tare[,4,]*tare[,5,] #arrayprod=apply(tare,c(1,3),prod)
    preda=t(arrayprod*repai)
    
    tmp=Y*log(t(preda)) + (1-Y) * log(t(1-preda)) #1000x30
    
    d0 <- colSums(tmp) + dbeta(pai, 1, 1, log=TRUE)
    
    pai1 <- runif(J, pai-0.052, pai+0.052)
    pai1 <- ifelse(pai1<0, 1e-8, pai1)
    pai1 <- ifelse(pai1>1, 1-1e-8, pai1)
    
    repai1 = t(matrix(rep(pai1,N), nrow=J))
    preda1 = t(arrayprod*repai1)
    
    tmp1 = Y*log(t(preda1))+(1-Y)*log(t(1-preda1))	#1000x30
    
    d1 <- colSums(tmp1)+ dbeta(pai1, 1, 1, log=TRUE)
    
    accept <- d1-d0
    accept <-ifelse(accept>0, 0, accept)
    accept <- ifelse(runif(J)<exp(accept), 1, 0) 
    accept <- which(accept==1)
    pai[accept] <- pai1[accept]
    pai.out = rbind(pai.out, pai)
  }
  #rat = rat/niter
  r.out = array(r.out, c(J, niter, K))
  out = list(alpha=alpha.out, pi=pi.out, pai=pai.out, r=r.out)
  class(out) = 'cdmcmc'
  out
}

#source('EstQMCMC.R')
system.time(out <- EstQMCMC(yy, K=5, q=Q, niter=7000))  # the system time reports how long it took #
################################
bn=2000

#apply(out$alpha[,-(1:bn),], c(1,3), mean) 
#apply(out$pai[-(1:bn),],2,mean)
################################
#Q*apply(out$r[,-(1:bn),], c(1,3), mean) #the -(1:200) discards the first 200 #
################################

N = 1000
K = 5
TOT = N*K
alpha.out = array(out$alpha, c(N, 7001,K))
#alpha.out[,1,]
alphaout = apply(alpha.out[,-(1:bn),], c(1,3), mean) 
#write.csv(alphaout, file="alphaout.casv")

#accr1=(5000-sum(abs(alphaout-a)))/5000	### accuracy rate
accr2=(TOT-sum(abs(round(alphaout, digits=0)-a)))/TOT	### accuracy rate
#print(accr1)
print(accr2)


}
