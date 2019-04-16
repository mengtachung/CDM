for(lll in 1:100){
  
  ##########################################
  ##########################################
  rm(list=ls())
  
  DINA.SIM <- function (Q, N, J, K, R){
    m = matrix(R, K, K)
    diag(m) = 1	#m <- ifelse(row(m)==col(m), 1, R)
    ch = chol(m)  # Choleski decomposition
    #u = matrix(runif(N*K), ncol=K)  	# a random matrix unif(0,1)
    u = matrix(rnorm(N*K), ncol=K)  	# a random matrix normal(0,1)
    uc = u %*% ch
    #cr = uc/max(uc)
    cr = pnorm(uc) 		# theta
    # alpha <- matrix(qbinom(cr,1,0.5), nrow=N)	# cutoff 0.5
    alpha <- matrix(ifelse(runif(N*K)>cr, 1, 0), nrow=N)	# inverse
    #alpha = matrix(0, N, K)	
    #for (i in 1:K) alpha[,i] = ifelse(cr[,i]>i/(K+1), 1, 0)   
    
    tm=alpha%*%t(Q) 
    at=matrix(rep(apply(Q, 1, sum), N), nrow=N, byrow=T)
    eta=ifelse(tm==at, 1, 0)
    y=ifelse(eta==1, 0.8, 0.2) ## s=0.2; g=0.2 ##
    comp=c(runif(N*J, 0, 1))
    y=ifelse(y>comp, 1, 0)
    
    out=list(cr=cr, alpha=alpha, y=y)
    out
  }
  
  ####### de la Torre ##################################################################
  #q <- matrix(c(1,0,0,0,0,1,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,
  #              0,1,0,0,0,0,1,0,0,0,1,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0,
  #              0,0,1,0,0,0,0,1,0,0,0,1,0,0,1,0,0,1,1,0,1,0,0,1,1,0,1,1,0,1,
  #              0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,1,0,1,0,1,0,1,0,1,0,1,1,0,1,1,
  #              0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,1,0,1,1,0,0,1,0,1,1,0,1,1,1), nrow=30)
  ########################################################################################

  ############ subtraction fraction ###################################
  q = matrix(c(0,1,0,1,0,
               0,1,0,0,1,
               0,0,1,1,0,
               0,0,1,0,1,
               0,0,0,1,1,
               1,1,1,0,0,
               1,1,0,1,0,
               1,1,0,0,1,
               1,0,1,1,0,
               1,0,1,0,1,
               1,0,0,1,1,
               0,1,1,1,1,
               0,1,1,0,1,
               0,0,0,1,1,
               0,0,1,1,1), byrow=TRUE, nrow=15, ncol=5)
  ##################################################################################
  out=DINA.SIM(Q=q, N=500, J=nrow(q), K=ncol(q), R=0.1)
  cor(out$cr)
  cor(out$alpha)
  yy = out$y
  ################################################################################
  ################################################################################
  ################################################################################
  ################################################################################
  
  as.binary <- function(x){
    ans <- NULL
    while(any(x!=0)){
      ans <- cbind(x%%2,ans)
      x <- floor(x/2)
    }
    ans
  }
  
  
  EstQMCMC = function(Y, K=NULL,q.start = NULL, g.start=NULL, s.start=NULL,pi.start=NULL, a.start, niter){
    
    N = nrow(Y)
    J = ncol(Y)
    Y = as.matrix(Y)
    
    if(is.null(K) & is.null(q.start))
      stop('User must supply either the number of attributes or a starting Q matrix!\n\n')
    
    if(is.null(q.start)){
      Q = matrix(rbinom(J*K, 1, 0.5), J, K)
      Q[which(apply(Q, 1, sum)==0), 1] = 1  # assign (1,0,0,0,0)
    }
    else{
      Q = q.start
    }
    
    if(is.null(g.start))
      g = runif(J, 0.1, 0.5)    # g = rep(.2,J)
    else
      g = g.start
    
    if(is.null(s.start))
      s = runif(J, 0.1, 0.5)    # s = rep(.2,J)
    else
      s = s.start
    
    if(is.null(pi.start)){
      pi = exp(rnorm(2^K))
      pi = pi/sum(pi)
    }
    else
      pi = pi.start
    
    all.a = as.binary(0:(2^K-1))
    natt = apply(Q,1,sum)
    
    Yt = t(Y)
    a = a.start
    p = matrix(runif(J*K), J, K)
    p.ls = p
    pi.out = pi; Q.out=Q; p.out=p; g.out=g; s.out=s
    
    for(ii in 1:niter){
      ############################
      # UPDATE individual ALPHAS #
      ############################
      etajm = tcrossprod(Q,all.a)
      natt = apply(Q,1,sum)
      etajm = (etajm == natt)      # determines which attribute patterns have required skills #
      pp = g*(1-etajm) + (1-s)*etajm
      ll = Y %*% log(pp) + (1-Y)%*%log(1-pp)
      ll = sweep(ll,2,log(pi),'+')
      pp = exp(ll)
      pp = apply(pp,1,cumsum)
      pp = sweep(pp,2,pp[2^K,],'/')
      u = runif(N)
      alpha = apply(sweep(pp,2,u,'<'),2,sum)
      alpha=as.binary(c(2^K-1,alpha))[-1,]
      
      ##############
      # Update pi  #
      ##############
      cc = as.vector(alpha%*%(2^((K-1):0)))
      cc = apply(outer(cc,0:(2^K-1),'=='),2,sum)
      
      pi = rgamma(2^K, 1+cc)
      pi = pi/sum(pi)
      pi.out = rbind(pi.out,pi)
      
      #########################
      # UPDATE Guess and Slip #
      #########################
      etaim = tcrossprod(Q,alpha)
      etaim = (etaim == natt)  # determines which individuals have required skills #
      ga = apply((1-etaim)*Yt,1,sum)
      gb = apply((1-etaim)*(1-Yt),1,sum)
      sa = apply(etaim*(1-Yt),1,sum)
      sb = apply(etaim*Yt,1,sum)
      
      g = qbeta(runif(J, 0,pbeta(1-s,1+ga,2+gb)),1+ga,2+gb)
      s = qbeta(runif(J,0,pbeta(1-g,1+sa,2+sb)),1+sa,2+sb)
      g.out = rbind(g.out,g)
      s.out = rbind(s.out,s)
      
      ####################
      # UPDATE Q matrix  #
      ####################
      Q = t(sapply(1:J,function(j){sample.Q(all.a[-1,],g[j],s[j],Y[,j],alpha,p[j,])}))
      Q.out = rbind(Q.out,Q)
      
      #############
      # Update P  #
      #############
      pa = a + Q; pb = a +1-Q
      ppp = rbeta(J*K, pa, pb)
      p = matrix(ppp, J, K)
      p.out = cbind(p.out, p)
      
      #############
      # relable p #
      #############
      #relout = reorder(J=15, K=5, a=p.ls, b=p)
      #p = relout$reorder.b
      
    }
    p.out = array(p.out, c(J,K,niter))
    Q.out = array(Q.out, c(J,niter,K))
    out = list(pi=pi.out, Q=Q.out, g=g.out, s=s.out, p=p.out)
    class(out) = 'cdmcmc'
    out
  }
  
  sample.Q = function(all.a, g, s, Y, alpha,pp){
    natt = apply(all.a,1,sum)
    cc = tcrossprod(all.a,alpha)
    etaim = (cc==natt)
    pp[pp<1e-8] = 1e-8
    pp[pp>1-1e-8] = 1-1e-8
    pp = all.a%*%log(pp) + (1-all.a)%*%log(1-pp)
    ga = (1-etaim)%*%Y
    gb = (1-etaim)%*%(1-Y)
    sa = etaim%*%(1-Y)
    sb = etaim%*%Y
    pm = ga*log(g) + gb*log(1-g) + sa*log(s) + sb*log(1-s)
    pm = pm + pp
    pm = pm -max(pm)
    pm = as.vector(exp(pm))
    pm = pm/sum(pm)
    kk = nrow(all.a)
    q = sample(1:kk,size=1,prob=pm)
    q = (as.binary(c(q,kk))[1,])
    q
  }
  ######################################################################
  
  #source('dina_0325_label.r')
  
  system.time(out <- EstQMCMC(yy, K=5, a.start=1, niter=22000))
  
  bn=2000
  
  apply(out$Q[,-(1:bn),], c(1,3), mean)
  qest =  apply(out$Q[,-(1:bn),], c(1,3), mean)
  
  apply(out$p[,,-(1:bn)],c(1,2),mean)
  
  apply(out$g[-(1:bn),],2,mean)
  
  apply(out$s[-(1:bn),],2,mean)
  
  apply(out$pi[-(1:bn),],2,mean)
  
  ### estimated Q-matrix
  #qest=round(apply(out$Q[,-(1:bn),], c(1,3), mean), digits=0) 
  
  
  ######################## Original Q-matrix ###########################################
  #q <- matrix(c(1,0,0,0,0,1,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,
  #              0,1,0,0,0,0,1,0,0,0,1,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0,
  #              0,0,1,0,0,0,0,1,0,0,0,1,0,0,1,0,0,1,1,0,1,0,0,1,1,0,1,1,0,1,
  #              0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,1,0,1,0,1,0,1,0,1,0,1,1,0,1,1,
  #              0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,1,0,1,1,0,0,1,0,1,1,0,1,1,1), nrow=30)
  #
  ######################################################################################
  q = matrix(c(0,1,0,1,0,
0,1,0,0,1,
0,0,1,1,0,
0,0,1,0,1,
0,0,0,1,1,
1,1,1,0,0,
1,1,0,1,0,
1,1,0,0,1,
1,0,1,1,0,
1,0,1,0,1,
1,0,0,1,1,
0,1,1,1,1,
0,1,1,0,1,
0,0,0,1,1,
0,0,1,1,1
), byrow=TRUE, nrow=15, ncol=5)
  ########################################################################################
  qqq = 0.6*q
  ############## Relable so that qest is closest to q ####################################
  permutations <- function(n){
    if(n==1){
      return(matrix(1))
    } else {
      sp <- permutations(n-1)
      p <- nrow(sp)
      A <- matrix(nrow=n*p,ncol=n)
      for(i in 1:n){
        A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
      }
      return(A)
    }
  }
  
  reorder <- function(J, K, a, b){
    vec.a = matrix(a, ncol=1) #as.vector(a)  # starting p matrix : p.ls
    vec.b = matrix(b, ncol=1) #as.vector(b)  # matrix of p
    pm = permutations(K) ## 5x4x3x2x1 = 120
    tpm =t(pm)
    vec.b.matrix = matrix(as.vector(b[,c(tpm[,1:factorial(K)])]), J*K, factorial(K))
    vec.bind.ab = as.matrix(cbind(vec.a, vec.b.matrix))
    dist.matrix = as.matrix(dist(t(vec.bind.ab), method="euclidean"))
    ds = dist.matrix[,1]
    min.value = (min(ds[ds>0]))
    matrix.number = as.numeric(which(ds == min.value))
    reorder.b = matrix(vec.b.matrix[, matrix.number-1], J, K)    ## result of reordered matrix  
    output = list(reorder.b=reorder.b)
    output$reorder.b ####################
  }
  ######################################################################
  #### original Q-matrix estimate; no relabeling yet ###################
  ######################################################################
  rqest=reorder(J=15, K=5, a=qqq, b=qest)   ### reorder qest
  
  diff1 = rqest-q
  delta1 = 1-(sum(abs(diff1))/(15*5))
  #print(delta1)   ###print(1-(sum(abs(diff1))/(15*5)))
  cat("delta1 (original) =", delta1, "\n")
  
  ######################################################################
  #### relabeling starts here ##########################################
  ######################################################################
  QA = array(dim=c(15,5,20000)) ## FOR later store
  QB = out$Q[,-(1:bn),]  ## Q-matrix estimates
  
  for (rr in 1:20000){
    QA[,,rr]=reorder(15, 5, qest, QB[,rr,])
  }
  Qest=apply(QA[,,(1:20000)], c(1,2), mean) ### USE THIS FOR LATER REFERENCE
  #Qest=round(Qest, digits=0)
  
  Rqest = reorder(J=15, K=5, a=qqq, b=Qest)
  diff2 = Rqest - q
  delta2 = 1-(sum(abs(diff2))/(15*5))
  #print(delta2)   ### print(1-(sum(abs(diff2))/(15*5)))  ### result after reorder
  cat("delta2=", delta2, "\n")
    
  ##################################################################################
  for(mmm in 1:100){    ### assuming 100 iterations are needed   
    for (rr in 1:20000){
      QA[,,rr]=reorder(15, 5, Qest, QA[,,rr])
    }
    RQest=apply(QA[,,(1:20000)], c(1,2), mean)  

########################## Snooping #############################################
    TRQest = reorder(J=15, K=5, a=qqq, b=RQest)  ### snooping
    TRQest = 1-(sum(abs(TRQest-q)))/(15*5)
    cat("delt = ", TRQest, "\n") ### snooping, comment out this line later
#################################################################################

    dif = (sum(abs(RQest-Qest)))/(15*5)    
    cat("dif=", dif, "\n")  
    if (dif < 0.0001){    
      RQest = reorder(J=15, K=5, a=qqq, b=RQest)
      diff3 = RQest - q
      delta3 = 1-(sum(abs(diff3))/(15*5))
      #print(delta3)
      cat("delta3 (final) =", delta3, "\n")
      break
    }
    
    Qest = RQest
  }
  ######################################################################################### ##########################################################################################
  
  ###################################################
  ### break to see result if diff2>diff1 ############
  ### Check Rqest (relabled) and rqest (original) ###
  ###################################################

  #if (delta3-delta1>0.15){

    ####################################################################################
    ######### If break, then check Kernel Density for [1,1] of the 150 entries #########
    ####################################################################################
  #  dst1=density(out$Q[1,-(1:bn),1], adjust=6)  ##bad
  #  plot(dst1)
  #  dst2=density(QA[1,1,(1:20000)], adjust=6) ## good
  #  plot(dst2)
  #  break
  #}
  
} ####### end for

