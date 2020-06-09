install.packages('devtools')
install.packages('roxygen2')
install.packages('usethis')
library(devtools)
library(roxygen2)

setwd('/Users/mostafa/Mostafa/UNC/Research/Proposal/OurRpackage')
create('ILFMM')

ILFMM<- function(Y, subj, t, func1, func2, covariates=NULL, comm.pen=TRUE,  
                 pentype='Ridge', L.user1=NULL,L.user2=NULL, f1_t=NULL,f2_t=NULL, Q1=NULL,Q2=Null, a=10^3,c=10^3,
                 se=FALSE, ...)
{
  require(nlme)
  require(magic)
  lobj=TRUE
  
  W<- as.matrix(func1)
  D<-as.matrix(func2)
  #norm<- sqrt(diag(W%*%t(W)))
  #W<- W/norm
  K<- ncol(W) 
  
  Y<- as.matrix(Y)
  id<- as.matrix(subj)
  t<-as.matrix(t)
  
  #Check 1:Making sure Y, subj and t have only 1 column
  if(dim(Y)[2]>1) return(cat("Error: No. of column for Y cannot be greater than 1. \nThe ILFMM() will not proceed further.\n"))
  if(dim(id)[2]>1) return(cat("Error: No. of column for subj cannot be greater than 1. \nThe ILFMM() will not proceed further.\n"))
  if(dim(t)[2]>1) return(cat("Error: No. of column for t cannot be greater than 1. \nThe ILFMM() will not proceed further.\n"))
  
  Yl<- dim(Y)[1]
  idl<- dim(id)[1]
  tl<- dim(t)[1]
  
  #Check 2: Do check for intercept in X matrix
  if (!is.null(covariates)){
    covariates<- as.matrix(covariates)
    X.chk<- apply(covariates, 2, sd)
    if(any(X.chk==0)) return(cat("Error: Drop intercept or equivalent to intercept term from covariate. \nThe ILFMM() will not proceed further.\n"))
  }
  X<-  cbind(1, covariates)
  
  #Check 3: Check the dimension of Y, id, t, W, D and X
  chk.eq<- ifelse(Yl==idl & idl==tl & tl==nrow(W) &tl==nrow(D), 0 ,1)
  if(chk.eq==1) return(cat("Error: At least one of (1) length of Y, (2) lenght of subj, (3) length of \nt, and (4) number of row of funcs are not equal.\n The ILFMM() will not proceed further.\n"))
  if(!is.null(covariates) & Yl!=nrow(cbind(X,X))) return(cat("Error: length of Y and number of rows of X is not equal.\n The ILFMM() will not proceed further.\n"))
  
  #Organizing f1(t) and f2(t)
  if(length(dim(f1_t))>0 ){
    f1_t<- f1_t
  } else
    if(length(f1_t)>0){
      f1_t<- matrix(f1_t, ncol=1)
    } else
    {
      f1_t<- matrix(rep(1, Yl), ncol=1)
    }
  f1_t<- f1_t[,which(apply(f1_t, 2, sd)>0)]
  f1_t<- cbind(1, f1_t)
  d1=ncol(f1_t)-1
  if(d1>5) warning("Only first 5 time components will be used", call. = FALSE)
  d1=min(d1,5)
  
  if(length(dim(f2_t))>0 ){
    f2_t<- f2_t
  } else
    if(length(f2_t)>0){
      f2_t<- matrix(f2_t, ncol=1)
    } else
    {
      f2_t<- matrix(rep(1, Yl), ncol=1)
    }
  f2_t<- f2_t[,which(apply(f2_t, 2, sd)>0)]
  f2_t<- cbind(1, f2_t)
  d2=ncol(f2_t)-1
  if(d2>5) warning("Only first 5 time components will be used", call. = FALSE)
  d2=min(d2,5)
  
  #Check 4: check in f1(t) and f2(t)
  if(dim(f1_t)[1]!=Yl) return(cat("Error: f1_t and Y are not compatible in dimension. \nThe ILFMM() will not proceed further.\n"))
  if(dim(f2_t)[1]!=Yl) return(cat("Error: f2_t and Y are not compatible in dimension. \nThe ILFMM() will not proceed further.\n"))
  
  #Sort the data by id and t & removal of missing and infinite observations
  tdata<- data.frame(id, t, Y, W,D, X, f1_t,f2_t)
  tdata<- tdata[which(apply(is.na(tdata), 1, sum)==0),]
  tdata<- tdata[which(apply(is.infinite(as.matrix(tdata)), 1, sum)==0),]
  tdata<- tdata[order(id, t), ]
  tdata<- tdata[!is.na(tdata$id) & !is.na(tdata$t) & !is.na(tdata$Y),]
  id<- tdata$id
  t<- tdata$t
  Y<- tdata$Y
  W.e<- dim(W)[2]+3; 
  W<- as.matrix(tdata[,c(4:W.e)])
  D.e<- dim(D)[2]+3; 
  D<- as.matrix(tdata[,c(4:D.e)])
  X.s<- W.e +D.e+ 1; 
  X.e<- dim(X)[2]+W.e+D.e; 
  X<- as.matrix(tdata[,c(X.s:X.e)])
  f_t.s<- X.e + 1; 
  f_t.e<- dim(f1_t)[2]+X.e; 
  f_t<- as.matrix(tdata[,c(f_t.s:f_t.e)])
  
  N<- length(unique(id))
  NT<- length(Y)
  
  #Check 5.1a.: Compatibility of Q1 and W matrix
  if(toupper(pentype)=='DECOMP' | toupper(pentype)=='DECOMPOSITION'){
    if (!comm.pen) 
    {
      if(ncol(Q1)!=(d+1)*ncol(W)) return(cat('Error: For different penalty, number of columns of Q1 need to be (d+1) \ntimes of number of columns of func1.\nThe ILFMM() will not proceed further.\n'))
    }
    if (comm.pen) 
    {
      if(ncol(Q1)!=ncol(W)) return(cat('Error: For common penalty, number of columns of pred for first time scale(i.e. W) and Q2 need to be equal.\nThe ILFMM() will not proceed further.\n'))
      Q1<- Q1
      for(i in 1:d) Q1<- cbind(Q1, Q1)
      Q1<- Q1
    }
    
    #Check 5.1b.: Compatibility of Q2 and D matrix
    if(toupper(pentype)=='DECOMP' | toupper(pentype)=='DECOMPOSITION'){
      if (!comm.pen) 
      {
        if(ncol(Q2)!=(d+1)*ncol(D)) return(cat('Error: For different penalty, number of columns of Q2 need to be (d+1) \ntimes of number of columns of func2.\nThe ILFMM() will not proceed further.\n'))
      }
      if (comm.pen) 
      {
        if(ncol(Q2)!=ncol(D)) return(cat('Error: For common penalty, number of columns of pred for secod time scale(i.e. D) and Q2 need to be equal.\nThe ILFMM() will not proceed further.\n'))
        Q2<- Q2
        for(i in 1:d) Q2<- cbind(Q2, Q2)
        Q2<- Q2
      }
      
    
    #Check 5.2: Removing rows containing missing and infinite values
    Q1<- Q1[which(apply(is.na(Q1), 1, sum)==0),]
    Q1<- Q1[which(apply(is.infinite(Q1), 1, sum)==0),]
    Q2<- Q2[which(apply(is.na(Q2), 1, sum)==0),]
    Q2<- Q2[which(apply(is.infinite(Q2), 1, sum)==0),]
    #Q<- matrix(Q, ncol=K)
    
    #Check 5.3: Singularity of Q1 and Q2 matrices
    Q1.eig<- abs(eigen(Q1 %*% t(Q1))$values)
    if(any(Q.eig<1e-12)) return(cat('Error: Q1 matrix is singular or near singular.\nThe ILFMM() will not proceed further.\n'))
    }
    
    Q2.eig<- abs(eigen(Q2 %*% t(Q2))$values)
    if(any(Q2.eig<1e-12)) return(cat('Error: Q2 matrix is singular or near singular.\nThe ILFMM() will not proceed further.\n'))
  }
  
  #Check 6.1a: Dimension of L1 matrix
  if(toupper(pentype)=='USER'){
    if (!comm.pen) 
    {
      if(ncol(L.user1)!=(d+1)*ncol(W)) return(cat('Error: For different penalty, number of columns of L1 need to be (d+1) \ntimes of number of columns of pred.\nThe ILFMM() will not proceed further.\n'))
    }
    if (comm.pen) 
    {
      if(ncol(L.user1)!=ncol(W)) return(cat('Error: For common penalty, number of columns of pred and L.user1 need to be equal.\nThe ILFMM() will not proceed further.\n'))
      L1.user<- L.user1
      for(i in 1:d) L1.user<- adiag(L1.user, L.user1)
      L.user1<- L1.user
    }
    
    L1<- L.user1
    
    
    #Check 6.1: Dimension of L matrix
    if(toupper(pentype)=='USER'){
      if (!comm.pen) 
      {
        if(ncol(L.user2)!=(d+1)*ncol(D)) return(cat('Error: For different penalty, number of columns of L2 need to be (d+1) \ntimes of number of columns of pred.\nThe ILFMM() will not proceed further.\n'))
      }
      if (comm.pen) 
      {
        if(ncol(L.user2)!=ncol(D)) return(cat('Error: For common penalty, number of columns of pred and L.user2 need to be equal.\nThe ILFMM() will not proceed further.\n'))
        L2.user<- L.user2
        for(i in 1:d) L2.user<- adiag(L2.user, L.user2)
        L.user2<- L2.user
      }
      
      L2<- L.user2
      
    #Check 6.2: Removing rows containing missing and infinite values
    L1<- L1[which(apply(is.na(L1), 1, sum)==0),]
    L1<- L[which(apply(is.infinite(L1), 1, sum)==0),]
    
    L2<- L2[which(apply(is.na(L2), 1, sum)==0),]
    L2<- L2[which(apply(is.infinite(L2), 1, sum)==0),]
    
    #Check 6.3: Singularity of L1'L1 and L2'L2 matrices
    L1L1<- t(L1)%*%L1
    L1L1.eig<- abs(eigen(L1L1 %*% t(L1L1))$values)
    if(any(L1L1.eig<1e-12)) return(cat("Error: L1'L1 matrix is singular or near singular.\nThe ILFMM() will not proceed further.\n"))
    }
    
    L2L2<- t(L2)%*%L2
    L2L2.eig<- abs(eigen(L2L2 %*% t(L2L2))$values)
    if(any(L2L2.eig<1e-12)) return(cat("Error: L2'L2 matrix is singular or near singular.\nThe ILFMM() will not proceed further.\n"))
  }
  
  #Generate L matrix for D2 penalty
  if(toupper(pentype)=='D2'){
    Left<- cbind(diag(rep(1,K-2)),rep(0,K-2),rep(0,K-2))
    Middle<- cbind(rep(0,K-2),diag(rep(-2,K-2)),rep(0,K-2))
    Right<- cbind(rep(0,K-2),rep(0,K-2),diag(rep(1,K-2)))
    D.2<- rbind(Left+Middle+Right, c(rep(0, K-2), 1, -2), c(rep(0, K-2), 0, 1))
  }
  
  #Generate W1 matrix (This is referred as in the paper W matrix)
  for(i in 0:d){
    if (i==0) W1<-data.matrix(W)*f_t[,(i+1)] 
    if (i>0) W1<- cbind(W1, data.matrix(W)*f_t[,(i+1)])
  }
  
  #Generate W* matrix
  for(i in 0:d){
    if(toupper(pentype)=='DECOMP' | toupper(pentype)=='DECOMPOSITION'){
      tQ1<- Q1[,(i*K+1):((i+1)*K)] 
      tP_Q1 <- t(tQ1) %*% solve(tQ1 %*% t(tQ1)) %*% tQ1
      tL_PEER1<- a*(diag(K)- tP_Q1) + 1*tP_Q1
      rm(tQ1); rm(tP_Q1)
    } else
      if(toupper(pentype)=='RIDGE'){
        tL_PEER1<- diag(K)
      } else
        if(toupper(pentype)=='D2'){
          tL_PEER1<- D.2
        } else
          if(toupper(pentype)=='USER'){
            tL_PEER1<- L1[(i*K+1):((i+1)*K),(i*K+1):((i+1)*K)]
          } 
    
    J1<- tL_PEER1
    J1.eig <- eigen(J1)
    J1.sqrt <- J1.eig$vectors %*% diag(sqrt(J1.eig$values)) %*% solve(J1.eig$vectors)
    v1<- diag(K)
    if(K>N) v1<-  svd((data.matrix(W)*f_t[,(i+1)])%*% solve(tL_PEER1))$v1
    assign(paste('W',i+1, '_PEER', sep=''), 
           (data.matrix(W)*f_t[,(i+1)])%*% solve(J1.sqrt) %*% v1)
    rm(tL_PEER1) ; rm(v1); rm(J1.sqrt); rm(J1); rm(J1.eig)
  }
  
  
  #Generate D1 matrix (This is referred as in the paper D matrix)
  for(i in 0:d){
    if (i==0) D1<-data.matrix(D)*f_t[,(i+1)] 
    if (i>0) D1<- cbind(D1, data.matrix(D)*f_t[,(i+1)])
  }
  
  #Generate D* matrix
  for(i in 0:d){
    if(toupper(pentype)=='DECOMP' | toupper(pentype)=='DECOMPOSITION'){
      tQ2<- Q2[,(i*K+1):((i+1)*K)] 
      tP_Q2 <- t(tQ2) %*% solve(tQ2 %*% t(tQ2)) %*% tQ2
      tL_PEER2<- c*(diag(K)- tP_Q2) + 1*tP_Q2
      rm(tQ2); rm(tP_Q2)
    } else
      if(toupper(pentype)=='RIDGE'){
        tL_PEER2<- diag(K)
      } else
        if(toupper(pentype)=='D2'){
          tL_PEER2<- D.2
        } else
          if(toupper(pentype)=='USER'){
            tL_PEER2<- L2[(i*K+1):((i+1)*K),(i*K+1):((i+1)*K)]
          } 
    
    J2<- tL_PEER2
    J2.eig <- eigen(J2)
    J2.sqrt <- J2.eig$vectors %*% diag(sqrt(J2.eig$values)) %*% solve(J2.eig$vectors)
    v2<- diag(K)
    if(K>N) v2<-  svd((data.matrix(D)*f_t[,(i+1)])%*% solve(tL_PEER))$v2
    assign(paste('D',i+1, '_PEER', sep=''), 
           (data.matrix(D)*f_t[,(i+1)])%*% solve(J2.sqrt) %*% v2)
    rm(tL_PEER2) ; rm(v2); rm(J2.sqrt); rm(J2); rm(J2.eig)
  }
  #Generate Z
  id.bd1<- factor(rep(1, NT))
  ni<- tapply(id, id, length)
  for(i in 1:N){
    if (i==1) Z<- matrix(1, nrow=ni[i])
    if (i>1) Z<- adiag(Z, matrix(1, nrow=ni[i]))
  }
  
  #Input for random argument of lme function 
  for(i in 0:d) assign(paste('pd', i+1, sep=''), 
                       pdIdent(form=as.formula(paste('~W+D', i+1, '_PEER -1', sep=''))))
  pdid<- pdIdent(~Z-1)
  
  if(d==0) tXX<- pdBlocked(list(pd1, pdid))
  if(d==1) tXX<- pdBlocked(list(pd1, pd2, pdid))
  if(d==2) tXX<- pdBlocked(list(pd1, pd2, pd3, pdid))
  if(d==3) tXX<- pdBlocked(list(pd1, pd2, pd3, pd4, pdid))
  if(d==4) tXX<- pdBlocked(list(pd1, pd2, pd3, pd4, pd5, pdid))
  if(d==5) tXX<- pdBlocked(list(pd1, pd2, pd3, pd4, pd5, pd6, pdid))
  
  #Fitting the model
  subject=as.matrix(rep(1:N,each=20))
  one <- rep(1,nrow(Y))
  
  out_ILFMM<- lme(fixed=Y~X+W+D-1, random=list(one=pdIdent(~subject-1)), ... ) 
  
  cat('The fit is successful.\n')
  
  
  
  #Gamma_PEER<-matrix(out_PEER$coeff$random$id.bd1, ncol=1)
  Coeff<-matrix(out_ILPEER$coeff$fixed, ncol=1)
  
  for(i in 0:d) {
    if(toupper(pentype)=='DECOMP' | toupper(pentype)=='DECOMPOSITION'){
      tQ1<- Q1[,(i*K+1):((i+1)*K)] 
      tP_Q1 <- t(tQ1) %*% solve(tQ1 %*% t(tQ1)) %*% tQ1
      tL_PEER1<- a*(diag(K)- tP_Q1) + 1*tP_Q1
      rm(tQ1); rm(tP_Q1)
      
      tQ2<- Q2[,(i*K+1):((i+1)*K)] 
      tP_Q2 <- t(tQ2) %*% solve(tQ2 %*% t(tQ2)) %*% tQ2
      tL_PEER2<- c*(diag(K)- tP_Q2) + 1*tP_Q2
      rm(tQ2); rm(tP_Q2)
    } else
      if(toupper(pentype)=='RIDGE'){
        tL_PEER1<- diag(K)
        tL_PEER2<- diag(K)
      } else
        if(toupper(pentype)=='D2'){
          tL_PEER1<- D.2
          tL_PEER2<- D.2
        } else
          if(toupper(pentype)=='USER'){
            tL_PEER1<- L1[(i*K+1):((i+1)*K),(i*K+1):((i+1)*K)]
            tL_PEER2<- L2[(i*K+1):((i+1)*K),(i*K+1):((i+1)*K)]
          } 
    
    v1<- diag(K)
    v2<-diag(K)
    if(K>N) v1<-  svd((data.matrix(W)*f_t[,(i+1)])%*% solve(tL_PEER1))$v1
    if(K>N) v2<-  svd((data.matrix(D)*f_t[,(i+1)])%*% solve(tL_PEER1))$v2
    K1<- ncol(v1)
    
    r1=K1+2
    r2=K1+3
    r3=2*K1+2
    r4=2*K1+3
    r5=3*K1+2
    r6=3*K1+3
    r7=4*K1+2
    
    beta.hat.ILPEER<-matrix(Coeff[1:2], ncol=1)
    Gamma1.ILPEER.hat<- matrix(Coeff[3:r1], ncol=1)
    Gamma2.ILPEER.hat<- matrix(Coeff[r2:r3], ncol=1)
    Eta1.ILPEER.hat<- matrix(Coeff[r4:r5], ncol=1)
    Eta2.ILPEER.hat<- matrix(Coeff[r6:r7], ncol=1)
    
    #tGamma.PEER.hat<- matrix(Gamma_PEER[(i*r+1):((i+1)*r)],ncol=1)
    J1<- tL_PEER1
    J1.eig <- eigen(J1)
    J1.sqrt <- J1.eig$vectors %*% diag(sqrt(J1.eig$values)) %*% solve(J1.eig$vectors)
    tGamma1Hat <- solve(J1.sqrt) %*% v1 %*%Gamma1.ILPEER.hat
    tGamma2Hat <- solve(J1.sqrt) %*% v1 %*%Gamma2.ILPEER.hat
    
    J2<- tL_PEER2
    J2.eig <- eigen(J2)
    J2.sqrt <- J2.eig$vectors %*% diag(sqrt(J2.eig$values)) %*% solve(J2.eig$vectors)
    tEta1Hat <- solve(J2.sqrt) %*% v2 %*%Eta1.ILPEER.hat
    tEta2Hat <- solve(J2.sqrt) %*% v2 %*%Eta2.ILPEER.hat
    
    if(i==0) 
      {
      Gamma1Hat<- matrix(tGamma1Hat ,ncol=1)
      Gamma2Hat<- matrix(tGamma2Hat ,ncol=1)
      Eta1Hat<- matrix(tEta1Hat ,ncol=1)
      Eta2Hat<- matrix(tEta2Hat ,ncol=1)
    }
      
    if(i>0) 
      {
      Gamma1Hat<- cbind(Gamma1Hat, tGamma1Hat)
      Gamma2Hat<- cbind(Gamma2Hat, tGamma2Hat)
      Eta1Hat<- cbind(Eta1Hat, tEta1Hat)
      Eta2Hat<- cbind(Eta2Hat, tEta2Hat)
    }
    rm(tL_PEER1); rm(v1); rm(tGamma1Hat); rm(J1.sqrt); rm(J1); rm(J1.eig)
    rm(tL_PEER2); rm(v2); rm(tGamma2Hat); rm(J2.sqrt); rm(J2); rm(J2.eig)
    rm(tEta1Hat);
    rm(tEta2Hat);
  }
  colnames(Gamma1Hat)<- paste('Gamma1', 0:d, sep='')
  colnames(Gamma2Hat)<- paste('Gamma2', 0:d, sep='')
  colnames(Eta1Hat)<- paste('Eta1', 0:d, sep='')
  colnames(Eta2Hat)<- paste('Eta2', 0:d, sep='')
  
  beta.hat<-matrix(Coeff[1:2], ncol=1)
  names(beta.hat)<- c('Intercept', 'slope', colnames(covariates))
  
  fitted.vals<- summary(out_ILFMM)$fitted
  logLik<- summary(out_ILFMM)$logLik
  #AIC<- summary(out_PEER)$AIC
  #BIC<- summary(out_PEER)$BIC
  
  
  
  tVarCorr<- nlme:::VarCorr(out_ILFMM, rdig=4)[,2]
  for(i in 0:d) assign(paste('lambda', i, sep=''), 
                       1/ as.numeric(unique(tVarCorr[(i*r+1):((i+1)*r)])))
  for(i in 0:d)
  {
    tLambda<- 1/ as.numeric(unique(tVarCorr[(i*r+1):((i+1)*r)]))
    if(i==0) lambda<- tLambda
    if(i>0) lambda<- c(lambda, tLambda)
  }
  names(lambda)<- paste('lambda', 0:d, sep='')
  
  print(lambda)
  
  
  sd_int.est<- as.numeric(unique(tVarCorr[((d+1)*r+1):((d+1)*r+N)]))
  sigma<- out_ILPEER$sigma
  
  ###---- Standard Error
  Sigma.u<- sd_int.est
  sigma.e<- sigma
  
  for(i in 0:d)
  {
    if(toupper(pentype)=='DECOMP' | toupper(pentype)=='DECOMPOSITION'){
      tQ1<- Q1[,(i*K+1):((i+1)*K)] 
      tP_Q1 <- t(tQ1) %*% solve(tQ1 %*% t(tQ1)) %*% tQ1
      tL_PEER1<- a*(diag(K)- tP_Q1) + 1*tP_Q1
      rm(tQ1); rm(tP_Q1)
      
      tQ2<- Q2[,(i*K+1):((i+1)*K)] 
      tP_Q2 <- t(tQ2) %*% solve(tQ2 %*% t(tQ2)) %*% tQ2
      tL_PEER2<- c*(diag(K)- tP_Q2) + 1*tP_Q2
      rm(tQ2); rm(tP_Q2)
    } else
      if(toupper(pentype)=='RIDGE'){
        tL_PEER1<- diag(K)
        tL_PEER2<- diag(K)
      } else
        if(toupper(pentype)=='D2'){
          tL_PEER1<- D.2
          tL_PEER2<- D.2
        } else
          if(toupper(pentype)=='USER'){
            tL_PEER1<- L1[(i*K+1):((i+1)*K),(i*K+1):((i+1)*K)]
            tL_PEER2<- L2[(i*K+1):((i+1)*K),(i*K+1):((i+1)*K)]
          } 
    tsigma<- as.numeric(unique(tVarCorr[(i*r+1):((i+1)*r)]))
    if(i==0)
      {
      L1L1.inv<- tsigma^2*solve(t(tL_PEER1)%*%tL_PEER1) 
      L2L2.inv<- tsigma^2*solve(t(tL_PEER2)%*%tL_PEER2) 
    }
    if(i>0) 
      {
      L1L1.inv<- adiag(L1L1.inv, tsigma^2*solve(t(tL_PEER1)%*%tL_PEER1))
      L2L2.inv<- adiag(L2L2.inv, tsigma^2*solve(t(tL_PEER2)%*%tL_PEER2))
    }
    if(i==0)
      {
      L1L1<- 1/(tsigma^2)*(t(tL_PEER1)%*%tL_PEER1) 
      L2L2<- 1/(tsigma^2)*(t(tL_PEER2)%*%tL_PEER2)
    }
    if(i>0) 
      {
      L1L1<- adiag(L1L1, 1/(tsigma^2)*(t(tL_PEER1)%*%tL_PEER1))
      L2L2<- adiag(L2L2, 1/(tsigma^2)*(t(tL_PEER2)%*%tL_PEER2))
    }
    
    rm(tsigma); rm(tL_PEER1)
    rm(tL_PEER2)
  }
  
  
  rand.var<- Sigma.u^2*(Z%*%t(Z))
  
  if(!se){
    status<- 0
    ret <- list(out_ILPEER, beta.hat,  fitted.vals,  
                Gamma1Hat, Gamma2Hat,Eta1Hat,Eta2Hat,logLik, lobj, 
                lambda, N, K, NT, Sigma.u, sigma.e, d, status)
    names(ret)<- c("fit", "BetaHat", "fitted.vals",
                   "Gamma1Hat", "Gamma2Hat","Eta1Hat","Eta2Hat","logLik", "lpeerobj",
                   "lambda", "N", "K", "TotalObs", "Sigma.u", "sigma", "d", "status")
    return(ret)
  }
  
  V.1<- W1%*%L1L1.inv%*%t(W1)+rand.var+sigma.e^2*diag(rep(1, NT))
  V.2<- D1%*%L2L2.inv%*%t(D1)+rand.var+sigma.e^2*diag(rep(1, NT))
  V<- rand.var+sigma.e^2*diag(rep(1, NT))
  Vconb<-cbind(V.1,V.2)
  X.Vinv.X.inv<- solve(t(X)%*%solve(Vconb)%*%X)
  X.Vinv.Y<- t(X)%*%solve(Vcob)%*%Y
  
  se.Beta<- sqrt(diag(X.Vinv.X.inv%*%t(X)%*%solve(Vconb)%*%V%*%solve(Vconb)%*%X%*%X.Vinv.X.inv))
  names(se.Beta)<- c('Intercept', colnames(covariates))
  Beta<- cbind(beta.hat, se.Beta)
  
  p1<- L1L1.inv%*%t(W1)%*%solve(V.1)
  p2<- L2L2.inv%*%t(D1)%*%solve(V.2)
  p3<- V.1 - X%*%X.Vinv.X.inv%*%t(X)
  p4<- V.2 - X%*%X.Vinv.X.inv%*%t(X)
  p5<- solve(V.1)
  p6<- solve(V.2)
  pGamma<- p1%*%p3%*%p5
  pEta<-p2%*%p4%*%p6
  
  SE.gamma<- sqrt(diag(pGamma%*%V%*%t(pGamma)))
  SE.eta<- sqrt(diag(pEta%*%V%*%t(pEta)))
  for(i in 0:d)
  {
    if(i==0) 
      {
      se.Gamma<- matrix(SE.gamma[(i*K+1):((i+1)*K)], ncol=1)
      se.Eta<- matrix(SE.eta[(i*K+1):((i+1)*K)], ncol=1)
    }
    if(i>0) 
      {
      se.Gamma<- cbind(se.Gamma, SE.gamma[(i*K+1):((i+1)*K)])
      se.Eta<- cbind(se.Eta, SE.eta[(i*K+1):((i+1)*K)])
    }
  } 
  colnames(se.Gamma)<- paste('Gamma', 0:d, sep='')
  colnames(se.Eta)<- paste('Eta', 0:d, sep='')
  #------- Calcualtion of AIC
  C=cbind(X, W1,D1)
  Xcol<- ncol(X)
  XCol=adiag(matrix(rep(0,Xcol^2),Xcol,Xcol), L1L1,L2L2)
  Smooth=solve(t(C)%*%solve(V)%*%C + XCol)%*%(t(C)%*%solve(V))
  df<- sum(diag(Smooth))
  resid<- residuals(out_ILPEER)
  RSS<- sum(resid^2)
  AIC<- log(RSS)+2*df/NT
  
  
  status<- 1
  ret <- list(out_ILPEER,beta.hat,  se.Beta, Beta, fitted.vals,  
              Gamma1Hat,Gamma2Hat,Eta1Hat,Eta2Hat, se.Gamma,se.Eta, logLik, V.1, V.2, lobj,
              V, lambda, N, K, NT, Sigma.u, sigma.e, d, status, 
              X, W1, D1,L1L1,L2L2, C, XCol, resid, RSS,
              Smooth, df, AIC)
  names(ret)<- c("fit", "BetaHat", "se.Beta", "Beta", "fitted.vals",
                 "Gamma1Hat", "Gamma2Hat","Eta1Hat","Eta2Hat","se.Gamma", "se.Eta","logLik", "V1", "V2","lpeerobj",
                 "V", "lambda", "N", "K", "TotalObs", "Sigma.u", "sigma", "d", "status",
                 "X", "W","D" ,"L1L1","L1L1", "C", "XCol", "resid", "RSS",
                 "Smooth", "df", "AIC")
  ret
}



### Function to plot estimated functional ceofficients
plot.ILFMM<- function(lfit, conf=0.95, ...){
  if(!exists("lfit")) return (cat("Error: The value specified in lfit argument is not an ILFMM object.\n"))
  if(!is.list(lfit)) return (cat("Error: The value specified in lfit argument is not an ILFMM object.\n"))
  if(is.na(match("ILFMMobj", names(lfit)))) return (cat("Error: The value specified in lfit argument is not an ILFMM object.\n"))
  if(conf>0.99 | conf<0.70) return (cat("Error: Confidence level should be within 0.70 and 0.99\n"))
  d<- lfit$d
  status<- lfit$status
  if(d==0) par(mfrow=c(1,1))
  if(d==1) par(mfrow=c(1,2))
  if(d>1) par(mfrow=c(2,2))
  for(i in 0:d)
  {
    estGamma1<- lfit$Gamma1Hat[,(i+1)]
    estGamma2<- lfit$Gamma2Hat[,(i+1)]
    estEta1<- lfit$Eta1Hat[,(i+1)]
    estEta2<- lfit$Eta2Hat[,(i+1)]
    if(status==0) 
      {
      matplot(estGamma1, type='l', main=paste('gamma1', i, sep=''), ...)
      matplot(estGamma2, type='l', main=paste('gamma2', i, sep=''), ...)
      matplot(estEta1, type='l', main=paste('eta1', i, sep=''), ...)
      matplot(estEta2, type='l', main=paste('eta2', i, sep=''), ...)
    }
    if(status==1){
      llGamma<- lfit$Gamma1Hat[,(i+1)] - qnorm(0.5+conf/2)*lfit$se.Gamma[,(i+1)]
      ulGamma<- lfit$Gamma1Hat[,(i+1)] + qnorm(0.5+conf/2)*lfit$se.Gamma[,(i+1)]
      llEta<- lfit$Eta1Hat[,(i+1)] - qnorm(0.5+conf/2)*lfit$se.Eta[,(i+1)]
      ulEta<- lfit$Eta1Hat[,(i+1)] + qnorm(0.5+conf/2)*lfit$se.Eta[,(i+1)]
      matplot(estGamma1, type='l', ylim=range(est, llGamma, ulGammma), 
              main=paste('gamma', i, sep=''), ...)
      matplot(estEta1, type='l', ylim=range(est, llEta, ulEta), 
              main=paste('eta', i, sep=''), ...)
      matplot(llGamma, type='l', add=T, lty=2, col=2)
      matplot(ulGamma, type='l', add=T, lty=2, col=2)
      matplot(llEta, type='l', add=T, lty=2, col=2)
      matplot(ulEta, type='l', add=T, lty=2, col=2)
    }
    abline(h=0)
  }
}

setwd('/Users/mostafa/Mostafa/UNC/Research/Proposal/OurRpackage/ILFMM')
document()

install.packages(‘devtools’)
library(devtools)
install_github('Mostafa2020-bit/NewRepository')

