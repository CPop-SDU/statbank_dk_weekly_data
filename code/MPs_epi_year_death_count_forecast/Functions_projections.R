rm(list = ls())

library(compositions)


######### Life table functions #########

### function starting from mx

#mx is a matrix of death rates by age (rows) and year(columns)
#r is the radix of the life table, generally 1 or 100,000
#nx is the width of the age-interval

LifeT_mx<-function(mx, r, nx){
  
  ax<-matrix(nx/2, nrow(mx), ncol(mx))
  
  #Adjust ax
  nY1<-which(mx[1,]< 0.01724)
  ax[1,nY1]<-0.14903 - 2.05527*mx[1,nY1]
  
  nY2<-which(mx[1,]>= 0.01724 & mx[1,]< 0.06891)
  ax[1,nY2]<-0.04667 + 3.88089*mx[1,nY2]
  
  nY3<-which(mx[1,]>= 0.06891)
  ax[1,nY3]<-0.31411
  
  ax[nrow(ax),]<-1/mx[nrow(ax),]
  
  
  #qx
  qx<- mx*nx/(1+(nx-ax)*mx)
  qx[nrow(qx),]<-1
  qx[qx>1]<-1
  
  #px
  px<-1-qx
  
  #lx
  lx<-rbind(rep(r, ncol(mx)),
            r*apply(px[1:(nrow(mx)-1),], 2, FUN = "cumprod"))
  
  #dx
  dx<-qx*lx
  
  #Lx
  Lx<- lx*nx - dx*(nx-ax)
  
  #Tx
  Tx<- apply(Lx[nrow(Lx):1,], 2, "cumsum")[nrow(Lx):1,]
  
  #ex
  ex<- Tx/lx
  
  colnames(qx)<- colnames(dx)<-colnames(lx)<-colnames(Lx)<-colnames(Tx)<-colnames(ex)<-colnames(mx)
  rownames(qx)<- rownames(dx)<-rownames(lx)<-rownames(Lx)<-rownames(Tx)<-rownames(ex)<-rownames(mx)
  return(list(mx=mx, qx=qx, dx=dx, lx=lx, Lx=Lx, Tx=Tx, ex=ex))
}



###function starting from dx

#dx is a matrix of life table deaths by age (rows) and year(columns)
#r is the radix of the life table, generally 1 or 100,000
#nx is the width of the age-interval

LifeT_dx<-function(dx, r, nx){
  
  ax<-matrix(nx/2, nrow(dx), ncol(dx))
  
  #dx
  dx<-dx*r
  
  #lx
  lx<-apply(dx[nrow(dx):1,], 2, cumsum)[nrow(dx):1,]
  
  #qx
  qx<-dx/lx
  qx[qx>1]<-1
  
  #px
  px<-1-qx
  
  #mx first estimate
  mx1<-qx/(nx-((nx-ax)*qx))
  
  #Adjust ax
  nY1<-which(mx1[1,]< 0.01724)
  ax[1,nY1]<-0.14903 - 2.05527*mx1[1,nY1]
  
  nY2<-which(mx1[1,]>= 0.01724 & mx1[1,]< 0.06891)
  ax[1,nY2]<-0.04667 + 3.88089*mx1[1,nY2]
  
  nY3<-which(mx1[1,]>= 0.06891)
  ax[1,nY3]<-0.31411
  
  ax[nrow(ax),]<-1/mx1[nrow(ax)-1,]
  
  #mx
  mx<-qx/(nx-((nx-ax)*qx))
  
  #Lx
  Lx<- lx*nx - dx*(nx-ax)
  
  #Tx
  Tx<- apply(Lx[nrow(Lx):1,], 2, "cumsum")[nrow(Lx):1,]
  
  #ex
  ex<- Tx/lx
  
  colnames(mx)<- colnames(qx)<-colnames(lx)<-colnames(Lx)<-colnames(Tx)<-colnames(ex)<-colnames(dx)
  rownames(mx)<- rownames(qx)<-rownames(lx)<-rownames(Lx)<-rownames(Tx)<-rownames(ex)<-rownames(dx)
  return(list(mx=mx, qx=qx, dx=dx, lx=lx, Lx=Lx, Tx=Tx, ex=ex))
}




######### Forecast with CoDA (mortality only) #########


### Forecast mortality using the CoDA approach
# dx is a matrix of life table deaths by age (rows) and year(columns)
# h is the forecast horizon
# nsim is the number of simulation to be used to calculate the prediction intervals

CoDA.fun<-function(dx, h, nsim){
  
  dx2<-t(dx)
  ax<-dx2[nrow(dx2),]
  tdx<-sweep(acomp(dx2), 2, acomp(ax), FUN="-")
  cdx<- clr(acomp(tdx))
  
  svd.dx<-svd(cdx)
  kt<-svd.dx$u[,1]
  bx<-svd.dx$v[,1]
  s<-svd.dx$d
  EV<- cumsum((s)^2/sum((s)^2))[1]
  
  model<- Arima(kt, order=c(0,1,0), include.drift = T)
  fcst<-forecast(model, h=h)
  kt.fit<-fcst$fitted
  kt.fcst<-fcst$mean
  bst<-replicate(nsim, simulate(model, nsim=h, bootstrap=T), simplify = "array")
  
  cdx.ff<- c(c(kt, kt.fcst)*s[1]) %*% t(bx)
  tdx.ff<-clrInv(cdx.ff)
  
  dx.ff<-sweep(tdx.ff, 2, acomp(ax), FUN="+")
  
  cdx.bst<- sapply(1:nsim, function(j) (c(bst[,j]*s[1]) %*% t(bx)), simplify = "array")
  rm(bst)
  cdx.PI<-apply(cdx.bst, c(1,2),
                function(cdx.bst) quantile(cdx.bst, prob=c(0.025, 0.5, 0.975), type=8))
  rm(cdx.bst)
  dx.lw<-sweep(clrInv(cdx.PI[1,,]), 2, acomp(ax), FUN="+")
  dx.up<-sweep(clrInv(cdx.PI[3,,]), 2, acomp(ax), FUN="+")
  
  
  LT<-LifeT_dx(unclass(t(dx.ff[(nrow(dx2)+1):nrow(dx.ff),])), r, nx)
  LT.up<-LifeT_dx(unclass(t(dx.up)), r, nx)
  LT.lw<-LifeT_dx(unclass(t(dx.lw)), r, nx)
  
  
  output<-list(kt=kt, bx=bx, ax=ax,  EV=EV,
               LT=LT, LT.up=LT.up, LT.lw=LT.lw)
  return(output)
  
}


### Forecast mortality using the CoDA-coherent approach
# dx is a matrix of life table deaths by age (rows) and year(columns)
# dx.ref is a matrix of life table deaths by age (rows) and year(columns) for the reference population
# h is the forecast horizon
# nsim is the number of simulation to be used to calculate the prediction intervals for the reference population
# nsim is the number of simulation to be used to calculate the prediction intervals for the population-specific deviation


CoDAC.fun<-function(dx, dx.ref, h, nsim1, nsim2){
  
  #Average
  dx2.ref<-t(dx.ref)
  ax.ref<-acomp(dx2.ref[nrow(dx2.ref),])
  tdx.ref<-sweep(acomp(dx2.ref), 2, ax.ref, FUN="-")
  cdx.ref<- clr(acomp(tdx.ref))
  
  svd.ref<-svd(cdx.ref)
  Kt<-svd.ref$u[,1]
  Bx<-svd.ref$v[,1]
  S<-svd.ref$d
  EV<- cumsum((S)^2/sum((S)^2))[1]
  
  model.ref<- Arima(Kt, order=c(0,1,0), include.drift = T)
  fcst.ref<-forecast(model.ref, h=h)
  Kt.fit<-fcst.ref$fitted
  Kt.fcst<-fcst.ref$mean
  bst.ref<-replicate(nsim1, simulate(model.ref, nsim=h, bootstrap=T), simplify = "array")
  
  cdx.ff.ref<- c(c(Kt, Kt.fcst)*S[1]) %*% t(Bx)
  tdx.ff.ref<-clrInv(cdx.ff.ref)
  
  cdx.bst.ref<- sapply(1:nsim1, function(j) (c(bst.ref[,j]*S[1]) %*% t(Bx)), simplify = "array")
  rm(bst.ref)
  
  #deviation
  dx2<-t(dx)
  ax<-acomp(dx2[nrow(dx2),])
  tdx<-sweep(acomp(dx2), 2, acomp(ax), FUN="-")
  tdx2<-tdx-tdx.ff.ref[1:nrow(dx2),]
  cdx<- clr(acomp(tdx2))
  
  svd.dx<-svd(cdx)
  kt<-svd.dx$u[,1]
  bx<-svd.dx$v[,1]
  s<-svd.dx$d
  EV<- cumsum((s)^2/sum((s)^2))[1]
  
  model<- Arima(Kt, order=c(1,0,0), include.drift = F, include.mean = T)
  fcst<-forecast(model, h=h)
  kt.fit<-fcst$fitted
  kt.fcst<-fcst$mean
  bst<-replicate(nsim2, simulate(model, nsim=h, bootstrap=T), simplify = "array")
  
  cdx.ff<- c(c(kt, kt.fcst)*s[1]) %*% t(bx)
  tdx.ff2<-clrInv(cdx.ff)
  tdx.ff<-tdx.ff2+tdx.ff.ref
  
  dx.ff<-sweep(tdx.ff, 2, ax, FUN="+")
  
  cdx.bst2<- sapply(1:nsim2, function(j) (c(bst[,j]*s[1]) %*% t(bx)), simplify = "array")
  rm(bst)
  
  kseq<-jseq<-1:nsim1
  cdx.bst<-array(0, dim = c(dim(cdx.bst2)[1:2], nsim1*nsim2))
  for(i in 1:nsim2){
    cdx.bst[,,kseq]<-sapply(jseq, function(j) cdx.bst2[ , ,i] + cdx.bst.ref[ , ,j], simplify = "array")
    kseq<-kseq+nsim1
  }
  rm(cdx.bst2)
  rm(cdx.bst.ref)
  cdx.PI<-apply(cdx.bst, c(1,2),
                function(cdx.bst) quantile(cdx.bst, prob=c(0.025, 0.5, 0.975), type=8))
  rm(cdx.bst)
  dx.lw<-sweep(clrInv(cdx.PI[1,,]), 2, acomp(ax), FUN="+")
  dx.up<-sweep(clrInv(cdx.PI[3,,]), 2, acomp(ax), FUN="+")
  
  
  LT<-LifeT_dx(unclass(t(dx.ff[(nrow(dx2)+1):nrow(dx.ff),])), r, nx)
  LT.up<-LifeT_dx(unclass(t(dx.up)), r, nx)
  LT.lw<-LifeT_dx(unclass(t(dx.lw)), r, nx)
  
  
  output<-list(kt=kt, bx=bx, ax=ax,  EV=EV, 
               LT=LT, LT.up=LT.up, LT.lw=LT.lw)
  return(output)
  
}





######### Forecast with the LC model (mortality and fertility) #########

#rx is a matrix of rates (mortality or fertility) by age (rows) and year(columns)
# h is the forecast horizon
# nsim is the number of simulation to be used to calculate the prediction intervals

LC.fun<-function(rx, h, nsim){
  
  #log-transformed the rates
  trx<-log(rx)
  
  #Calculate the age-specific mean
  ax<-trx[,ncol(trx)]
  
  #Center the matrix
  crx<- sweep(trx, 1, ax, FUN="-")
  
  #Calculate the SVD
  svd.rx<-svd(crx)
  kt<-svd.rx$v[,1]
  bx<-svd.rx$u[,1]
  s<-svd.rx$d
  EV<- cumsum((s)^2/sum((s)^2))[1]
  
  #Forecast kt
  model<- Arima(kt, order=c(0,1,0), include.drift = T)
  fcst<-forecast(model, h=h)
  kt.fit<-fcst$fitted
  kt.fcst<-fcst$mean
  bst<-replicate(nsim, simulate(model, nsim=h, bootstrap=T), simplify = "array")
  
  #Reconstructing the matrix
  crx.ff<- c(bx) %*% t(c(kt, kt.fcst)*s[1])
  rx.ff<-exp(sweep(crx.ff, 1, ax, FUN="+"))
  
  #Reconstructing the PI
  cdx.bst<- sapply(1:nsim, function(j) (c(bx) %*% t(bst[,j]*s[1])), simplify = "array")
  rm(bst)
  cdx.PI<-apply(cdx.bst, c(1,2),
                function(cdx.bst) quantile(cdx.bst, prob=c(0.025, 0.5, 0.975), type=8))
  rm(cdx.bst)
  rx.lw<-exp(sweep(cdx.PI[1,,], 1, ax, FUN="+"))
  rx.up<-exp(sweep(cdx.PI[3,,], 1, ax, FUN="+"))
  
  output<-list(kt=kt, bx=bx, ax=ax,
               rx.fcst=rx.ff[,(ncol(rx)+1):ncol(rx.ff)], 
               rx.fit=rx.ff[,1:ncol(rx)], 
               rx.lw = rx.lw,
               rx.up = rx.up)
  return(output)
  
}



######### Cohort-component forecasts, based on the Leslie matrix #########

#Lx is the life table exposure forecast for females (LxF) and males (LxM)
# Tx is the cumulated exposure from age age
# lx is the number of survivors to age x in the life table
# Fx are the forecast fertility rates
# pop are the population counts by age at the last obversed year


CC_fun <- function(LxF, TxF, lxF,
                  LxM, TxM, lxM, 
                  Fx, 
                  popF, popM){
  nh<-ncol(LxF)
  nr<-nrow(LxF)
  projF<-matrix(0, nr, nh)
  projM<-matrix(0, nr, nh)
  
  ageF<-as.numeric(rownames(Fx))
  
  for(i in 1:nh){
    
    matF<-matrix(0, nr, nr)
    matM<-matrix(0, nr, nr)
    
    survF<-LxF[2:nr,i]/LxF[1:(nr-1),i]
    survF[survF>1]<-0.5
    matF[row(matF)==col(matF)+1] <- survF
    
    survM<-LxM[2:nr,i]/LxM[1:(nr-1),i]
    survM[survM>1]<-0.5
    matM[row(matM)==col(matM)+1] <- survM
    
    ASFR<-rep(0, nr)
    ASFR[ageF+1]<-Fx[,i]
    
    
    constF<- (1/(1+1.05))*(LxF[2,i]/(2*lxF[1,i]))
    constM<- (1.05/(1+1.05))*(LxM[2,i]/(2*lxM[1,i]))
    
    firstrowF <- constF*(ASFR[1:(nr-1)]+(ASFR[2:nr]*survF[1:(nr-1)]))
    firstrowM <- constM*(ASFR[1:(nr-1)]+(ASFR[2:nr]*survF[1:(nr-1)]))
    
    matF[1,] <- c(firstrowF,0)
    matM[1,] <- c(firstrowM,0)
    
    
    matF[nr,nr]<-(TxF[nr,i]/TxF[nr-1,i])
    matM[nr,nr]<-(TxM[nr,i]/TxM[nr-1,i])
    
    if(i==1){
      projF[,i] <- matF%*%popF}
    if(i>1){ 
      projF[,i]<-matF%*%projF[,i-1]}
    
    if(i==1){
      projM[,i] <- matM%*%popM
      projM[1,i]<-(matM%*%popF)[1]}
    if(i>1){
      projM[,i]<-matM%*%projM[,i-1]
      projM[1,i]<-(matM%*%projF[,i-1])[1]}
    
  }
  output<-list(projF=projF, projM=projM)
  return(output)
}

