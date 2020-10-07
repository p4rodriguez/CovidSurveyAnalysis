#some functions to help get change-points (piecewise)
# and granger causality

#from internet
find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  bdi=which(pks<m)
  if (length(bdi)>0) pks=pks[-bdi]
  pks
}
#install.packages('quantmod')
require(quantmod)
require(segmented)

mygetseg=function(x,pl=0){
  t=seq(1:length(x))
  xd=data.frame(x=x,t=t)
  out.lm = lm(x~t,data=xd)
  seg.control(display = FALSE,fix.npsi=FALSE,it.max=999)
  o<-segmented(out.lm) #1 breakpoint for x
  if (pl) {plot(o);points(dIY);}
  return(o$psi[2])
}

#-------------------------------
#my grangertest but dont take all lags
mygrtest = function(x,y,l2do){ 
  xl    =lag(x,n=l2do)
  xldat =na.contiguous(cbind(x,xl))
  colnames(xldat)=c('x','xl')
  
  yl     =lag(y,n=l2do)
  yldat =na.contiguous(cbind(y,yl))
  colnames(yldat)=c('y','yl')
  
  gdat=as.data.frame(cbind(xldat,yldat)) #xlag data, etc..
  fmU=lm(y~yl+xl,data=gdat)
  fmR=lm(y~yl,data=gdat)
  wres=waldtest(fmU,fmR)
  return(wres[2,4])
}
#take 2 lags
mygrtest2 = function(x,y,l2do){ 
  xl    =lag(x,n=l2do)
  xl1    =lag(x,n=l2do+1)
  xldat =na.contiguous(cbind(x,xl,xl1))
  colnames(xldat)=c('x','xl','xl1')
  
  yl     =lag(y,n=l2do)
  yl1    =lag(y,n=l2do+1)
  yldat =na.contiguous(cbind(y,yl,yl1))
  colnames(yldat)=c('y','yl','yl1')
  
  gdat=as.data.frame(cbind(xldat,yldat)) #xlag data, etc..
  fmU=lm(y~yl+yl1+xl+xl1,data=gdat)
  fmR=lm(y~yl+yl1,data=gdat)
  wres=waldtest(fmU,fmR)
  return(wres[2,4])
}

mygrtestopt= function(x,y){
  gres=matrix(1,7,3)
  for (ii in seq(1,7,1)){
    gres[ii,]=c(ii,mygrtest12(x,y,ii)) 
  }
#print(gres)
    bestgi=which.min(gres[,2])
  return(c(gres[bestgi,]))
}

mygrtest1 = function(x,y,l2do){ 
  xl    =lag(x,n=l2do)
  xldat =cbind(x,xl)
  colnames(xldat)=c('x','xl')
  
  yl    =lag(y,n=1) #does x lagged-at-n have info for y~y1
  yldat =cbind(y,yl)
  colnames(yldat)=c('y','yl')
  
  gdat=as.data.frame(na.contiguous(
    cbind(xldat,yldat))) #xlag data, etc..
  fmU=lm(y~yl+xl,data=gdat)
  fmR=lm(y~yl,data=gdat)
  stmp=summary(fmU)
  wres=waldtest(fmU,fmR)
  if (stmp$r.squared<.999 && abs(fmU$coefficients)[3]>.001) {  #beware machine precsn
    return(wres[2,4])
  } else { return(1) }
}

mygrtest12 = function(x,y,l2do){ 
  xl    =lag(x,n=l2do)
  xldat =cbind(x,xl)
  colnames(xldat)=c('x','xl')
  
  yl    =lag(y,n=1) #does x lagged-at-n have info for y~y1
  yldat =cbind(y,yl)
  colnames(yldat)=c('y','yl')
  
  gdat=as.data.frame(na.contiguous(
    cbind(xldat,yldat))) #xlag data, etc..
  fmU=lm(y~yl+xl,data=gdat)
  fmR=lm(y~yl,data=gdat)
  stmp=summary(fmU)
  wres=waldtest(fmU,fmR)
  if (stmp$r.squared<.999 && abs(fmU$coefficients)[3]>.001) {  #beware machine precsn
    return(c(wres[2,4],fmU$coefficients[3]))
  } else { return(c(1,0)) }
}



