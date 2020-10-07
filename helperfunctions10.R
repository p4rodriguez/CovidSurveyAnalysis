#functions
library(NlinTS)
library (timeSeries) # to extract time series
library(corrplot)
library(lmtest)  #grangertest

library(readr)
library(purrr)
library(dplyr)
library(data.table)
library(collections)

require(lubridate)
require(magrittr)

# ===========================
get_exdeaths=function(){
 ex_all=read.csv("Excess_Deaths_Associated_With_Covid-19.csv",stringsAsFactors = FALSE)
 gdi    =grep('exclud', ex_all[,'Outcome'])
 ex_allx=ex_all[gdi,]
 
 tmpd   =as.Date("1900-01-01")+(ex_allx[,1]-2)
 tmpdnum=tmpd-as.Date(begdate)
 ex_all2=data.frame("date"=tmpd,
                   "datenum"=tmpdnum,
                   "State"=as.character(ex_allx[,"State"]),
                   "exlow"=as.numeric(ex_allx[,"ExcessLowerEstimate"]),
                   "exhigh"=as.numeric(ex_allx[,"ExcessHigherEstimate"]),
                   stringsAsFactors = FALSE)
 ex_all2[is.na(ex_all2)]=0
 gdi     =which(tmpdnum>=0)  #as.Date("2020-04-01"))
 exdeaths=ex_all2[gdi,]
 return(exdeaths)
}
# ========================================
get_seroprev=function(){
  sero     =read.csv("Sero_agexstate0630.csv",stringsAsFactors = FALSE)
  sero     =sero[order(sero[,1]),]
  st2doset    =sero[order(sero[,1]),1]
  sero[,"reg"]=as.factor(sero[,"reg"])

  stmax=length(st2doset)

  stcodes=read.csv("StateCodes.csv",stringsAsFactors = FALSE)
  colnames(stcodes)[1]="StateU"
  seroxcodes  =inner_join(sero,stcodes,by=c("State"="Low"))
  return(seroxcodes)
}

#-----------------
#each state has different 'starting points' so
# to align them figure out which index is this date 
get_ind4date=function(enidate){
  enixst   =matrix(0,stmax,1)
  for (si in 1:stmax){
    tmpi=which(DsetOv[[si]][,'date']==as.character(enidate))
    if (length(tmpi)==0) {
      tmpi=which(DsetOv[[si]][,'date']==as.character(enidate+1))
    }
    enixst[si]=tmpi
  }
  return(enixst)
}

get_mindate=function(){
  enims   = numeric( stmax )
  class(enims) = "Date"  
  for (si in 1:stmax){
    enims[si]=DsetOv[[si]][1,'date']
  }
  return(enims)
}


#fill in missing hospitalizations,  not enough data
fillhosps=function(X2mod){
  xh  =X2mod[,'hosps']
  yh  =subset(X2mod, select = - c(reg,hosps))
  bdi =which(is.na(xh))
  xh2 =xh[-bdi]
  yh2 =yh[-bdi,]
  misslm    = lm(xh2~.,data=yh2)
  predh     = as.matrix(cbind(matrix(1,nrow(X2mod),1), yh)) %*% misslm$coefficients
  xh[bdi]   = predh[bdi]
  if (any(xh<=0)){
    zi   = which(xh<=0)
    xh[zi]=mean(xh[-bdi])
  }
  return(xh)
}

#--------------
myviewlines=function(Xl,stm){
  outfile=paste(".\\Images\\St6_",stem,"_",stm,".png",sep="")
  png(outfile) 
  for (fi in 1:length(F_ppf)){
    ms=max(Xl[[fi]])
    x =fi+Xl[[fi]]*1/ms #scale to 1
    #    print(range(x))
    #    print(summary(x))
    if (fi==1){
      plot(x,type='l',ylim=c(0,length(F_ppf)+1))
    } else { points(x,type='l')}
  }
  dev.off()
}



#cut matrix to high/low values
mycutoff=function(X,th,side){
  gdiN=which(X<th)
  gdiP=which(X>th)
  gdiN2=which(X< (-th))
  gdiP2=which(X> (-th))
  X=0*X
  if (side=='3')
  { X[gdiN2]=-1 
  X[gdiP]=1}
  if (side=='2')
  { X[gdiN]=-1 
  X[gdiP]=1}
  if (side=='L')
  {X[gdiN]=-1}
  if (side=='R')
  {X[gdiP]=1}
  return(X)
}
#remove outliers by quantiles, replace with median
#eg for bad data entry, median of >0 values
myremoiqr=function(x){
  spd=IQR(x)
  c  =3
  mn =median(x[which(x>0)])
  if (spd>0) { st=spd
  } else { st=sd(x)}
  out1p=which(x>(mn+c*st))
  out1n=which(x<(mn-c*st))
  x[out1p]=mn
  x[out1n]=mn
  return(x)
}
#for real values, replace with mean maybe viz median
myremoiqrmn=function(x){
  spd=IQR(x)
  c=3
  mn=mean(x)
  st=spd
  out1p=which(x>(mn+c*st))
  out1n=which(x<(mn-c*st))
  x[out1p]=mn
  x[out1n]=mn
  return(x)
}
#for normal dist. outliers..not used
myremo3st=function(x){
  c=3
  mn=mean(x)
  st=sd(x)
  out1p=which(x>(mn+c*st))
  out1n=which(x<(mn-c*st))
  x[out1p]=mn
  x[out1n]=mn
  return(x)
}
mysm=function(x,t){
  ma={}
  if (t>0)  
  { for (i in seq(t,length(x)))
  {  bei=i-t+1;  eni=i;
  #     print(paste('i',i,' b:',bei,' en',eni,' v',x[eni]))
  ma=append(ma,mean(x[bei:eni]))
  }
  }  else {ma=x;}
  return(ma)
}
#smooth but keep first indices, good for anextra smooth
mysm2=function(x,t){  
  ma={}
  if (t>0)  
  { for (i in seq(1,length(x)))
  {  bei=max(1,i-t+1);  eni=i;
  #     print(paste('i',i,' b:',bei,' en',eni,' v',x[eni]))
  ma=append(ma,mean(x[bei:eni]))
  }
  }  else {ma=x;}
  return(ma)
}



#smooth out 0  , only at 0 pts
mysm0=function(x,t=3){
  ma={}
  for (i in seq(1,length(x)))
  { if (x[i]<=0)
  {bei=max(1,i-t);eni=min(length(x),i+t);
  ma =append(ma,mean(x[c(bei:(i-1),(i+1):eni)]))
  }  else {ma=append(ma,x[i]);}
  }
  return(ma)
}

#-------------------------------
# a rolling summation, looking backward
myrlngsum=function(x,t1=7,t2=14){  
  ma              ={}
  hf_wts          =exp(-1*c(1:t2)/t2)
  w2use           =matrix(1,t2,1)
  w2use[(t1+1):t2]=hf_wts[1:(t2-t1)]
  t=t2 #end t
  for (i in seq(1,length(x)))
  {  bei=max(1,i-t+1);  eni=i;
  #     print(paste('i',i,' b:',bei,' en',eni,' v',x[eni]))
#  ma=append(ma,sum(x[bei:eni]))
  ma=append(ma,sum(x[eni:bei]*w2use[1:(eni-bei+1)] ))
  }
  return(ma)
}

