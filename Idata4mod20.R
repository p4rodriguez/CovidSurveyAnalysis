#run this after Imod...R  (which gets the infectns I to model)
#
# here, get a new data matrix with surv and covid track data
# first take all columns, and engineered columns
# later trim to make model*** FOR each STATE  ***
# 
#  for the following:
#IrlxT  estimate of I pool at t, this can give TPrate est.
#dIxT estimate of new I

doderivs   =1
eniset2    =eniset  #   
Imod_d1ind =eniset2[1]-eniset[1]+1  
            #the Imod data 1st entry is an earlier date

dI_dates   =as.Date(D_startdate)+eniset2

doStmeta   =0;
c2sv       =c(2:20,21,22:27,28,29,33,34)  #set of column indices to save from 
                              #the Dset list of matrices (survey&covid data)

#====================
#   1  build data matrices for each state
#     it will have more columns than we need
#     so that we can explore what we need/want
#====================
#Get a first data matrix
eni2use = eniset2[1]
source('getX2sv4sero_15.R') #for each state 1 time point
X4I = X2sv;   

#Initialize list of state data matrices 
X4Iset={}  
tmpX  = matrix(0,length(eniset2),ncol(X4I))  #time is going down
for (sti in 1:stmax){
  tmpX[1,]     =X4I[sti,]
  X4Iset[[sti]]=tmpX;
}

#now do rest of states, adding 1 row time point each loop
rowi  = 1 #row index to add (eg next time point)
for (rowi in 2:length(eniset2)){
  eni2use =eniset2[rowi]
  source('getX2sv4sero_15.R') #for each state 1 time point
  X4I     = X2sv;  #this is 1 col 4 final set of data matrices 
  for (sti in 1:stmax){
     X4Iset[[sti]][rowi,]=X4I[sti,]
     }
}

#====================
#   2 now get X data matrix
#====================

#now check out columns
#cols to use for a model 

pcw_c2use=c("datenum0","pcw0",
            "pcw1","pcw2") 

basic_c2use=c("datenum0","pcw0",      
              "pcw1","pcw2",
              "pas2_1smf0",
              "pas2_1smf1","pas2_1smf2"
#              "lidf0","lidf1",   "lidf2"
              )

#        "posinc0","ppv0") 

wder_c2use=c("datenum0","pcw0", "pcw1","pcw2",#"lidf0",
        "p_avoid0","p_avoid1","p_avoid2",
        "pas2_1smf0",
        "pas2_1smf1","pas2_1smf2",
        "posinc0","posinc1","posinc2",
        "ppv0","ppv1","ppv2","totTrinc0")
#c2use=c(2,5:(ncol(X4I)-4))
c2use=c(wder_c2use)
#c2use=pcw_c2use

XtoYcorres=matrix(0,stmax,length(c2use)) #save corrltn rsults
lm4I_set  ={}
dIYset    = dIxT[,Imod_d1ind:ncol(dIxT)]
          #  matrix(0,stmax,ncol(dIxT)-Imod_d1ind+1)
IrlYset   = IrlxT[,Imod_d1ind:ncol(dIxT)]
X4Imodset ={} #the model data matrix is subset of full data matrix

#now for each state's data matrix X
#  save some correltn and pvalues wrt daily new-case estimate (dI)
for (sti in 1:stmax){
  dIY  =dIYset[sti,] #dxT[sti,Imod_d1ind:ncol(dIxT)] #get Y values to model
  dIY = c(dIY[2:length(dIY)],dIY[length(dIY)])
  
  IrlY  = IrlYset[sti,];
  tmpX               =X4Iset[[sti]]
  colnames(tmpX)=xnms;
  X4Imod             =tmpX[,c2use]  #take just these colmns
  #X4Imod  = as.data.frame(X4Imod)
  X4Imodset[[sti]]=X4Imod  #save for later
  XtoYcorres[sti,]=cor(X4Imod,IrlY) #dIY)
#  lm4I_set[[sti]]=lm(dIY~.,data=as.data.frame(X4Imod)) #rum a lin model too
  lm4I_set[[sti]]=lm(IrlY~.,data=as.data.frame(X4Imod)) #rum a lin model too
}
#now XtoYcorres is corr of each data col to dIY (the newcases estimate)
rownames(XtoYcorres)=st2doset
colnames(XtoYcorres) =c(c2use)
corrplot(XtoYcorres,tl.cex=.4)

#------------- to get summary if lin reg coefficients
#  a clustering analysis showed 3 or 4 patter of 
#  coefficents, for example, of how survey data predicts
#  this I new case estimate, 
#Now get snapshot of each states model to the newcases (dIY)
if (0){
 lm4i_pvs =matrix(0,stmax,length(c2use)+1)
 lm4i_coef=matrix(0,stmax,length(c2use)+1)
 lm4i_fp  =matrix(0,stmax,length(c2use)+1)
 for (sti in 1:stmax){
  lmsum           =summary(lm4I_set[[sti]])
  lm4i_pvs[sti,]  =lmsum$coefficients[,4] #save pvalues
  lm4i_coef[sti,] =lmsum$coefficients[,1] #save coeff itself
  lm4i_fp[sti,1]=pf(lmsum$fstatistic[1],lmsum$fstatistic[2],lmsum$fstatistic[3],lower.tail=FALSE)
  lm4i_fp[sti,2]=mean(abs(lm4I_set[[sti]]$fitted.values-dIYset[sti,]))
  lm4i_fp[sti,3]=max(abs( (lm4I_set[[sti]]$fitted.values-dIYset[sti,])/
                            1))
 }

#rownames(lm4i_pvs)=st2doset
#colnames(lm4i_pvs)=c("Inter",c2use)
#corrplot(1-lm4i_pvs,tl.cex=.4)
mean(lm4i_fp[,2])
mean(lm4i_fp[,3])
ncT=ncol(dIxT)
plot(lm4I_set[[10]]$fitted.values[2:(ncT)],IrlYset[10,1:(ncT-1)])
}

#----------
# For some phase plots
# aka, 2 cols at time
# ------------------------------------
source('segment10.R')
#nc=ncol(X4Imod)


#1
# to get correlation of pcw,ppv at 1 pt in time
# across states (cross sectional snap shot)
# at different days (di) date indices
Yc=c("ppv0")
Xc=c("pcw0")
for (di in c(34,48,62,76,90,104)){
  xy=matrix(0,stmax,2)
  for (si in seq(1,stmax)){
    tmpX  =X4Imodset[[si]]
#    dIY   =dIYset[si,] 
#    IrlY  = IrlYset[sti,];
    
    X     =tmpX[,Xc]
    Y     =tmpX[,Yc]
    xy[si,1:2]=c(X[di],Y[di])  #get time slice, 
  }
  md2=substr(as.character(lmr_dates[di]),6,10)
  plot(xy[,1],xy[,2],cex=.25,ylim=c(0,.30),xlim=c(0,1.5),
       ylab="PosRate",xlab="PCLI Wtd",
       main=paste("cor=",round(cor(xy[,1],xy[,2]),2)," on ",md2))
  text(xy[,1],xy[,2],st2doset,cex=.75,col='red')
} #end di


#2
#for diff states get phase plots
Yc  =c("ppv0") #or posinc0"
Xc  =c("pcw0")
Yc  =c("posinc0") #or posinc0"
Xc  =c("totTrinc0")
bei =25;eni=length(eniset2); #
d2i  =which(lmr_dates=="2020-07-01")
d3i  =which(lmr_dates=="2020-08-01")
#d1i  =which(lmr_dates=="2020-06-01")+bei-1
dts2plt=c(bei,eni) 
for (si in 1:10) { #1:stmax) { # seq(1,stmax)){
  dts2plt=c(bei,eni) 
  tmpX  =X4Imodset[[si]]
  dIY   =dIYset[si,]
  IrlY  = IrlYset[si,bei:eni];

  #combine 2 peak/alley finding functions
  pks1    =findPeaks(IrlY,0)
  vls1    =findValleys(IrlY,0)
  pks2    =find_peaks(IrlY,7)+1
  vls2    =find_peaks(-1*IrlY,7)+1

  pks = intersect(pks1,pks2) #pks 
  vls = intersect(vls1,vls2)
  
  if (length(pks)==0) pks=which.max(IrlY)
  if (length(vls)==0) vls=which.min(IrlY)
  
  stem=paste("./Images/Gestimate_tprs/",st2doset[si],"_",si,"_IrlY2.png",sep="")
  png(stem)  
  plot(IrlY,type='l',main="Inftd. pool estimate as %",
       ylim=c(0,max(IrlY)+.01));
  text(pks,pks*0+.01,"+",col='red',cex=1.5)
  text(vls,vls*0+.01,"V",col='red',cex=.75)
  for (i in pks){
    md =substr(as.character(lmr_dates[bei-1+i]),7,10)
    text(i+1,0.06,md,col='blue',cex=.75)  
  }
  for (i in vls){
    md =substr(as.character(lmr_dates[bei-1+i]),7,10)
    text(i+1,.06,md,col='blue',cex=.75)  
  }
  dev.off()
  
  if (Yc=="ppv0")
    { Y     =tmpX[bei:eni,'ppv0']
    } else {
     Y     =tmpX[bei:eni,"posinc0"]/((irl[bei:eni]/100)*Pset[si])
              #tprate, TruePos / estimate of Infct pool 
    }
  X     =tmpX[bei:eni,Xc]
  
  stem=paste("./Images/Gestimate_tprs/",st2doset[si],"_nTest2.png",sep="")
  png(stem)  
  nTest =tmpX[bei:eni,"totTrinc0"]/1000
  plot(nTest,type='l',main="Num Test Results (1k) ");
  dev.off()
    
  stem=paste("./Images/Gestimate_tprs/",st2doset[si],"_phase_2.png",sep="")
  png(stem)  
  plot(X,Y,type='l',ylab=Yc,xlab=Xc,ylim=c(0,max(Y)+.01),
       main=paste(st2doset[si],round(cor(Y,X),2)))
  dts2plt=c(dts2plt,pks+bei-1,vls+bei-1)
  for (di in 1:length(dts2plt)) {
    i4dt =dts2plt[di]
    md   =substr(as.character(lmr_dates[i4dt]),7,10)
    i4xy =i4dt-bei+1
    text(X[i4xy]+0,Y[i4xy]+.005,md,col='blue',cex=.75)  
    }
  #md=substr(as.character(lmr_dates[bei]),7,10)
  #text(X[1],Y[1],md,col='blue',cex=.75)  
  for (i in pks){
    text(X[i],Y[i],"+",col='red',cex=1.5)  
  }
  for (i in vls){
    text(X[i],Y[i],"V",col='red',cex=1)  }
  dev.off()
}
#if you are only testing pcw, corr shold be high
#
