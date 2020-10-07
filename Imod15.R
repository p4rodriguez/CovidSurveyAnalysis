# ----------------------------------
#now run this after seromod..R
# Get weekly infection amount using 
# lmrps  (lin mod result predictions for cumulative Infections)
#   but actually randfor model is used
# 1 row per state

#for eg
#max(cor(lmrps1[-outi,],Y))
#cor(Y,lmrps3[-outi,])) max is .94 vs .93 .90 vos lmr2 lmr1

#Note; eniset is set of days to track I, 
# first day is D_startdate+ eniset[1]
#lmr_dates  =as.Date(D_startdate)+eniset, done in seromod...R

#Save I estimates
IrlxT        =matrix(0,stmax,ncol(lmrps2))  #estimate of I pool at t
dIxT        =matrix(0,stmax,ncol(lmrps2)) #estimate of new I case at t
rec_numdays=7  # averge num days to recover, or still infcted
last_rec   =14 #but upto 14 
               #(this ignores deaths removal from active pool)

totIxt =matrix(0,stmax,ncol(lmrps2))
for (sti in 1:stmax) {
    I         = lmrps2p[sti,]  #this is cum total I adjst smthd
#    if (any(I<0))    { I=I-min(I)}  #shouldnt need this
    totIxt[sti,] = I 
    dI         = I-lag(I,1)    #the new infections is change of totalI
    dI[1]      = dI[2]         #1st  change is a jmp
    if (any(dI<0))    { dI=dI-min(dI)} 
    dIsm       =mysm2(dI,7)
    dIxT[sti,]  =dIsm
    
    #rolling sum, to simulate actively infectious (ie recovry time)
    dIrl       =myrlngsum(dI,rec_numdays,last_rec)
    if (any(dIrl<0,na.rm=TRUE)) { dIrl=dIrl-min(dIrl) } 
        #occasionly need this, esp. for states wout seroprev. data

    dIrlsm       =mysm2(dIrl,7)
    IrlxT[sti,]  =dIrlsm
}    
if (0){  #when ready run this
  for (sti in 1:stmax){
    
    stem=paste("PoolIcaseest_",st2doset[sti],"_",dateofsp,sep="")
    png(paste("./Images/PoolIEst_figs/Fig_",stem,".png",sep=""))
    plot(IrlxT[sti,],type='l',main=stem,ylab="I pool %",xlabel="Time")
    dev.off()
  }
}

if (0){  #when ready run this
  for (sti in 1:stmax){
    
    stem=paste("NewIcaseest_",st2doset[sti],"_",dateofsp,sep="")
    png(paste("./Images/NewIEst_figs/Fig_",stem,".png",sep=""))
    plot(dIxT[sti,],type='l',main=stem,ylab="new I %",xlabel="Time")
    dev.off()
  }
}

#Now get a new data matrix with surv and covid track data
# first take all columns, and engineered columns
# later trim to make model
# 

