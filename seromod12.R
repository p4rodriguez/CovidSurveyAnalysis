
# ========================================
#   Build and model data for seroprevalance
#  (ie the cumulative number of infections)
# ========================================
#     Start here, the first script in pipeline
#  this calls GetData to build matrix of data for
#  each state from covid tracking and cmu survey
# data  (at state aggragate level).
# This takes set of data matrices, uses the 
#  'getX2sv4sero..R' script to get data matrix
#   to save ('X2sv') prior to modeling, it extracts
#  data across states at 1 time point
#  and uses random forest (and tries an OLS model)
#  to estimate cumulative I (infctns) for all time points.
#  After this script, use the Imod..R and Idata4mod..R scripts
#  to estmiate I  pool at all time point (infectious pool)
#  ------------------------------------------------------
# ========================================
source("helperfunctions10.R")

#-------get seroprevalence by state (from ref Stanford art.)
seroxcodes=get_seroprev()
st2dosetUp  =seroxcodes[,'StateU'] #uppercase state abbrev
st2doset    =seroxcodes[,1] #list of states to run 
stmax= length(st2doset)

pop         =read.csv('State Populations.csv',stringsAsFactors = FALSE)
seroxcdxst  =inner_join(seroxcodes,pop,by=c("StateU"="State"))
Pset        =seroxcdxst[,dim(seroxcdxst)[2]]

# ======== gather up survey and covid tracking data
if (1){
ageval2do = "overall"
begdate    ="2020-04-01"  #only take data after this
D_startdate='2020-04-15'  #line up data wrt this date
fstdatenum =as.Date(D_startdate)-as.Date(begdate)
source("./GetD10.R")
DsetOv    = Dset
NsetOv    = Nset
n2do=min(as.numeric(NsetOv[,5]))

#======= get excess death counts from CDC csv file
exdeaths=get_exdeaths()
}


# ================
# these col inds for get func, the cols to get from final joins
# in GetD...R script
# (see nms for the names)
c2sv =c(1,2:4,29,33)  #12,18:20,22:27,29,33)  
#c2sv =c(1,2:3,29,33)  #12,18:20,22:27,29,33)  
outi =c(22,27,31,40,47) #27,31,47,22) #c(1,12,27,31,47,22)  #states not 2use,ith row
gdi  =seq(1:51)
gdi  =gdi[-outi]

Y           =seroxcodes[gdi,"sp"]

#Now last part, run models --------------------
#  eni is an 'end' index, eg take data upto this date
enidate  ="2020-06-15"  
      #1st try i did 6-15 but SeroPrev;ABodies take 1-3wks to apper
enidate  ="2020-06-22"  #so maybe 2 weeks back from 2ndwkjuly?

eni2use  = as.Date(enidate)-as.Date(D_startdate)
doderivs=1  
doStmeta=1
   #eg for getting lag info, use for Infctn case models
source('getX2sv4sero_15.R')

#cols to use for a model
c2use=c(#"p_avoidS","p_avoid0","p_avoid1","p_avoid2",
  #"pcwS",
#    "pcw1","pcw2", #not really help much...mmm
  #"ppvS",
  #"ppv1","ppv2",
  #"pcw2",
#  "stringindex0",
        #"pmmedconS",
        #"pas2_1smf0",
        #"pas2_1smf1",
        #"pas2_1smf2",
        #"lidf1",
        #"lidf2", 
        "a19to34","a65p",
        #"ppv","ppv2",
        "pdens",
        "cap",
        "reg",
        "exdeaths",
#        "hosps",
        "deaths")
#age,pdens,cap,regfactor,exdeaths,deaths
# fit is .89 cor to Y; still few states with some neg or funny
# slopes far from 6-30.. but ok

X2mod =X2sv[,c2use]
X2mod        = as.data.frame(X2mod) #%>% mutate(tm=1:47)

#-------------------- clean up X
#make reg5 indicator
X2mod[,"reg"]= as.factor(X2mod[,"reg"])
X2mod        =X2mod[-outi,]

#r5i     =which(X2mod[,"reg"]==5)
#r5ind   =matrix(0,nrow(X2mod),1)
#r5ind[r5i]=1  
#X2mod[,"reg"]=r5ind
#X2mod[,'hosps']=fillhosps(X2mod)

X2mod_eni=X2mod
lmr      =lm(Y~.,data=X2mod_eni)
print(summary(lmr))
#print(cor(subset(X2mod_eni,select=-c(reg)),Y))
library(randomForest)
rf1=randomForest(X2mod_eni,Y,ntree=10000,importance = TRUE)

print("****** corr lmr, rf1 fitted *******")
print(cor(lmr$fitted.values,Y))
print(cor(predict(rf1,X2mod_eni),Y))

#do a LOO test
if (1){
  pred_yout=matrix(0,nrow(X2mod_eni),1)
  lmpred_yout=matrix(0,nrow(X2mod_eni),1)
  for (i in 1:nrow(X2mod_eni)){
  testx=X2mod_eni[i,]
  testy=Y[i]
  trnx =X2mod_eni[-i,]
  trny =Y[-i]
  rftst         =randomForest(trnx,trny,ntree=10000)
  pred_yout[i]  =predict(rftst,testx)
  lmtst         =lm(trny~.,data=trnx)
  lmpred_yout[i]=predict.lm(lmtst,testx)
}  
  
print("****** corr LOO rf, lm, mn abs diff *******")
print(cor(pred_yout,Y))
print(cor(lmpred_yout,Y))
print(mean(abs(pred_yout-Y)))
print(mean(abs(lmpred_yout-Y)))
#varImpPlot(rf1)
} #end do LOO test 
# --------------------------------------------
#Now get predictions of infected amount
# based on deaths and other columns for all 'end-indices'
#
eniset  =seq(21,n2do,1); #2 is Dstart+5 ~4/20,
                        #data should be avail for stuff
#corxeni =matrix(0,length(eniset),length(c2use)-1) #reg is a factor
li=0
lmrps1   =matrix(0,stmax,length(eniset)) #length(Y),length(eniset))
lmrps2   =matrix(0,stmax,length(eniset)) #length(Y),length(eniset))
lmrps2p   =matrix(0,stmax,length(eniset)) #length(Y),length(eniset))
lmrps3   =matrix(0,stmax,length(eniset)) #length(Y),length(eniset))
for (eni4d in eniset){
  #enidate  =as.Date("2020-04-01")+eni4d #90; #6/15 #6/30 for deaths
    #gathre eni for each state (end index) that is for 6/30
  #enixst   =get_ind4date(enidate)
  eni2use=eni4d
  source('getX2sv4sero_15.R')  #X2sv=getXfromDset(enidate)
  
  #X2sv   =getXfromDset(eni)
    X2mod  =X2sv[,c2use]
    X2mod  = as.data.frame(X2mod) 
    X2mod[,"reg"]= as.factor(X2mod[,"reg"])
#    X2mod[,"reg"]=r5ind
#    X2mod[,'hosps']=fillhosps(X2mod)

    lmr2   =predict(lmr,X2mod)
    lmrf   =predict(rf1,X2mod)
    li    = li+1
    lmrps3[,li]=(lmr2+lmrf)/2
    lmrps1[,li]=lmr2
    lmrps2[,li]=lmrf
    #   corxeni[li,]=cor(X2mod[c(1:4,6:7)],Y)
}
#seroxcdxst[51,"X2018.Population"]*.7/100
lmr_dates  =as.Date(D_startdate)+eniset

#now make a cumul sum adjusting for negative changes
#eg  make sure a negative change in cumul. infectino estimates
#(ie total that has been infected) is NOT cnted as decrease

for (sti in 1:stmax){
   Yxt    = lmrps2[sti,]  
   Yxt_l1 = lag(Yxt,1)
   d1     = Yxt-Yxt_l1
   d1[1]  = Yxt[1]
   d1[which(d1<0)]=0
   lmrps2p[sti,]=mysm2(mysm2(cumsum(d1),7),7) 
      #rand for. tends to make sharp decisions, this helps
}

#48th date is 6/22; so 
#test cor of
print("****** corr for xtime predict *******")
dateofsp = which(lmr_dates=="2020-06-22")
print(cor(lmrps2[-outi,dateofsp],Y))
print(cor(lmrps2p[-outi,dateofsp],Y))

if (0){  #when ready run this
for (sti in 1:stmax){
  
  stem=paste("cumIest_",st2doset[sti],"_",dateofsp,sep="")
  png(paste("./Images/IEstimates_figs/Fig_",stem,".png",sep=""))
  plot(lmrps2p[sti,],type='l',main=stem,ylab="total I %",xlabel="Time")
  dev.off()
}
}