#Load and join data sets

#dofig=0  #to save to file
#dofigsum=1
#th2do  =.6    #is correlatn 'high' as this
#gth2do =.95   #is the 1-GCI pvalue sig at this level 
#ord2do= 1 #for checking lag of x vs y,and lag1 of y
#1 

#Alabama, Arkansas, Delaware, District of Columbia, Florida, Georgia,
#Kentucky, Louisiana, Maryland, Mississippi, North Carolina, Oklahoma, 
#South Carolina, Tennessee, Texas, Virginia, and West Virginia
#DC no 55+?

#st2doset = unique(D[,"state_code"])[c(1:10,11:50,57)]
#st2doset=  c("al","ar","de","fl") 
#             #,"ga","ky","la","md","ms","nc","ok","sc","tn",
#"tx","va","wv")
#st2doset=c("fl")

#-------------------------------
D0     =read.csv('overall-state.csv')
#D0     =read.csv('overall-state-smoothed.csv')
#begdate="2020-04-16"  #for itally "2020-02-23"
tmpd   =as.Date(D0[,'date'])-as.Date(begdate)
D0[,'datenum']=tmpd
badi          = which(D0[,'datenum']<=0)
if (length(badi)>0)
{  D             = D0[-badi,]
 } else { D=D0 }

S0  =read.csv('https://covidtracking.com/api/v1/states/daily.csv',stringsAsFactors=FALSE)
tmpd=as.Date(as.character(S0[,'date']),'%Y%m%d')-as.Date(begdate)
S0[,'datenum']=tmpd
badi          = which(S0[,'datenum']<=0)
S             = S0[-badi,]

R0            = read.csv('../oxford-govtresp/USA-covid-policy-master/data/OxCGRT_US_latest.csv')
tmpd          = as.Date(as.character(R0[,'Date']),'%Y%m%d')-as.Date(begdate)
R0[,'datenum']= tmpd
badi          = which(R0[,'datenum']<=-50)
R             = R0[-badi,]
#-------------------------------


#-------------------------------
sml   =7
sml2  =7
stmax =length(st2doset)
nsi   =0
Nset  =matrix(0,100+stmax+1,8)
Dset  ={}
#-------------------------------

#-------------------------------
#-------------------------------

for (ageval in ageval2do){
  #c("55+")) {   #"overall")) {
  for (genval in c("overall")){
  #     c("female","male","overall")) {

    DFsum =0
    CRsum =0
    
for (stcnt in seq(1:stmax)){
  
  st2do=st2doset[stcnt]
  print(paste('st:',stcnt,st2do))
  
  #1 String. Index -----------------------
  Rxst=filter(R,RegionCode==paste("US_",toupper(st2do),sep=""))
  NR = nrow(Rxst)
  R2      =data.frame(
    'datenum'   = Rxst[sml:NR,'datenum'],
    'stringindex' = Rxst[sml:NR,'StringencyIndexForDisplay'])
  #-------------------------------
  
  #2 CovTrack State counts  
  #print('starting S data..')
  Sxst=filter(S,state==toupper(st2do)) 
  Sxst=Sxst[order(Sxst[,'date']),]
  N   =nrow(Sxst)
  pos        =Sxst[,'positive']
  posinc     =Sxst[,'positiveIncrease']
  deathsinc   =Sxst[,'deathIncrease']
  #deathsincF  =lag(deathsinc,-21) #mv it ahead 21 days
  hospsinc    =Sxst[,'hospitalizedIncrease']
  deaths      =Sxst[,'death']
  hosps       =Sxst[,'hospitalized']

  pos[is.na(pos)]      =0
  pos[which(pos<0)]    =0
  posinc[is.na(posinc)]=0
  rec                  =Sxst[,'recovered']
  rec[is.na(rec)]      =0
  if (0) {pos_mr=pos-rec}
  pos_sm    = mysm(pos,sml)
  #posmr_sm  = mysm(pos,sml)
  #_diff  = diff(pos_sm)
  totTresinc = Sxst[,'totalTestResultsIncrease']
  #totTresinc[which(totTresinc<0)]=1
  
  
 # print('..  totTres part')
  #totTres_sm       = mysm(Sxst[,'totalTestResults'],sml)
  
  #maybe for some states do this and use it?
  totTres      = Sxst[,'totalTestResults']
  if(st2do=="tx")
  { totTresincd = Sxst[,'totalTestResultsIncrease']
    badi=which(totTresincd<0)
    totTresincd[badi]=mean(totTresincd[c(badi-1,badi+1)])
    } else {
     Naft = length(totTres)
     totTresincd         = totTres
     totTresincd[2:Naft] = totTres[2:Naft]-totTres[1:(Naft-1)]
      totTresincd[1]      = totTresincd[2]
    }
  #posinc =pos[2:N]-pos[1:(N-1)]
  #testinc=tests[2:N]-tests[1:(N-1)]
  
  #smooth the big negative increase of GA with
  #  linear interpolation (b/c it should be steady changes
  #   and the next smoothing will make it ~differentlble)
  if (any(totTresincd<0,na.rm=TRUE)) #st2do=='ga' )
       { negi=which(totTresincd<0)
        negimid=negi[ceiling(length(negi)/2)]
        bbi=max(1,negi[1]-10)
        bei=min(length(totTresincd),negi[1]+10)
        dfi=totTresincd[bbi]-totTresincd[bei]
        for (bbii in bbi:bei){
           adj=(bbii-bbi)*dfi/(bei-bbi) #slope X index
           totTresincd[bbii]=totTresincd[bbi]+adj
        }}
  
  i0 = any(totTresincd<=0)
  if (i0) totTresincd=mysm0(totTresinc)
  i0 = any(posinc<=0)
  if (i0) posinc    =mysm0(posinc)

  posinc_sm       = mysm2(mysm(posinc,sml),sml) #doulbe smooth
  totTresinc_sm   = mysm(totTresincd,sml)
#  ppv         =mysm(posinc_sm/totTresinc_sm,sml)
  ppv         =mysm2(posinc_sm/totTresinc_sm,sml)
  if (st2do=="tx") ppv = mysm2(ppv,15) #
#  print('..ppv part') 
  ppv[is.na(ppv)]=0    #just in case
  ppv[which(ppv>1)]=1  

  #get deathsinc forward by 21 days
  deathsinc_sm=mysm(deathsinc,sml)
  dincfwd     =deathsinc_sm[21:length(deathsinc_sm)]
  dincfwd[length(dincfwd):length(deathsinc_sm)]=0
   
  S2      =data.frame(
   'datenum'     =Sxst[sml:N,'datenum'],
   'posinc'   =posinc_sm,
   'totTrinc' =totTresinc_sm,
   'ppv'         =ppv,
   'pos'         =pos_sm,
   'totTr'       =mysm(totTres,sml),
   'deaths'      =mysm(deaths,sml),
   'deathsinc'   =mysm(deathsinc,sml),
   'deathsincF'  =dincfwd,
   'hosps'       =mysm(hosps,sml),
   'hospinc'     =mysm(hospsinc,sml) )
 
 badi=which(S2[,"ppv"]==Inf)
 S2[badi,"ppv"]=mean(S2[-badi,"ppv"])
 badi=which(S2[,"ppv"]<0)
 S2[badi,"ppv"]=mean(S2[-badi,"ppv"])
 
 #-------------------------------
 #SurveyData
 print('starting D data..')
 dxst=filter(D,state_code==st2do,
               gender==genval,     #"overall",
               age_bucket==ageval) #"overall") 

 N2  =nrow(dxst)

 tmp_soccont=myremoiqr(dxst[,"mean_outside_hh_contact_in_social_gatherings_ct_weighted"])
 tmp_cmclict=myremoiqr(dxst[,"mean_cmnty_cli_ct_weighted"])
 tmp_pcw    =myremoiqr(dxst[,"pct_cli_weighted"])
 #tmp_pcw    =tmp_pcw+1/dxst[,'n']
 print('..pcw part')
 zfl  =0
 difst=1 #use sml-sml2+1 if the smooth amounts are different
 pcwtmp  =mysm(tmp_pcw,sml2)
 if (any(pcwtmp==0)==TRUE){
   zfl     =length(which(pcwtmp==0))
   pcwtmp  =mysm0(pcwtmp,5) #smooth out 0s, still trailing but from1:N
 }
 phwtmp  =mysm(dxst[difst:N2,'pct_hh_cli_weighted'],sml2)
 if (any(phwtmp==0)==TRUE){
   phwtmp  =mysm0(phwtmp,5) #smooth out 0s, still trailing but from1:N
 }

 print('..dxst2 part')
 
 mxsml=max(sml,sml2)
 dxst2=data.frame(
  'date'     =dxst[sml2:N2,'date'],
  'datenum'  =as.numeric(dxst[sml2:N2,'datenum']),
  'nsize'    =mysm(dxst[difst:N2,'n'],sml2),
  'pcw'      =pcwtmp, #mysm(dxst[difst:N2,'pct_cli_weighted'],sml2),
  'piw'      =mysm(dxst[difst:N2,'pct_ili_weighted'],sml2),
 'phhc'      =phwtmp, #mysm(dxst[difst:N2,'pct_hh_cli_weighted'],sml2),
 'paaw'      =mysm(dxst[difst:N2,"pct_cli_anosmia_ageusia_weighted"],sml2),
 'pcmc_kn'   =mysm(dxst[difst:N2,'pct_cmnty_cli_weighted'],sml2),
 'pmmedcon'  =mysm(dxst[difst:N2,'pct_multiple_medical_conditions'],sml2),
 'psfever'   =mysm(dxst[difst:N2,'pct_self_fever_weighted'],sml2),
 'psnosym1'  =100-mysm(dxst[difst:N2,'pct_self_none_of_above_weighted'],sml2),
 'psmsym'    =mysm(dxst[difst:N2,'pct_self_multiple_symptoms_weighted'],sml2),
 'pcontpos'  =mysm(dxst[difst:N2,"pct_contact_covid_positive_weighted"],sml2), 
 #'mnhhout'=mysm(dxst[difst:N2,"mean_outside_hh_contact_in_social_gatherings_ct_weighted"],sml2),
 'ptpos'     =mysm(dxst[difst:N2,"pct_tested_and_positive_weighted"],sml2),
 'ptneg'     =mysm(dxst[difst:N2,"pct_tested_and_negative_weighted"],sml2),
 'pttry'     =mysm(dxst[difst:N2,"pct_could_not_get_tested_weighted"],sml2),
 'ptntry'    =mysm(dxst[difst:N2,"pct_did_not_try_to_get_tested_weighted"],sml2),
 'mncmcliw'  =mysm(tmp_cmclict[difst:N2],sml2),
 'mnouthh'   =mysm(tmp_soccont[difst:N2],sml2), #maybe notsmooth to find singular days
 'p_avoid'   =mysm(dxst[difst:N2,"pct_avoid_contact_all_or_most_time_weighted"],sml2),
 stringsAsFactors = FALSE
  )
#======= feat eng ============================= 
#mild/asymptom evidence?  pct_contact_covid_positive
#  maybe sub off ptpos+ptneg, to get asym and untested?
#
  ptsum =mysm(dxst[difst:N2,"pct_tested_no_result"],sml2)+
          dxst2[,'ptpos']+dxst2[,'ptneg']
  dxst2[,'ptsum']=ptsum
  
  as1     = dxst[difst:N2,'pct_self_multiple_symptoms_weighted']-
   dxst[difst:N2,'pct_self_fever_weighted']
 dxst2[,'pas1']=mysm(as1,sml2)
#another attempt, at least 1 sym - fever 
 as2     = ( (100-dxst[difst:N2,'pct_self_none_of_above_weighted'])
                   -dxst[difst:N2,'pct_self_fever_weighted'])
 dxst2[,'pas2_1smf']=mysm(as2,sml2)

 as3    = dxst[difst:N2,'pct_self_multiple_symptoms_weighted']-
          dxst[difst:N2,"pct_hh_cli_weighted"]                
 dxst2[,'pas3_mmh']=mysm(as3,sml2)
 as4    = dxst[difst:N2,'pct_self_multiple_symptoms_weighted']-
          dxst[difst:N2,"pct_cli_weighted"]                
 dxst2[,'pas4_mmc']=mysm(as4,sml2)
 if (any(as3<0) || any(as2<0)|| any(as1<0))
 { print("some as lt 0...")}

 lidf    = dxst[difst:N2,"pct_ili_weighted"]-
          dxst[difst:N2,"pct_cli_weighted"]                
 dxst2[,'lidf']=mysm(lidf,sml2)
 
 
 #================================ 
 dxst2           =dxst2[order(dxst2[,"datenum"]),]
 
 #-------------------------------joins
 DxS  = inner_join(dxst2,S2,by=c("datenum"="datenum"))
 DxSxR=inner_join(DxS,R2,by=c("datenum"="datenum"))
 
 #line up everything to 15 days after 0401
 daycut = as.Date(D_startdate)-as.Date(begdate)
 DxSxR=DxSxR[which(DxSxR[,'datenum']>daycut),]  
    #so all data starts same day
 Dset[[stcnt]]=DxSxR
 
 nc  =ncol(DxSxR)
 nr  =nrow(DxSxR)
 #-------------------------------
 nsi       =nsi+1
 Nset[nsi,]=c(as.character(st2do),stcnt,
              ageval,genval,nr,range(DxSxR[,'ppv']),zfl)
 

 }   #end stateloop

}} #end age,gen loops

Nset=Nset[1:stmax,]
