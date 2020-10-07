#Get a data matrix, X2sv (X matrix to save with all potential
# columns)  of state information to model
# sero prevalence on 6/15/20, then apply it across time
# (see seromod..R script, 

#this is used to get a data matrix to model the estimate of
#  Infectious amount across time
# Use survey and covid tracking data, doderivs will 
#  also get changes in values as new engineeered features
#  but requires 20 time steps to look back


  X    =matrix(0,stmax,n2do)
  nms  =colnames(DsetOv[[1]])
  c2i  =0 #index which of colset 
  
#  X2sv =matrix(0,length(gdi),length(c2sv)*4) #build a data matrix to 
  X2sv =matrix(0,stmax,length(c2sv)*5+15) #build a data matrix to 
  xnms={}  #store col names
  
  for (ci in c2sv){
    #print(ci)
    for (si in seq(1,stmax)){
      #eni =enixst[si] 
      #eni  =enid-15   
      #line up to index of interest
      begi=1  #+(eni-enimin) #start here to align all states
      endi=n2do      #

      x1 =(DsetOv[[si]][,ci])   #/100)  *Pset[si]
      if (ci==which(nms=="deathsinc") || 
          ci==which(nms=="nsize")) {
              #ci==which(nms=="posinc"))  {   #normalize to percnts
             x1=100*x1/Pset[si]  }
      X[si,begi:endi]=x1[begi:endi]
      
    } #end stats loop
    #now for this col, each row is 1 state's valuesXtime
    if (ci==which(nms=="deathsinc")) { 
        #|| 
        #ci==which(nms=="nsize") 
        #||
        #ci==which(nms=="posinc")) {
        #ci==which(nms== "pcw"))  {   
          #take a row sum upto enimin index
          rs90 =rowSums(X[,1:eni2use]) #all state are aligned 
          c2i       =c2i+1
          X2sv[,c2i]=rs90
          xnms      =append(xnms,nms[ci])
     } else { #take last 2 weeks avgs and trend
      #bp=c(12,8,4) first tests
       rsmS=rowSums(X[,1:eni2use])
       c2i =c2i+1
       X2sv[,c2i]=rsmS
       xnms=append(xnms,paste(nms[ci],"S",sep=""))
   if (doderivs) {
      rsm3=rowMeans(X[,(eni2use-20):(eni2use-14)])
      rsm2=rowMeans(X[,(eni2use-13):(eni2use-7)])
      rsm1=rowMeans(X[,(eni2use-6):eni2use])
      c2i =c2i+1
      X2sv[,c2i]=rsm1
      xnms=append(xnms,paste(nms[ci],"0",sep="")) #0th deriv
      c2i =c2i+1
      X2sv[,c2i]=rsm1-rsm2
      xnms=append(xnms,paste(nms[ci],"1",sep=""))#1st deriv (diff)
      c2i =c2i+1
      X2sv[,c2i]=(rsm1-rsm2)-(rsm2-rsm3) #2nd der
      xnms=append(xnms,paste(nms[ci],"2",sep=""))
       } #end do derivs
      } #else take avgs
  }  #column loop

  #print('adding fields..')   now add other fields
 if (doStmeta) {
  c2i       =c2i+1
  X2sv[,c2i]=seroxcodes[,3]
  xnms      =append(xnms,"a19to34")
  c2i       =c2i+1
  X2sv[,c2i]=seroxcodes[,4]
  xnms      =append(xnms,"a65p")
  c2i       =c2i+1
  X2sv[,c2i]=seroxcodes[,5]
  xnms      =append(xnms,"pdens")
  c2i       =c2i+1
  X2sv[,c2i]=seroxcodes[,6]
  xnms      =append(xnms,"reg")
  c2i       =c2i+1
  X2sv[,c2i]=seroxcodes[,7]
  xnms      =append(xnms,"cap")

  #==========get these fields too,but just 1 data pt
  # per state, no need to take cumsum
  rs90=matrix(0,stmax,1)  
  #get each states death cnt on eni-th day
  for (si in seq(1,stmax)){
    x1      =(100*DsetOv[[si]][,"hosps"])/Pset[si]   #/100)  *Pset[si]
    rs90[si]=x1[eni2use]
  }
  c2i       =c2i+1
  X2sv[,c2i]=rs90[]
  xnms      =append(xnms,"hosps")
  
  rs90=matrix(0,stmax,1)  
  #get each states death cnt on eni-th day
  for (si in seq(1,stmax)){
    x1      =(100*DsetOv[[si]][,"deaths"])/Pset[si]   #/100)  *Pset[si]
    rs90[si]=x1[eni2use]
  }
  c2i       =c2i+1
  X2sv[,c2i]=rs90[]
  xnms      =append(xnms,"deaths")
  
  rs90=rs90*0
  for (si in seq(1,stmax)){
    exst  = exdeaths %>% filter(State==st2dosetUp[si]) 
    exst  = exst[order(exst[,'datenum']),]
    exeni = which(exst[,'datenum']<=eni2use & 
                  exst[,'datenum']>fstdatenum)
    tmp = matrix(0,eni2use,1)
    tmp[exst[exeni,'datenum']]=exst[exeni,"exhigh"] 
      #fill in wks values, exlow or exhigh?
    tmps    =smooth.spline(mysm2(tmp,14))
    tmpsy = tmps$y  #mysm2(tmps$y,14)
#    plot(mysm(tmps$y,1))
    rs90[si]=100*sum(tmpsy[1:eni2use ])/Pset[si]
    #sum(exst[exeni,"exhigh"])/Pset[si]
  }
  c2i       =c2i+1
  X2sv[,c2i]=rs90[]
  xnms      =append(xnms,"exdeaths")
 } #end of do state meta data columns  
  #c2mod =c(4,10,11,12,26,29,33)  #set of cols to get
  X2sv          =X2sv[,1:c2i]
  colnames(X2sv)=xnms
#  return(X2sv)
#}