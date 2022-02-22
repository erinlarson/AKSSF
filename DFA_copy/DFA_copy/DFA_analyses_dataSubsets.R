#####
##  Dynamic Factor Analyses for AKSSF
##  Timothy Cline
#####

setwd('/home/aq_user/aq_virtual_machine/AKSSF/DFA')

#Load required libraries
library(dplyr)
library(TMB)
library(doParallel)
library(foreach)
library(MARSS)

#function for standardizing time series
Zscore <- function(x){return((x-mean(x,na.rm=T))/sd(x,na.rm=T))}

#source TMB DFA code
source('DFA_code.R')

#load metadata for region id
MetaDat <- read.csv('../../data/md_2022-02-08.csv',header=T,stringsAsFactors = F)

head(MetaDat)


#Sequence of days between June 1 and Aug 31
DOYseq<-seq(152,243)

#load stream temperature data
#create some new date/time columns
#add region column
#limit timeframe to chosen date range
Tdat <- read.csv('../../data/summer_data_wair2022-02-12.csv',header=T,stringsAsFactors=F) %>% 
  mutate(sampleDate = as.Date(sampleDate), DOY = format(sampleDate,'%j') %>% as.numeric(),Year=format(sampleDate,'%Y')%>% as.numeric()) %>%
  mutate(Region = MetaDat$Region[match(SiteID,MetaDat$Site)]) %>%
  filter(DOY %in% DOYseq)

head(Tdat)

if(sum(is.na(Tdat$airDT))>0) print('Error: Missing air temperature data')
#head(MetaDat)
#head(Tdat)

#build time series matrix of mean daily water temp
MeanDWT_TS <- Tdat %>% filter(Year==2019) %>% xtabs(meanDT ~ SiteID + DOY,data=.)
MeanDWT_TS[MeanDWT_TS==0] <- NA

#dump sites with <70% of data
MeanDWT_TS <- MeanDWT_TS[which(rowSums(!is.na(MeanDWT_TS)) >= round(0.8*length(DOYseq)) ),]

#build time series of matching air temp
MeanDAT_TS <- Tdat %>% filter(Year==2019, SiteID %in% row.names(MeanDWT_TS)) %>% xtabs(airDT ~ SiteID + DOY,data=.)

#Zscore response and air temp covariate
zMeanDWT_TS <- apply(MeanDWT_TS,1,FUN=Zscore) %>% t()
zMeanDAT_TS <- apply(MeanDAT_TS,1,FUN=Zscore) %>% t()

#set.seed(15)
#randomSet<-sample.int(nrow(zMeanDWT_TS),10)
#tIn<-zMeanDWT_TS[randomSet,]
#aIn<-zMeanDAT_TS[randomSet,]
#m1 <- MARSS(tIn,model=list(R='diagonal and unequal',m=1,D='diagonal and unequal'),covariates=aIn,form='dfa')

Sub1_Sites <- MetaDat %>% filter(subset1==1) %>% pull(Site) %>% unique()
Sub2_Sites <- MetaDat %>% filter(subset2==1) %>% pull(Site) %>% unique()
Sub3_Sites <- MetaDat %>% filter(subset3==1) %>% pull(Site) %>% unique()
  

#setup cluster for parallel processing
n.cores <- parallel::detectCores()
cl <- makeCluster(n.cores)
registerDoParallel(cl)

#Parallel loop to run all years and store in a list
#Each parallel instance is like its own R session
#We have to pass the proper packages to run the functions and load dynamic model codes 
Subset1_DFA <- foreach(y = 2011:2020,.packages=c('dplyr','TMB')) %dopar% {
  dyn.load(dynlib("dfa1tmb"))
  Sub1 <- Tdat %>% filter(SiteID %in% Sub1_Sites)
  MeanDWT_TS <- Sub1 %>% filter(Year==y) %>% xtabs(meanDT ~ SiteID + DOY,data=.)
  MeanDWT_TS[MeanDWT_TS==0] <- NA #xtabs puts zeros in for missing dates in the dataset, we set those exact zeros to NA
  
  #keep sites with >80% of data
  MeanDWT_TS <- MeanDWT_TS[which(rowSums(!is.na(MeanDWT_TS)) >= round(0.8*length(DOYseq)) ),]
  
  #build time series of matching air temp
  MeanDAT_TS <- Tdat %>% filter(Year==y, SiteID %in% row.names(MeanDWT_TS)) %>% xtabs(airDT ~ SiteID + DOY,data=.)
  
  #Zscore response and air temp covariate
  zMeanDWT_TS <- apply(MeanDWT_TS,1,FUN=Zscore) %>% t()
  zMeanDAT_TS <- apply(MeanDAT_TS,1,FUN=Zscore) %>% t()
  
  m1 <- runDFA(zMeanDWT_TS,NumStates=1,ErrStruc='DUE',EstCovar=T,indivCovar=T,Covars = zMeanDAT_TS)
  m1$ydat <- zMeanDWT_TS
  m1$xdat <- zMeanDAT_TS
  
  m2 <- runDFA(zMeanDWT_TS,NumStates=2,ErrStruc='DUE',EstCovar=T,indivCovar=T,Covars = zMeanDAT_TS)
  m2$ydat <- zMeanDWT_TS
  m2$xdat <- zMeanDAT_TS
  
  m3 <- runDFA(zMeanDWT_TS,NumStates=3,ErrStruc='DUE',EstCovar=T,indivCovar=T,Covars = zMeanDAT_TS)
  m3$ydat <- zMeanDWT_TS
  m3$xdat <- zMeanDAT_TS
  
  list(m1,m2,m3)
}

Subset2_DFA <- foreach(y = 2011:2020,.packages=c('dplyr','TMB')) %dopar% {
  dyn.load(dynlib("dfa1tmb"))
  Sub2 <- Tdat %>% filter(SiteID %in% Sub2_Sites)
  MeanDWT_TS <- Sub2 %>% filter(Year==y) %>% xtabs(meanDT ~ SiteID + DOY,data=.)
  MeanDWT_TS[MeanDWT_TS==0] <- NA #xtabs puts zeros in for missing dates in the dataset, we set those exact zeros to NA
  
  #keep sites with >80% of data
  MeanDWT_TS <- MeanDWT_TS[which(rowSums(!is.na(MeanDWT_TS)) >= round(0.8*length(DOYseq)) ),]
  
  #build time series of matching air temp
  MeanDAT_TS <- Tdat %>% filter(Year==y, SiteID %in% row.names(MeanDWT_TS)) %>% xtabs(airDT ~ SiteID + DOY,data=.)
  
  #Zscore response and air temp covariate
  zMeanDWT_TS <- apply(MeanDWT_TS,1,FUN=Zscore) %>% t()
  zMeanDAT_TS <- apply(MeanDAT_TS,1,FUN=Zscore) %>% t()
  
  m1 <- runDFA(zMeanDWT_TS,NumStates=1,ErrStruc='DUE',EstCovar=T,indivCovar=T,Covars = zMeanDAT_TS)
  m1$ydat <- zMeanDWT_TS
  m1$xdat <- zMeanDAT_TS
  
  m2 <- runDFA(zMeanDWT_TS,NumStates=2,ErrStruc='DUE',EstCovar=T,indivCovar=T,Covars = zMeanDAT_TS)
  m2$ydat <- zMeanDWT_TS
  m2$xdat <- zMeanDAT_TS
  
  m3 <- runDFA(zMeanDWT_TS,NumStates=3,ErrStruc='DUE',EstCovar=T,indivCovar=T,Covars = zMeanDAT_TS)
  m3$ydat <- zMeanDWT_TS
  m3$xdat <- zMeanDAT_TS
  
  list(m1,m2,m3)
}

Subset3_DFA <- foreach(y = 2011:2020,.packages=c('dplyr','TMB')) %dopar% {
  dyn.load(dynlib("dfa1tmb"))
  Sub3 <- Tdat %>% filter(SiteID %in% Sub3_Sites)
  MeanDWT_TS <- Sub3 %>% filter(Year==y) %>% xtabs(meanDT ~ SiteID + DOY,data=.)
  MeanDWT_TS[MeanDWT_TS==0] <- NA #xtabs puts zeros in for missing dates in the dataset, we set those exact zeros to NA
  
  #keep sites with >80% of data
  MeanDWT_TS <- MeanDWT_TS[which(rowSums(!is.na(MeanDWT_TS)) >= round(0.8*length(DOYseq)) ),]
  
  #build time series of matching air temp
  MeanDAT_TS <- Tdat %>% filter(Year==y, SiteID %in% row.names(MeanDWT_TS)) %>% xtabs(airDT ~ SiteID + DOY,data=.)
  
  #Zscore response and air temp covariate
  zMeanDWT_TS <- apply(MeanDWT_TS,1,FUN=Zscore) %>% t()
  zMeanDAT_TS <- apply(MeanDAT_TS,1,FUN=Zscore) %>% t()
  
  m1 <- runDFA(zMeanDWT_TS,NumStates=1,ErrStruc='DUE',EstCovar=T,indivCovar=T,Covars = zMeanDAT_TS)
  m1$ydat <- zMeanDWT_TS
  m1$xdat <- zMeanDAT_TS
  
  m2 <- runDFA(zMeanDWT_TS,NumStates=2,ErrStruc='DUE',EstCovar=T,indivCovar=T,Covars = zMeanDAT_TS)
  m2$ydat <- zMeanDWT_TS
  m2$xdat <- zMeanDAT_TS
  
  m3 <- runDFA(zMeanDWT_TS,NumStates=3,ErrStruc='DUE',EstCovar=T,indivCovar=T,Covars = zMeanDAT_TS)
  m3$ydat <- zMeanDWT_TS
  m3$xdat <- zMeanDAT_TS
  
  list(m1,m2,m3)
}

#End cluster
stopCluster(cl)

#Extract AIC for multiple trends

lapply(1:length(Subset1_DFA),FUN=function(x){
  SS<-Subset1_DFA[[x]]
  AIC<-lapply(SS,FUN=function(y){
    return(y$AIC)
  }) %>% bind_cols()
  colnames(AIC) <- c('T1','T2','T3')
  return(AIC)
}) %>% bind_rows() %>% apply(1,FUN=function(xx){return(xx-min(xx))}) %>% t()

#Build model result tables
##### ALL SUBSETS
lapply(1:3,FUN=function(sub.ind){
  #Set subset data.frame
  SS <- switch(sub.ind,Subset1_DFA,Subset2_DFA,Subset3_DFA)
  
  GlobalResults <- lapply(1:length(SS),FUN=function(x){
    #x<-1
    Sites <- SS[[x]][[1]]$ydat %>% row.names()
    Year <- rep(2010 + x,length(Sites))
    TempSens <- lapply(SS[[x]],FUN=function(yy){
      diag(yy$Estimates$D) %>% as.vector()
    }) %>% bind_cols()
    colnames(TempSens) <- paste0('TempSens_',seq(1,3))
    
    TrendLoad <- lapply(SS[[x]],FUN=function(yy){
      #yy<-Subset1_DFA[[x]][[2]]
      TLoad<-yy$Estimates$Z
      TLoad_out <- lapply(1:ncol(TLoad),FUN=function(zz){return(data.frame(TLoad[,zz]))})
      return(TLoad_out)
    }) %>% bind_cols()
    colnames(TrendLoad) <- paste0('TrendLoad_',c(1,rep(2,2),rep(3,3)),rep('_',6),c(1,1:2,1:3))
    
    
    return(data.frame(SiteID = Sites,Year,TempSens,TrendLoad))
  }) %>% bind_rows()
  GlobalResults$Region <- Tdat$Region[match(GlobalResults$SiteID,Tdat$SiteID)]
  #Extract common trends from each year
  GlobalTrends <- lapply(1:length(SS),FUN=function(x){
    #x<-1
    DOY <- SS[[x]][[1]]$ydat %>% colnames() %>% as.numeric()
    Year <- rep(2010+x,length(DOY))
    Trend <- lapply(SS[[x]],FUN=function(yy){
      U<-yy$Estimates$u
      U_out<-lapply(1:nrow(U),FUN=function(zz){return(data.frame(U[zz,]))})
      return(U_out)
    }) %>% bind_cols()
    colnames(Trend) <- paste0('Trend_',c(1,rep(2,2),rep(3,3)),rep('_',6),c(1,1:2,1:3))
    return(data.frame(Year=Year,DOY=DOY,Trend))
  }) %>% bind_rows()
  
  assign(paste0('GlobalResults_Subset',sub.ind),value=GlobalResults,envir=.GlobalEnv)
  assign(paste0('GlobalTrends_Subset',sub.ind),value=GlobalTrends,envir=.GlobalEnv)
  return(paste0('Subset ',sub.ind,' has written to Global environment'))
})


lapply(1:3,FUN=function(sub.ind){
  GlobalResults <- switch(sub.ind,GlobalResults_Subset1,GlobalResults_Subset2,GlobalResults_Subset3)
  GlobalTrends <- switch(sub.ind,GlobalTrends_Subset1,GlobalTrends_Subset2,GlobalTrends_Subset3)
  Ntrends <- 3
  
  pdf(paste0('DFAresults_Subset',sub.ind,'.pdf'),width=8.5,height=11)
  
  lapply(1:Ntrends,FUN=function(TRD){
    TLs <- switch(TRD,GlobalResults %>% select(TrendLoad_1_1),GlobalResults %>% select(TrendLoad_2_1,TrendLoad_2_2),GlobalResults %>% select(TrendLoad_3_1,TrendLoad_3_2,TrendLoad_3_3))
    Ts <- switch(TRD,GlobalTrends %>% select(Trend_1_1),GlobalTrends %>% select(Trend_2_1,Trend_2_2),GlobalTrends %>% select(Trend_3_1,Trend_3_2,Trend_3_3))
    
    UniYear <- sort(unique(GlobalResults$Year))
    Nyear <- length(UniYear)
     lapply(1:Nyear,FUN=function(yr){
       #yr<-1
       TL_y <- TLs %>% filter(GlobalResults$Year == UniYear[yr])
       T_y <- Ts %>% filter(GlobalTrends$Year == UniYear[yr])
       DOY_y <- GlobalTrends %>% filter(GlobalTrends$Year == UniYear[yr]) %>% pull(DOY)
       ColChoice <- c('#000000','#0000FF','#FF0000')
       if(yr <= 5){pyr <- yr}else{pyr <- yr - 5}
       
       if(yr %in% c(1,6)){
         par(fig=c(0,0.6,1-0.2*pyr,1-0.2*(pyr-1)))
       }else{
         par(fig=c(0,0.6,1-0.2*pyr,1-0.2*(pyr-1)),new=T)
       }
       par(mar=c(1,1,1,1),oma=c(3,3,0,0))
       #Plot Trends
       plot(DOY_y,T_y[,1],type='l',col=ColChoice[1],lwd=2,ylim=c(-3,3),xlab='DOY',ylab='Trend')
       if(TRD>1){
         for(TT in 2:ncol(T_y)){
           points(DOY_y,T_y[,TT],type='l',col=ColChoice[TT],lwd=2)
         }
       }
       abline(h=0,lty=3)
       mtext(UniYear[yr],3,line=-1.5,cex=1)
       
       #Plot loadings histograms
       par(fig=c(0.6,1,1-0.2*pyr,1-0.2*(pyr-1)),new=T)
       plot(density(TL_y[,1]),type='l',col=ColChoice[1],lwd=2,xlim=c(-3,3),main='')
       polygon(density(TL_y[,1]), col=paste0(ColChoice[1],'55'), border=ColChoice[1])
       if(TRD>1){
         for(TT in 2:ncol(T_y)){
           points(density(TL_y[,TT]),type='l',col=ColChoice[TT],lwd=2)
           polygon(density(TL_y[,TT]), col=paste0(ColChoice[TT],'55'), border=ColChoice[TT])
         }
       }
       abline(v=0,lty=3)
       
       
     }) #End Plotting over Year
     
     
     })#End Apply over trends
  dev.off() #Turnoff PDF
}) # End apply over subsets


lapply(1:3,FUN=function(sub.ind){
  GlobalResults <- switch(sub.ind,GlobalResults_Subset1,GlobalResults_Subset2,GlobalResults_Subset3)  GlobalTrends <- switch(sub.ind,GlobalTrends_Subset1,GlobalTrends_Subset2,GlobalTrends_Subset3)
  Ntrends <- 3
  
  pdf(paste0('DFAresults_Subset',sub.ind,'_withTempSens.pdf'),width=8.5,height=11)
  
  lapply(1:Ntrends,FUN=function(TRD){
    TLs <- switch(TRD,GlobalResults %>% select(TrendLoad_1_1),GlobalResults %>% select(TrendLoad_2_1,TrendLoad_2_2),GlobalResults %>% select(TrendLoad_3_1,TrendLoad_3_2,TrendLoad_3_3))
    Ts <- switch(TRD,GlobalTrends %>% select(Trend_1_1),GlobalTrends %>% select(Trend_2_1,Trend_2_2),GlobalTrends %>% select(Trend_3_1,Trend_3_2,Trend_3_3))
    
    UniYear <- sort(unique(GlobalResults$Year))
    Nyear <- length(UniYear)
    lapply(1:Nyear,FUN=function(yr){
      #yr<-1
      TL_y <- TLs %>% filter(GlobalResults$Year == UniYear[yr])
      T_y <- Ts %>% filter(GlobalTrends$Year == UniYear[yr])
      DOY_y <- GlobalTrends %>% filter(GlobalTrends$Year == UniYear[yr]) %>% pull(DOY)
      ColChoice <- c('#000000','#0000FF','#FF0000')
      if(yr <= 5){pyr <- yr}else{pyr <- yr - 5}
      
      if(yr %in% c(1,6)){
        par(fig=c(0,0.6,1-0.2*pyr,1-0.2*(pyr-1)))
      }else{
        par(fig=c(0,0.6,1-0.2*pyr,1-0.2*(pyr-1)),new=T)
      }
      par(mar=c(1,1,1,1),oma=c(3,3,0,0))
      #Plot Trends
      plot(DOY_y,T_y[,1],type='l',col=ColChoice[1],lwd=2,ylim=c(-3,3),xlab='DOY',ylab='Trend')
      if(TRD>1){
        for(TT in 2:ncol(T_y)){
          points(DOY_y,T_y[,TT],type='l',col=ColChoice[TT],lwd=2)
        }
      }
      abline(h=0,lty=3)
      mtext(UniYear[yr],3,line=-1.5,cex=1)
      
      #Plot loadings histograms
      par(fig=c(0.6,1,1-0.2*pyr,1-0.2*(pyr-1)),new=T)
      plot(density(TL_y[,1]),type='l',col=ColChoice[1],lwd=2,xlim=c(-3,3),main='')
      polygon(density(TL_y[,1]), col=paste0(ColChoice[1],'55'), border=ColChoice[1])
      if(TRD>1){
        for(TT in 2:ncol(T_y)){
          points(density(TL_y[,TT]),type='l',col=ColChoice[TT],lwd=2)
          polygon(density(TL_y[,TT]), col=paste0(ColChoice[TT],'55'), border=ColChoice[TT])
        }
      }
      abline(v=0,lty=3)
      
      
    }) #End Plotting over Year
    
    
  })#End Apply over trends
  dev.off() #Turnoff PDF
}) # End apply over subsets


save(list=ls(),file = 'DFAresults_AllSubsets.Rdata')


# Trends_TS <- xtabs(Trend ~ Year + DOY,data=GlobalTrends)
# Trends_TS[Trends_TS==0] <- NA
# 
# library(RColorBrewer)
# Cols <- brewer.pal(9,'Set1')
# plot(Trends_TS[1,],type='l',lwd=2,col='black',ylim=c(-3,3))
# for(i in 2:10){
#   points(Trends_TS[i,],type='l',lwd=2,col=Cols[i-1])
# }
# 
# plot(GlobalResults$TrendLoad)
# 
# lmer1<-lmer(TempSens ~ TrendLoad + (1|Region) + (1|Year) + (1|SiteID),data=GlobalResults)
# summary(lmer1)
# 
# ranef(lmer1)

######
######
######
###### Below is code to run DFA's by region. Just in case we want it.

n.cores <- parallel::detectCores()
cl <- makeCluster(n.cores)
registerDoParallel(cl)

UniModels <- expand.grid(unique(Tdat$Region),seq(2011,2020))

RegionalDFA <- foreach(y = 1:nrow(UniModels),.packages=c('dplyr','TMB')) %dopar% {
  dyn.load(dynlib("dfa1tmb"))
  #y<-10
  MeanDWT_TS <- Tdat %>% filter(Year==UniModels[y,2],Region == UniModels[y,1])
  if(length(unique(MeanDWT_TS$SiteID))>5){
    MeanDWT_TS <- MeanDWT_TS %>% xtabs(meanDT ~ SiteID + DOY,data=.)
    MeanDWT_TS[MeanDWT_TS==0] <- NA #xtabs puts zeros in for missing dates in the dataset, we set those exact zeros to NA
  
    #keep sites with >80% of data
    MeanDWT_TS <- MeanDWT_TS[which(rowSums(!is.na(MeanDWT_TS)) >= round(0.8*length(DOYseq)) ),]
  
  
    #build time series of matching air temp
    MeanDAT_TS <- Tdat %>% filter(Year==UniModels[y,2], SiteID %in% row.names(MeanDWT_TS)) %>% xtabs(airDT ~ SiteID + DOY,data=.)
    
    #Zscore response and air temp covariate
    zMeanDWT_TS <- apply(MeanDWT_TS,1,FUN=Zscore) %>% t()
    zMeanDAT_TS <- apply(MeanDAT_TS,1,FUN=Zscore) %>% t()
    
    m1 <- runDFA(zMeanDWT_TS,NumStates=1,ErrStruc='DUE',EstCovar=T,indivCovar=T,Covars = zMeanDAT_TS)
    m1$ydat <- zMeanDWT_TS
    m1$xdat <- zMeanDAT_TS
    
  }else{
    m1 <- NA
  }
  m1
}

names(RegionalDFA) <- lapply(1:nrow(UniModels),FUN=function(x){paste(UniModels[x,1],UniModels[x,2],sep='.')})

#End cluster
stopCluster(cl)

library(stringr)

RegionalResults <- lapply(1:length(RegionalDFA),FUN=function(x){
  if(!is.na(RegionalDFA[[x]][1])){
    #print(x)
    Sites <- RegionalDFA[[x]]$ydat %>% row.names()
    TempSens <- diag(RegionalDFA[[x]]$Estimates$D)
    TrendLoad <- RegionalDFA[[x]]$Estimates$Z[,1]
    
    Year <- rep(str_split(names(RegionalDFA)[x],pattern='\\.')[[1]][2] %>% as.numeric(),length(Sites))
    Region <- rep(str_split(names(RegionalDFA)[x],pattern='\\.')[[1]][1],length(Sites))
    return(data.frame(Region, SiteID = Sites,Year,TempSens,TrendLoad))
  }else{
    return(data.frame(Region=NA,SiteID=NA,Year=NA,TempSens=NA,TrendLoad=NA))
  }
}) %>% bind_rows() %>% na.omit()

RegionalTrends <- lapply(1:length(RegionalDFA),FUN=function(x){
  if(!is.na(RegionalDFA[[x]][1])){
    DOY <- RegionalDFA[[x]]$ydat %>% colnames() %>% as.numeric()
    Year <- rep(str_split(names(RegionalDFA)[x],pattern='\\.')[[1]][2],length(DOY))
    Region <- rep(str_split(names(RegionalDFA)[x],pattern='\\.')[[1]][1],length(DOY))
    Trend <- RegionalDFA[[x]]$Estimates$u %>% as.vector()
  
    return(data.frame(Region,Year,DOY,Trend))
  }else{
    return(data.frame(Region=NA,Year=NA,DOY=NA,Trend=NA))
  }
}) %>% bind_rows() %>% na.omit()


MergedModels <- merge(GlobalResults,RegionalResults,by=c('SiteID','Year'))
# 
plot(MergedModels$TempSens.x,MergedModels$TempSens.y)
summary(lm(MergedModels$TempSens.x~MergedModels$TempSens.y))

TodaysDate <- paste(str_split(as.Date(Sys.time()),pattern='-')[[1]],collapse='.')

save(list=ls(),file=paste0('DFAresults_global_and_regional_',TodaysDate,'.Rdata'))

