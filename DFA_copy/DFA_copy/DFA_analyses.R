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
MetaDat <- readRDS('../data/md.rds')

#Sequence of days between June 1 and Aug 31
DOYseq<-seq(152,243)

#load stream temperature data
#create some new date/time columns
#add region column
#limit timeframe to chosen date range
Tdat <- readRDS('../../data/summer_data_wair2021-11-23.rds') %>% 
  mutate(sampleDate = as.Date(sampleDate), DOY = format(sampleDate,'%j') %>% as.numeric(),Year=format(sampleDate,'%Y')%>% as.numeric()) %>%
  mutate(Region = MetaDat$Region[match(SiteID,MetaDat$SiteID)]) %>%
  filter(DOY %in% DOYseq)


which(is.na(Tdat$airDT))
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

#setup cluster for parallel processing
n.cores <- parallel::detectCores()
cl <- makeCluster(n.cores)
registerDoParallel(cl)

#Parallel loop to run all years and store in a list
#Each parallel instance is like its own R session
#We have to pass the proper packages to run the functions and load dynamic model codes 
GlobalDFA <- foreach(y = 2011:2020,.packages=c('dplyr','TMB')) %dopar% {
  dyn.load(dynlib("dfa1tmb"))
  MeanDWT_TS <- Tdat %>% filter(Year==y) %>% xtabs(meanDT ~ SiteID + DOY,data=.)
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
  
  m1
  
}

#End cluster
stopCluster(cl)

#Build model result tables
GlobalResults <- lapply(1:length(GlobalDFA),FUN=function(x){
  #x<-1
  Sites <- GlobalDFA[[x]]$ydat %>% row.names()
  TempSens <- diag(GlobalDFA[[x]]$Estimates$D)
  TrendLoad <- GlobalDFA[[x]]$Estimates$Z[,1]
  Year <- rep(2010 + x,length(Sites))
  
  data.frame(SiteID = Sites,Year=Year,TempSens=TempSens,TrendLoad)
}) %>% bind_rows()
GlobalResults$Region <- Tdat$Region[match(GlobalResults$SiteID,Tdat$SiteID)]

#Extract common trends from each year
GlobalTrends <- lapply(1:length(GlobalDFA),FUN=function(x){
  #x<-1
  DOY <- GlobalDFA[[x]]$ydat %>% colnames() %>% as.numeric()
  Year <- rep(2010+x,length(DOY))
  Trend <- GlobalDFA[[x]]$Estimates$u %>% as.vector()
  
  data.frame(Year=Year,DOY=DOY,Trend=Trend)
}) %>% bind_rows()

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




