rm(list=ls())
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
Tdat <- read.csv('../../data/summer_data_wair_dayl2022-03-17.csv',header=T,stringsAsFactors=F) %>% 
  mutate(sampleDate = as.Date(sampleDate), DOY = format(sampleDate,'%j') %>% as.numeric(),Year=format(sampleDate,'%Y')%>% as.numeric()) %>%
  mutate(Region = MetaDat$Region[match(SiteID,MetaDat$Site)]) %>%
  filter(DOY %in% DOYseq)

head(Tdat)

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
  #y<-2019
  dyn.load(dynlib("dfa1tmb"))
  Sub1 <- Tdat %>% filter(SiteID %in% Sub1_Sites)
  MeanDWT_TS <- Sub1 %>% filter(Year==y) %>% xtabs(meanDT ~ SiteID + DOY,data=.)
  MeanDWT_TS[MeanDWT_TS==0] <- NA #xtabs puts zeros in for missing dates in the dataset, we set those exact zeros to NA

  #keep sites with >80% of data
  MeanDWT_TS <- MeanDWT_TS[which(rowSums(!is.na(MeanDWT_TS)) >= round(0.8*length(DOYseq)) ),]
  
  #build time series of matching daylength variable
  DayLength_TS <- Tdat %>% filter(Year==y,SiteID %in% row.names(MeanDWT_TS)) %>% xtabs(global_dayl ~ SiteID + DOY,data=.)
  #build time series of matching air temp
  MeanDAT_TS <- Tdat %>% filter(Year==y, SiteID %in% row.names(MeanDWT_TS)) %>% xtabs(airDT ~ SiteID + DOY,data=.)
  
  #Zscore tables
  zMeanDWT_TS<-apply(MeanDWT_TS,1,FUN=Zscore) %>% t()
  zDayLength_TS<-apply(DayLength_TS,1,FUN=Zscore) %>% t()
  zMeanDAT_TS<-apply(MeanDAT_TS,1,FUN=Zscore) %>% t()
  
  m1.air <- runDFA(zMeanDWT_TS,NumStates=1,ErrStruc='DUE',EstCovar=T,indivCovar=T,Covars = zMeanDAT_TS)
  m1.air$ydat <- zMeanDWT_TS
  m1.air$xdat <- zMeanDAT_TS
  
  m1.day <- runDFA(zMeanDWT_TS,NumStates=1,ErrStruc='DUE',EstCovar=T,indivCovar=T,Covars = zDayLength_TS)
  m1.day$ydat <- zMeanDWT_TS
  m1.day$xdat <- zDayLength_TS
  
  #Custom D parameter matrix for two covariates per site
  Dp1 <- matrix(0,nrow=nrow(zMeanDWT_TS),ncol=nrow(zMeanDWT_TS))
  Df1 <- diag(x=1:nrow(zMeanDWT_TS))
  Df2 <- diag(x=(nrow(zMeanDWT_TS)+1):(2*nrow(zMeanDWT_TS)))
  if(!identical(dim(Dp1),dim(Df1))){break;print('Error')}
  Dmat2 <- cbind(Dp1,Dp1)
  Dfac2 <- cbind(Df1,Df2)
  Dfac2[Dfac2 == 0] <- NA
  Dfac2 <- as.factor(Dfac2)
  
  m1.both <- runDFA(zMeanDWT_TS,NumStates=1,ErrStruc='DUE',EstCovar=T,indivCovar=T,Dmat=Dmat2,Dfac=Dfac2,Covars = rbind(zMeanDAT_TS,zDayLength_TS))
  m1.both$ydat <- zMeanDWT_TS
  m1.both$xdat <- rbind(zMeanDAT_TS,zDayLength_TS)
  
  #Custom D parameter matrix for three covariates per site (air, day, interaction)
  Dp1 <- matrix(0,nrow=nrow(zMeanDWT_TS),ncol=nrow(zMeanDWT_TS))
  Df1 <- diag(x=1:nrow(zMeanDWT_TS))
  Df2 <- diag(x=(nrow(zMeanDWT_TS)+1):(2*nrow(zMeanDWT_TS)))
  Df3 <- diag(x=(2*nrow(zMeanDWT_TS)+1):(3*nrow(zMeanDWT_TS)))
  if(!identical(dim(Dp1),dim(Df1))){break;print('Error')}
  Dmat2 <- cbind(Dp1,Dp1,Dp1)
  Dfac2 <- cbind(Df1,Df2,Df3)
  Dfac2[Dfac2 == 0] <- NA
  Dfac2 <- as.factor(Dfac2)
  
  IntVar <- zMeanDAT_TS*zDayLength_TS #apply(MeanDAT_TS*DayLength_TS,1,FUN=Zscore) %>% t()
  # if(F){index<-52
  # plot(as.numeric(colnames(MeanDAT_TS)),zMeanDAT_TS[index,],main=row.names(MeanDAT_TS)[index],type='l',col='blue',ylab='air temp',xlab='DOY')
  # points(as.numeric(colnames(MeanDAT_TS)),zDayLength_TS[index,],type='l',col='red')
  # points(as.numeric(colnames(MeanDAT_TS)),IntVar[index,],type='l',col='purple')}
  
  m1.int <- runDFA(zMeanDWT_TS,NumStates=1,ErrStruc='DUE',EstCovar=T,indivCovar=T,Dmat=Dmat2,Dfac=Dfac2,Covars = rbind(zMeanDAT_TS,zDayLength_TS,IntVar))
  m1.int$ydat <- zMeanDWT_TS
  m1.int$xdat <- rbind(zMeanDAT_TS,zDayLength_TS,IntVar)
  
  #lapply(1:nrow(zMeanDAT_TS),FUN=function(x){cor(zMeanDAT_TS[x,],zDayLength_TS[x,])}) %>% unlist() %>% hist()
  
  # m1.air$AIC
  # m1.day$AIC
  # m1.both$AIC
  # m1.int$AIC
  
  # summary(lm(as.vector(zMeanDWT_TS) ~ as.vector(m1.air$Fits)))
  # summary(lm(as.vector(zMeanDWT_TS) ~ as.vector(m1.day$Fits)))
  # summary(lm(as.vector(zMeanDWT_TS) ~ as.vector(m1.both$Fits)))
  # summary(lm(as.vector(zMeanDWT_TS) ~ as.vector(m1.int$Fits)))
  # 
  # hist(diag(m1.int$Estimates$D[,-c(1:(2*nrow(zMeanDWT_TS)))]))
  # diag(m1.int$Estimates$D[,-c(1:(2*nrow(zMeanDWT_TS)))]) %>% mean()
  # lapply(1:nrow(zMeanDWT_TS),FUN=function(i){
  #   summary(lm(zMeanDWT_TS[i,] ~ m1.int$Fits[i,]))$r.squared
  # }) %>% unlist() %>% hist()
  
  
  
  #m1.air$Optimization
  #m1.both$Optimization
  
  list(m1.air=m1.air,m1.day=m1.day,m1.both=m1.both,m1.int=m1.int)
}

Subset2_DFA <- foreach(y = 2011:2020,.packages=c('dplyr','TMB')) %dopar% {
  #y<-2019
  dyn.load(dynlib("dfa1tmb"))
  Sub1 <- Tdat %>% filter(SiteID %in% Sub2_Sites)
  MeanDWT_TS <- Sub1 %>% filter(Year==y) %>% xtabs(meanDT ~ SiteID + DOY,data=.)
  MeanDWT_TS[MeanDWT_TS==0] <- NA #xtabs puts zeros in for missing dates in the dataset, we set those exact zeros to NA
  
  #keep sites with >80% of data
  MeanDWT_TS <- MeanDWT_TS[which(rowSums(!is.na(MeanDWT_TS)) >= round(0.8*length(DOYseq)) ),]
  
  #build time series of matching daylength variable
  DayLength_TS <- Tdat %>% filter(Year==y,SiteID %in% row.names(MeanDWT_TS)) %>% xtabs(global_dayl ~ SiteID + DOY,data=.)
  #build time series of matching air temp
  MeanDAT_TS <- Tdat %>% filter(Year==y, SiteID %in% row.names(MeanDWT_TS)) %>% xtabs(airDT ~ SiteID + DOY,data=.)
  
  #Zscore tables
  zMeanDWT_TS<-apply(MeanDWT_TS,1,FUN=Zscore) %>% t()
  zDayLength_TS<-apply(DayLength_TS,1,FUN=Zscore) %>% t()
  zMeanDAT_TS<-apply(MeanDAT_TS,1,FUN=Zscore) %>% t()
  
  m1.air <- runDFA(zMeanDWT_TS,NumStates=1,ErrStruc='DUE',EstCovar=T,indivCovar=T,Covars = zMeanDAT_TS)
  m1.air$ydat <- zMeanDWT_TS
  m1.air$xdat <- zMeanDAT_TS
  
  m1.day <- runDFA(zMeanDWT_TS,NumStates=1,ErrStruc='DUE',EstCovar=T,indivCovar=T,Covars = zDayLength_TS)
  m1.day$ydat <- zMeanDWT_TS
  m1.day$xdat <- zDayLength_TS
  
  #Custom D parameter matrix for two covariates per site
  Dp1 <- matrix(0,nrow=nrow(zMeanDWT_TS),ncol=nrow(zMeanDWT_TS))
  Df1 <- diag(x=1:nrow(zMeanDWT_TS))
  Df2 <- diag(x=(nrow(zMeanDWT_TS)+1):(2*nrow(zMeanDWT_TS)))
  if(!identical(dim(Dp1),dim(Df1))){break;print('Error')}
  Dmat2 <- cbind(Dp1,Dp1)
  Dfac2 <- cbind(Df1,Df2)
  Dfac2[Dfac2 == 0] <- NA
  Dfac2 <- as.factor(Dfac2)
  
  m1.both <- runDFA(zMeanDWT_TS,NumStates=1,ErrStruc='DUE',EstCovar=T,indivCovar=T,Dmat=Dmat2,Dfac=Dfac2,Covars = rbind(zMeanDAT_TS,zDayLength_TS))
  m1.both$ydat <- zMeanDWT_TS
  m1.both$xdat <- rbind(zMeanDAT_TS,zDayLength_TS)
  
  #Custom D parameter matrix for three covariates per site (air, day, interaction)
  Dp1 <- matrix(0,nrow=nrow(zMeanDWT_TS),ncol=nrow(zMeanDWT_TS))
  Df1 <- diag(x=1:nrow(zMeanDWT_TS))
  Df2 <- diag(x=(nrow(zMeanDWT_TS)+1):(2*nrow(zMeanDWT_TS)))
  Df3 <- diag(x=(2*nrow(zMeanDWT_TS)+1):(3*nrow(zMeanDWT_TS)))
  if(!identical(dim(Dp1),dim(Df1))){break;print('Error')}
  Dmat2 <- cbind(Dp1,Dp1,Dp1)
  Dfac2 <- cbind(Df1,Df2,Df3)
  Dfac2[Dfac2 == 0] <- NA
  Dfac2 <- as.factor(Dfac2)
  
  IntVar <- zMeanDAT_TS*zDayLength_TS #apply(MeanDAT_TS*DayLength_TS,1,FUN=Zscore) %>% t()
  # if(F){index<-52
  # plot(as.numeric(colnames(MeanDAT_TS)),zMeanDAT_TS[index,],main=row.names(MeanDAT_TS)[index],type='l',col='blue',ylab='air temp',xlab='DOY')
  # points(as.numeric(colnames(MeanDAT_TS)),zDayLength_TS[index,],type='l',col='red')
  # points(as.numeric(colnames(MeanDAT_TS)),IntVar[index,],type='l',col='purple')}
  
  m1.int <- runDFA(zMeanDWT_TS,NumStates=1,ErrStruc='DUE',EstCovar=T,indivCovar=T,Dmat=Dmat2,Dfac=Dfac2,Covars = rbind(zMeanDAT_TS,zDayLength_TS,IntVar))
  m1.int$ydat <- zMeanDWT_TS
  m1.int$xdat <- rbind(zMeanDAT_TS,zDayLength_TS,IntVar)
  
  #lapply(1:nrow(zMeanDAT_TS),FUN=function(x){cor(zMeanDAT_TS[x,],zDayLength_TS[x,])}) %>% unlist() %>% hist()
  
  # m1.air$AIC
  # m1.day$AIC
  # m1.both$AIC
  # m1.int$AIC
  
  # summary(lm(as.vector(zMeanDWT_TS) ~ as.vector(m1.air$Fits)))
  # summary(lm(as.vector(zMeanDWT_TS) ~ as.vector(m1.day$Fits)))
  # summary(lm(as.vector(zMeanDWT_TS) ~ as.vector(m1.both$Fits)))
  # summary(lm(as.vector(zMeanDWT_TS) ~ as.vector(m1.int$Fits)))
  # 
  # hist(diag(m1.int$Estimates$D[,-c(1:(2*nrow(zMeanDWT_TS)))]))
  # diag(m1.int$Estimates$D[,-c(1:(2*nrow(zMeanDWT_TS)))]) %>% mean()
  # lapply(1:nrow(zMeanDWT_TS),FUN=function(i){
  #   summary(lm(zMeanDWT_TS[i,] ~ m1.int$Fits[i,]))$r.squared
  # }) %>% unlist() %>% hist()
  
  
  
  #m1.air$Optimization
  #m1.both$Optimization
  
  list(m1.air=m1.air,m1.day=m1.day,m1.both=m1.both,m1.int=m1.int)
}

Subset3_DFA <- foreach(y = 2011:2020,.packages=c('dplyr','TMB')) %dopar% {
  #y<-2019
  dyn.load(dynlib("dfa1tmb"))
  Sub1 <- Tdat %>% filter(SiteID %in% Sub3_Sites)
  MeanDWT_TS <- Sub1 %>% filter(Year==y) %>% xtabs(meanDT ~ SiteID + DOY,data=.)
  MeanDWT_TS[MeanDWT_TS==0] <- NA #xtabs puts zeros in for missing dates in the dataset, we set those exact zeros to NA
  
  #keep sites with >80% of data
  MeanDWT_TS <- MeanDWT_TS[which(rowSums(!is.na(MeanDWT_TS)) >= round(0.8*length(DOYseq)) ),]
  
  #build time series of matching daylength variable
  DayLength_TS <- Tdat %>% filter(Year==y,SiteID %in% row.names(MeanDWT_TS)) %>% xtabs(global_dayl ~ SiteID + DOY,data=.)
  #build time series of matching air temp
  MeanDAT_TS <- Tdat %>% filter(Year==y, SiteID %in% row.names(MeanDWT_TS)) %>% xtabs(airDT ~ SiteID + DOY,data=.)
  
  #Zscore tables
  zMeanDWT_TS<-apply(MeanDWT_TS,1,FUN=Zscore) %>% t()
  zDayLength_TS<-apply(DayLength_TS,1,FUN=Zscore) %>% t()
  zMeanDAT_TS<-apply(MeanDAT_TS,1,FUN=Zscore) %>% t()
  
  m1.air <- runDFA(zMeanDWT_TS,NumStates=1,ErrStruc='DUE',EstCovar=T,indivCovar=T,Covars = zMeanDAT_TS)
  m1.air$ydat <- zMeanDWT_TS
  m1.air$xdat <- zMeanDAT_TS
  
  m1.day <- runDFA(zMeanDWT_TS,NumStates=1,ErrStruc='DUE',EstCovar=T,indivCovar=T,Covars = zDayLength_TS)
  m1.day$ydat <- zMeanDWT_TS
  m1.day$xdat <- zDayLength_TS
  
  #Custom D parameter matrix for two covariates per site
  Dp1 <- matrix(0,nrow=nrow(zMeanDWT_TS),ncol=nrow(zMeanDWT_TS))
  Df1 <- diag(x=1:nrow(zMeanDWT_TS))
  Df2 <- diag(x=(nrow(zMeanDWT_TS)+1):(2*nrow(zMeanDWT_TS)))
  if(!identical(dim(Dp1),dim(Df1))){break;print('Error')}
  Dmat2 <- cbind(Dp1,Dp1)
  Dfac2 <- cbind(Df1,Df2)
  Dfac2[Dfac2 == 0] <- NA
  Dfac2 <- as.factor(Dfac2)
  
  m1.both <- runDFA(zMeanDWT_TS,NumStates=1,ErrStruc='DUE',EstCovar=T,indivCovar=T,Dmat=Dmat2,Dfac=Dfac2,Covars = rbind(zMeanDAT_TS,zDayLength_TS))
  m1.both$ydat <- zMeanDWT_TS
  m1.both$xdat <- rbind(zMeanDAT_TS,zDayLength_TS)
  
  #Custom D parameter matrix for three covariates per site (air, day, interaction)
  Dp1 <- matrix(0,nrow=nrow(zMeanDWT_TS),ncol=nrow(zMeanDWT_TS))
  Df1 <- diag(x=1:nrow(zMeanDWT_TS))
  Df2 <- diag(x=(nrow(zMeanDWT_TS)+1):(2*nrow(zMeanDWT_TS)))
  Df3 <- diag(x=(2*nrow(zMeanDWT_TS)+1):(3*nrow(zMeanDWT_TS)))
  if(!identical(dim(Dp1),dim(Df1))){break;print('Error')}
  Dmat2 <- cbind(Dp1,Dp1,Dp1)
  Dfac2 <- cbind(Df1,Df2,Df3)
  Dfac2[Dfac2 == 0] <- NA
  Dfac2 <- as.factor(Dfac2)
  
  IntVar <- zMeanDAT_TS*zDayLength_TS #apply(MeanDAT_TS*DayLength_TS,1,FUN=Zscore) %>% t()
  # if(F){index<-52
  # plot(as.numeric(colnames(MeanDAT_TS)),zMeanDAT_TS[index,],main=row.names(MeanDAT_TS)[index],type='l',col='blue',ylab='air temp',xlab='DOY')
  # points(as.numeric(colnames(MeanDAT_TS)),zDayLength_TS[index,],type='l',col='red')
  # points(as.numeric(colnames(MeanDAT_TS)),IntVar[index,],type='l',col='purple')}
  
  m1.int <- runDFA(zMeanDWT_TS,NumStates=1,ErrStruc='DUE',EstCovar=T,indivCovar=T,Dmat=Dmat2,Dfac=Dfac2,Covars = rbind(zMeanDAT_TS,zDayLength_TS,IntVar))
  m1.int$ydat <- zMeanDWT_TS
  m1.int$xdat <- rbind(zMeanDAT_TS,zDayLength_TS,IntVar)
  
  #lapply(1:nrow(zMeanDAT_TS),FUN=function(x){cor(zMeanDAT_TS[x,],zDayLength_TS[x,])}) %>% unlist() %>% hist()
  
  # m1.air$AIC
  # m1.day$AIC
  # m1.both$AIC
  # m1.int$AIC
  
  # summary(lm(as.vector(zMeanDWT_TS) ~ as.vector(m1.air$Fits)))
  # summary(lm(as.vector(zMeanDWT_TS) ~ as.vector(m1.day$Fits)))
  # summary(lm(as.vector(zMeanDWT_TS) ~ as.vector(m1.both$Fits)))
  # summary(lm(as.vector(zMeanDWT_TS) ~ as.vector(m1.int$Fits)))
  # 
  # hist(diag(m1.int$Estimates$D[,-c(1:(2*nrow(zMeanDWT_TS)))]))
  # diag(m1.int$Estimates$D[,-c(1:(2*nrow(zMeanDWT_TS)))]) %>% mean()
  # lapply(1:nrow(zMeanDWT_TS),FUN=function(i){
  #   summary(lm(zMeanDWT_TS[i,] ~ m1.int$Fits[i,]))$r.squared
  # }) %>% unlist() %>% hist()
  
  
  
  #m1.air$Optimization
  #m1.both$Optimization
  
  list(m1.air=m1.air,m1.day=m1.day,m1.both=m1.both,m1.int=m1.int)
}

#End cluster
stopCluster(cl)

#Extract AIC for multiple covariates
lapply(1:length(Subset1_DFA),FUN=function(x){
  SS<-Subset1_DFA[[x]]
  AIC<-lapply(SS,FUN=function(y){
    return(y$AIC)
  }) %>% bind_cols()
  colnames(AIC) <- c('Air','DayLen','Air.DayLen','Air.DayLen.Int')
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
    
    TempSens <- lapply(1:length(SS[[x]]),FUN=function(yy){
      YY<-SS[[x]][[yy]]
      if(names(SS[[x]])[yy]!='m1.day'){
        diag(YY$Estimates$D[1:length(Sites),1:length(Sites)]) %>% as.vector()
      }else{
        rep(NA,length(Sites))
      }
    }) %>% bind_cols()
    colnames(TempSens) <- paste0('TempSensZ_',c('Air','DayLen','Air.DayLen','Air.DayLen.Int'))
    TempSens <- TempSens %>% select(TempSensZ_Air,TempSensZ_Air.DayLen,TempSensZ_Air.DayLen.Int)
    
    DayLenSens <- lapply(1:length(SS[[x]]),FUN=function(yy){
      YY<-SS[[x]][[yy]]
      if(names(SS[[x]])[yy]=='m1.day'){
        diag(YY$Estimates$D[1:length(Sites),1:length(Sites)]) %>% as.vector()
      }else if(names(SS[[x]])[yy] %in% c('m1.both','m1.int')){
        diag(YY$Estimates$D[1:length(Sites),(length(Sites)+1):(2*length(Sites))]) %>% as.vector()
      }else{  
        rep(NA,length(Sites))
      }
    }) %>% bind_cols()
    colnames(DayLenSens) <- paste0('DayLenSensZ_',c('Air','DayLen','Air.DayLen','Air.DayLen.Int'))
    DayLenSens <- DayLenSens %>% select(DayLenSensZ_DayLen,DayLenSensZ_Air.DayLen,DayLenSensZ_Air.DayLen.Int)
    
    Interaction <- lapply(1:length(SS[[x]]),FUN=function(yy){
      YY<-SS[[x]][[yy]]
      if(names(SS[[x]])[yy]=='m1.int'){
        diag(YY$Estimates$D[1:length(Sites),(2*length(Sites)+1):(3*length(Sites))]) %>% as.vector()
      }else{  
        rep(NA,length(Sites))
      }
    }) %>% bind_cols()
    colnames(Interaction) <- paste0('InteractionZ_',c('Air','DayLen','Air.DayLen','Air.DayLen.Int'))
    Interaction <- Interaction %>% select(InteractionZ_Air.DayLen.Int)
    
    TrendLoad <- lapply(SS[[x]],FUN=function(yy){
      #yy<-Subset1_DFA[[x]][[2]]
      TLoad<-yy$Estimates$Z
      TLoad_out <- lapply(1:ncol(TLoad),FUN=function(zz){return(data.frame(TLoad[,zz]))})
      return(TLoad_out)
    }) %>% bind_cols()
    colnames(TrendLoad) <- paste0('TrendLoad_',c('Air','DayLen','Air.DayLen','Air.DayLen.Int'))
    
    R2 <- lapply(1:length(SS[[x]]),FUN=function(yy){
      YY<-SS[[x]][[yy]]
      pdf(paste0('DFAfits_',x+2010,'_',names(SS[[x]])[yy],'_Subset',sub.ind,'.pdf'))
      R2yy<-lapply(1:nrow(YY$ydat),FUN=function(ii){
        plot(as.numeric(colnames(YY$ydat)),YY$ydat[ii,],pch=16,ylab='Z-score water temp',xlab='DOY',type='o',main=row.names(YY$ydat)[ii])
        points(as.numeric(colnames(YY$ydat)),YY$Fits[ii,],type='l',col='blue',lwd=3)
        r2ii <- summary(lm(YY$ydat[ii,]~YY$Fits[ii,]))$r.squared
        return(r2ii)
      }) %>% unlist()
      dev.off()
      
      return(R2yy)
    }) %>% unlist()
    
    return(data.frame(SiteID = Sites,Year,TempSens,DayLenSens,Interaction,TrendLoad,R2))
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
    colnames(Trend) <- paste0('Trend_',c('Air','DayLen','Air.DayLen'))
    return(data.frame(Year=Year,DOY=DOY,Trend))
  }) %>% bind_rows()
  
  assign(paste0('GlobalResults_DayLen_Subset',sub.ind),value=GlobalResults,envir=.GlobalEnv)
  assign(paste0('GlobalTrends_DayLen_Subset',sub.ind),value=GlobalTrends,envir=.GlobalEnv)
  return(paste0('Subset ',sub.ind,' has written to Global environment'))
})

#Code to convert from Zscore to degree C space
lapply(1:3,FUN=function(sub.ind){
  #Set subset data.frame
  SS <- switch(sub.ind,GlobalResults_DayLen_Subset1,GlobalResults_DayLen_Subset2,GlobalResults_DayLen_Subset3)
  
  ZscoreConv <- lapply(1:nrow(SS),FUN=function(i){
    #i<-1
    #SSin <- SS
    datin <- Tdat %>% filter(SiteID == SS$SiteID[i], Year==SS$Year[i])
    aT<-datin %>% pull(airDT)
    sT<-datin %>% pull(meanDT)
    dayl <- datin %>% pull(global_dayl)
  
    ZSO<-data.frame(TempSens_Air = SS$TempSensZ_Air[i] * sd(sT,na.rm=T) * (1/sd(aT,na.rm=T)),
    TempSens_Air.DayLen = SS$TempSensZ_Air.DayLen[i] * sd(sT,na.rm=T) * (1/sd(aT,na.rm=T)),
    DayLenSens_DayLen = SS$DayLenSensZ_DayLen[i] * sd(sT,na.rm=T) * (1/sd(dayl,na.rm=T)),
    DayLenSens_Air.DayLen = SS$DayLenSensZ_Air.DayLen[i] * sd(sT,na.rm=T) * (1/sd(dayl,na.rm=T)))
    return(ZSO)
  }) %>% bind_rows()
  
  assign(paste0('GlobalResults_DayLen_Subset',sub.ind),value=SS %>% bind_cols(ZscoreConv),envir=.GlobalEnv)
  return('Zscores converted')
  })

#TempSens for interaction model
lapply(1:3,FUN=function(sub.ind){
  #Set subset data.frame
  SS <- switch(sub.ind,GlobalResults_DayLen_Subset1,GlobalResults_DayLen_Subset2,GlobalResults_DayLen_Subset3)
  
  ZscoreConv <- lapply(1:nrow(SS),FUN=function(i){
    #i<-1
    datin <- Tdat %>% filter(SiteID == SS$SiteID[i], Year==SS$Year[i])
    DOY <- datin$DOY
    aT<-datin %>% pull(airDT)
    sT<-datin %>% pull(meanDT)
    dayl <- datin %>% pull(daylen_min)
    TempSens_Air.DayLen.Int_byDOY <- (SS$TempSensZ_Air.DayLen.Int[i] + SS$InteractionZ_Air.DayLen.Int[i]*Zscore(dayl) ) * (sd(sT,na.rm=T) * (1/sd(aT,na.rm=T)))
    return(data.frame(Year=rep(SS$Year[i],length(DOY)),SiteID=rep(SS$SiteID[i],length(DOY)),DOY,TempSens_Int=TempSens_Air.DayLen.Int_byDOY))
  }) %>% bind_rows()
  
  assign(paste0('TempSens_Interaction_Subset',sub.ind),value=ZscoreConv,envir=.GlobalEnv)
  return('Zscores converted')
})

TempSens_Interaction_Subset1

save(list=ls(),file = 'DFAresults_GlobalDayLength_wInteraction_AllSubsets.Rdata')


# plot(GlobalTrends_DayLen_Subset1$Trend_Air,type='l')
# points(GlobalTrends_DayLen_Subset1$Trend_DayLen,type='l',col='red')
# points(GlobalTrends_DayLen_Subset1$Trend_Air.DayLen,type='l',col='blue')
# 
# head(GlobalResults_DayLen_Subset1)
# 
# plot(GlobalResults_DayLen_Subset1$TempSens_Air,GlobalResults_DayLen_Subset1$TempSens_Air.DayLen)
# abline(0,1)
# plot(GlobalResults_DayLen_Subset1$DayLenSens_DayLen,GlobalResults_DayLen_Subset1$DayLenSens_Air.DayLen)
# abline(0,1)
# 
# plot(GlobalResults_DayLen_Subset1$TempSens_Air.DayLen,GlobalResults_DayLen_Subset1$DayLenSens_Air.DayLen)
# abline(0,1)
# 
# hist(GlobalResults_DayLen_Subset1$TempSens_Air)
# hist(GlobalResults_DayLen_Subset1$TempSens_Air.DayLen,add=T,col='lightblue')
# 
# 
# 
# 
# 
# lapply(1:3,FUN=function(sub.ind){
#   GlobalResults <- switch(sub.ind,GlobalResults_Subset1,GlobalResults_Subset2,GlobalResults_Subset3)
#   GlobalTrends <- switch(sub.ind,GlobalTrends_Subset1,GlobalTrends_Subset2,GlobalTrends_Subset3)
#   Ntrends <- 3
#   
#   pdf(paste0('DFAresults_Subset',sub.ind,'.pdf'),width=8.5,height=11)
#   
#   lapply(1:Ntrends,FUN=function(TRD){
#     TLs <- switch(TRD,GlobalResults %>% select(TrendLoad_1_1),GlobalResults %>% select(TrendLoad_2_1,TrendLoad_2_2),GlobalResults %>% select(TrendLoad_3_1,TrendLoad_3_2,TrendLoad_3_3))
#     Ts <- switch(TRD,GlobalTrends %>% select(Trend_1_1),GlobalTrends %>% select(Trend_2_1,Trend_2_2),GlobalTrends %>% select(Trend_3_1,Trend_3_2,Trend_3_3))
#     
#     UniYear <- sort(unique(GlobalResults$Year))
#     Nyear <- length(UniYear)
#      lapply(1:Nyear,FUN=function(yr){
#        #yr<-1
#        TL_y <- TLs %>% filter(GlobalResults$Year == UniYear[yr])
#        T_y <- Ts %>% filter(GlobalTrends$Year == UniYear[yr])
#        DOY_y <- GlobalTrends %>% filter(GlobalTrends$Year == UniYear[yr]) %>% pull(DOY)
#        ColChoice <- c('#000000','#0000FF','#FF0000')
#        if(yr <= 5){pyr <- yr}else{pyr <- yr - 5}
#        
#        if(yr %in% c(1,6)){
#          par(fig=c(0,0.6,1-0.2*pyr,1-0.2*(pyr-1)))
#        }else{
#          par(fig=c(0,0.6,1-0.2*pyr,1-0.2*(pyr-1)),new=T)
#        }
#        par(mar=c(1,1,1,1),oma=c(3,3,0,0))
#        #Plot Trends
#        plot(DOY_y,T_y[,1],type='l',col=ColChoice[1],lwd=2,ylim=c(-3,3),xlab='DOY',ylab='Trend')
#        if(TRD>1){
#          for(TT in 2:ncol(T_y)){
#            points(DOY_y,T_y[,TT],type='l',col=ColChoice[TT],lwd=2)
#          }
#        }
#        abline(h=0,lty=3)
#        mtext(UniYear[yr],3,line=-1.5,cex=1)
#        
#        #Plot loadings histograms
#        par(fig=c(0.6,1,1-0.2*pyr,1-0.2*(pyr-1)),new=T)
#        plot(density(TL_y[,1]),type='l',col=ColChoice[1],lwd=2,xlim=c(-3,3),main='')
#        polygon(density(TL_y[,1]), col=paste0(ColChoice[1],'55'), border=ColChoice[1])
#        if(TRD>1){
#          for(TT in 2:ncol(T_y)){
#            points(density(TL_y[,TT]),type='l',col=ColChoice[TT],lwd=2)
#            polygon(density(TL_y[,TT]), col=paste0(ColChoice[TT],'55'), border=ColChoice[TT])
#          }
#        }
#        abline(v=0,lty=3)
#        
#        
#      }) #End Plotting over Year
#      
#      
#      })#End Apply over trends
#   dev.off() #Turnoff PDF
# }) # End apply over subsets
# 
# 
# lapply(1:3,FUN=function(sub.ind){
#   GlobalResults <- switch(sub.ind,GlobalResults_Subset1,GlobalResults_Subset2,GlobalResults_Subset3)  GlobalTrends <- switch(sub.ind,GlobalTrends_Subset1,GlobalTrends_Subset2,GlobalTrends_Subset3)
#   Ntrends <- 3
#   
#   pdf(paste0('DFAresults_Subset',sub.ind,'_withTempSens.pdf'),width=8.5,height=11)
#   
#   lapply(1:Ntrends,FUN=function(TRD){
#     TLs <- switch(TRD,GlobalResults %>% select(TrendLoad_1_1),GlobalResults %>% select(TrendLoad_2_1,TrendLoad_2_2),GlobalResults %>% select(TrendLoad_3_1,TrendLoad_3_2,TrendLoad_3_3))
#     Ts <- switch(TRD,GlobalTrends %>% select(Trend_1_1),GlobalTrends %>% select(Trend_2_1,Trend_2_2),GlobalTrends %>% select(Trend_3_1,Trend_3_2,Trend_3_3))
#     
#     UniYear <- sort(unique(GlobalResults$Year))
#     Nyear <- length(UniYear)
#     lapply(1:Nyear,FUN=function(yr){
#       #yr<-1
#       TL_y <- TLs %>% filter(GlobalResults$Year == UniYear[yr])
#       T_y <- Ts %>% filter(GlobalTrends$Year == UniYear[yr])
#       DOY_y <- GlobalTrends %>% filter(GlobalTrends$Year == UniYear[yr]) %>% pull(DOY)
#       ColChoice <- c('#000000','#0000FF','#FF0000')
#       if(yr <= 5){pyr <- yr}else{pyr <- yr - 5}
#       
#       if(yr %in% c(1,6)){
#         par(fig=c(0,0.6,1-0.2*pyr,1-0.2*(pyr-1)))
#       }else{
#         par(fig=c(0,0.6,1-0.2*pyr,1-0.2*(pyr-1)),new=T)
#       }
#       par(mar=c(1,1,1,1),oma=c(3,3,0,0))
#       #Plot Trends
#       plot(DOY_y,T_y[,1],type='l',col=ColChoice[1],lwd=2,ylim=c(-3,3),xlab='DOY',ylab='Trend')
#       if(TRD>1){
#         for(TT in 2:ncol(T_y)){
#           points(DOY_y,T_y[,TT],type='l',col=ColChoice[TT],lwd=2)
#         }
#       }
#       abline(h=0,lty=3)
#       mtext(UniYear[yr],3,line=-1.5,cex=1)
#       
#       #Plot loadings histograms
#       par(fig=c(0.6,1,1-0.2*pyr,1-0.2*(pyr-1)),new=T)
#       plot(density(TL_y[,1]),type='l',col=ColChoice[1],lwd=2,xlim=c(-3,3),main='')
#       polygon(density(TL_y[,1]), col=paste0(ColChoice[1],'55'), border=ColChoice[1])
#       if(TRD>1){
#         for(TT in 2:ncol(T_y)){
#           points(density(TL_y[,TT]),type='l',col=ColChoice[TT],lwd=2)
#           polygon(density(TL_y[,TT]), col=paste0(ColChoice[TT],'55'), border=ColChoice[TT])
#         }
#       }
#       abline(v=0,lty=3)
#       
#       
#     }) #End Plotting over Year
#     
#     
#   })#End Apply over trends
#   dev.off() #Turnoff PDF
# }) # End apply over subsets

#saveRDS(SiteSpecificAICtable,'SiteSpecificAICtable.RDS')
#saveRDS(GlobalAICTable,'GlobalAICtable.RDS')


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

