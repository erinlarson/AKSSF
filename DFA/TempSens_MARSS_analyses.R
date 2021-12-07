#####
##  Dynamic Factor Analyses for AKSSF
##  Timothy Cline
#####

setwd('/home/aq_user/aq_virtual_machine/AKSSF/DFA')

#Load required libraries
library(dplyr)
#library(TMB)
library(doParallel)
library(foreach)
library(MARSS)

#function for standardizing time series
Zscore <- function(x){return((x-mean(x,na.rm=T))/sd(x,na.rm=T))}

#load metadata for region id
MetaDat <- readRDS('../../data/md.rds')

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
TempSens_MARSS <- lapply(2011:2020,FUN=function(y){
  #y<-2011
  MeanDWT_TS <- Tdat %>% filter(Year==y) %>% xtabs(meanDT ~ SiteID + DOY,data=.)
  MeanDWT_TS[MeanDWT_TS==0] <- NA #xtabs puts zeros in for missing dates in the dataset, we set those exact zeros to NA
  
  #keep sites with >80% of data
  MeanDWT_TS <- MeanDWT_TS[which(rowSums(!is.na(MeanDWT_TS)) >= round(0.8*length(DOYseq)) ),]
  
  #build time series of matching air temp
  MeanDAT_TS <- Tdat %>% filter(Year==y, SiteID %in% row.names(MeanDWT_TS)) %>% xtabs(airDT ~ SiteID + DOY,data=.)
  
  TempSens <- foreach(s = 1:nrow(MeanDWT_TS),.packages=c('MARSS')) %dopar% {
    m1 <- MARSS(MeanDWT_TS[s,],model=list(R='zero',Q='diagonal and unequal',B='diagonal and unequal',A='zero',U='zero',C='unequal',c=matrix(MeanDAT_TS[s,],nrow=1)))
    df1 <- data.frame(Site=row.names(MeanDWT_TS)[s],Year=y,TempSens=coef(m1,type='matrix')$C,AR1=coef(m1,type='matrix')$B,Err=coef(m1,type='matrix')$Q)
    df1
  } %>% bind_rows()
  
  return(TempSens)
}) %>% bind_rows()

TempSens_ARIMA <- lapply(2011:2020,FUN=function(y){
  #y<-2011
  MeanDWT_TS <- Tdat %>% filter(Year==y) %>% xtabs(meanDT ~ SiteID + DOY,data=.)
  MeanDWT_TS[MeanDWT_TS==0] <- NA #xtabs puts zeros in for missing dates in the dataset, we set those exact zeros to NA
  
  #keep sites with >80% of data
  MeanDWT_TS <- MeanDWT_TS[which(rowSums(!is.na(MeanDWT_TS)) >= round(0.8*length(DOYseq)) ),]
  
  #build time series of matching air temp
  MeanDAT_TS <- Tdat %>% filter(Year==y, SiteID %in% row.names(MeanDWT_TS)) %>% xtabs(airDT ~ SiteID + DOY,data=.)
  
  TempSens <- foreach(s = 1:nrow(MeanDWT_TS)) %dopar% {
    m1 <- arima(x=MeanDWT_TS[s,],order=c(1,0,0),xreg=MeanDAT_TS[s,])
    df1 <- data.frame(Site=row.names(MeanDWT_TS)[s],Year=y,TempSens=m1$coef[3],AR1=m1$coef[1],Err=m1$sigma2)
    df1
  } %>% bind_rows()
  
  return(TempSens)
}) %>% bind_rows()



#End cluster
stopCluster(cl)

plot(TempSens_MARSS$TempSens ~ TempSens_ARIMA$TempSens)

library(stringr)

TodaysDate <- paste(str_split(as.Date(Sys.time()),pattern='-')[[1]],collapse='.')

saveRDS(TempSens_ARIMA,file=paste0('TempSens_Arima_Results',TodaysDate,'.RDS'))




