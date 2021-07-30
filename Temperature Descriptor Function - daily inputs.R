

library(caTools); library(lubridate); library(XLConnect)

# rm(list=ls(all=TRUE))
# library(xlsReadWrite)
# library(zoo)
# library(xts)
###############################
# Temperature Descriptor Function #
###############################

#Create separate function that screens data and outputs final data frame used for metrics

tempscreen <- function(datainput) {
  library(caTools); library(lubridate); library(XLConnect)
  #Add in screening for a complete year, 90% of days, before calculating metrics.
  #Start with JJA. Create a vector of site-years that meet 90% cutoff.

  datainput$site.year<-paste(datainput$site_name,year(datainput$date),sep=".")
  datainput <- na.omit(datainput) #need to put this here bc some usgs data still problematic
  mo.ct<-as.data.frame(table(list(month(datainput$date),datainput$site.year)))
  colnames(mo.ct)<-c("month_no","site.year","count")
  mo.ct<-reshape(mo.ct,v.names="count",direction="wide",timevar="month_no",idvar="site.year")
  mo.ct$JJA<-rowSums(mo.ct[,names(mo.ct) %in% c("count.6","count.7","count.8")])
  sy.sum<-as.character(mo.ct[mo.ct$JJA>82,"site.year"]) 
  
  #get summer data for complete site-years
  datainput.mod<-datainput[datainput$site.year %in% sy.sum & month(datainput$date) %in% 6:8,]  
  
  #  mo.ct$JA<-rowSums(mo.ct[,names(mo.ct) %in% c("count.7","count.8")])
  #  mo.ct[mo.ct$JA>55,"site.year"] 
  
  return(datainput.mod)
}

tempscreen.JA <- function(datainput) {
  library(caTools); library(lubridate); library(XLConnect)
  #Add in screening for a complete year, 90% of days, before calculating metrics.
  #New function that outputs just July/August data. Create a vector of site-years that meet 90% cutoff.
  
  datainput$site.year<-paste(datainput$site_name,year(datainput$date),sep=".")
  datainput <- na.omit(datainput) #need to put this here bc some usgs data still problematic
  mo.ct<-as.data.frame(table(list(month(datainput$date),datainput$site.year)))
  colnames(mo.ct)<-c("month_no","site.year","count")
  mo.ct<-reshape(mo.ct,v.names="count",direction="wide",timevar="month_no",idvar="site.year")
  mo.ct$JA<-rowSums(mo.ct[,names(mo.ct) %in% c("count.7","count.8")])
  sy.sum<-as.character(mo.ct[mo.ct$JA>55,"site.year"]) 
  
  #get summer data for complete site-years
  datainput.mod<-datainput[datainput$site.year %in% sy.sum & month(datainput$date) %in% 7:8,]  
  
  return(datainput.mod)
}




tempmetrics <- function(datainput,outtable) {

  library(caTools); library(lubridate); library(XLConnect)
  
  tab <- data.frame(site.year=as.character(),MxDAT=as.numeric(),MnDAT=as.numeric(),MxDMT=as.numeric(),MnDMT=as.numeric(),
                          MnDNT=as.numeric(),NDNT=as.numeric(),MxWMT=as.numeric(),MxWAT=as.numeric(),CDD=as.numeric(),WAT5=as.numeric(),
                          WAT25=as.numeric(),WAT50=as.numeric(),WAT75=as.numeric(),WAT95=as.numeric(),MA7d_DMT=as.numeric(),
                          MA7d_DAT=as.numeric(),DELTA_MAX=as.numeric(),SIGMA_MN=as.numeric(),SIGMA_MIN=as.numeric(),SIGMA_MAX=as.numeric(),
                          CV_MN=as.numeric(),CV_MIN=as.numeric(),CV_MAX=as.numeric(),RNG=as.numeric(),SUM_13=as.numeric(),SUM_18=as.numeric(),
                          SUM_20=as.numeric(),SUM7d_13=as.numeric(),SUM7d_18=as.numeric(),SUM_event13=as.numeric(),SUM_event18=as.numeric(),
                          SUM_MA7d_13=as.numeric(),SUM_MA7d_18=as.numeric(),SUM_MA7d_20=as.numeric(),
                          DUR_mx13=as.numeric(),DUR_mn13=as.numeric(),DUR_mx18=as.numeric(),DUR_mn18=as.numeric(),
                          DUR_mx20=as.numeric(),DUR_mn20=as.numeric(),MxDMT_jd=as.numeric(),MA7d_DMT_jd=as.numeric(),
                          WDMT_jd=as.numeric(),D_RNG=as.numeric(),stringsAsFactors=FALSE)
  
#Setup looping to go over all site-years in data frame.  datainput organized as:
# site, date, max, min, mean.
  sy.sum <- unique(datainput$site.year)
  for (j in 1:length(sy.sum)) {
    in.yr <- datainput[datainput$site.year == sy.sum[j],]
    in.yr <- in.yr[month(in.yr$date) %in% 6:8,]
    tab[j,"site.year"]<-sy.sum[j]

    #Magnitude
    tab[j,"MxDAT"] <- max(in.yr$mean) 
    tab[j,"MnDAT"] <- mean(in.yr$mean) 
    tab[j,"MxDMT"] <- max(in.yr$max) 
    tab[j,"MnDMT"] <- mean(in.yr$max) 
    tab[j,"MnDNT"] <- mean(in.yr$min) 
    tab[j,"NDNT"] <- min(in.yr$min) 

    in.yr$week <- as.numeric(format(in.yr$date,"%U"))
    wk.ct <- with(in.yr,table(week))            #number of days in each week.
    wk <- with(in.yr,tapply(mean,week,mean))    #vector of weekly means.
    wk <- wk[wk.ct>5]                           #screen for weeks with 6 or more days.
    wk.mx <- with(in.yr,tapply(max,week,mean))  #vector of weekly means based on daily max.
    wk.mx <- wk.mx[wk.ct>5]
    tab[j,"MxWMT"] <- max(wk.mx)
    tab[j,"MxWAT"] <- max(wk)
  
    tab[j,"CDD"] <- sum(in.yr$mean) #CDD for JJA

    tab[j,"WAT5"] <- as.numeric(quantile(wk,0.05)) #5th pctile mean weekly temps
    tab[j,"WAT25"] <- as.numeric(quantile(wk,0.25)) #25th pctile mean weekly temps
    tab[j,"WAT50"] <- as.numeric(quantile(wk,0.50)) #50th pctile mean weekly temps
    tab[j,"WAT75"] <- as.numeric(quantile(wk,0.75)) #75th pctile mean weekly temps
    tab[j,"WAT95"] <- as.numeric(quantile(wk,0.95)) #95th pctile mean weekly temps
    
    MA_max <- runmean(in.yr$max,k=7,endrule="NA",align="center") #7-day moving average of daily max temps
    tab[j,"MA7d_DMT"] <- max(MA_max, na.rm=TRUE) 
    tab[j,"MA7d_DAT"] <- max(runmean(in.yr$mean,k=7,endrule="NA",align="center"),na.rm=TRUE)
    
    #Variability
    tab[j,"DELTA_MAX"] <- max(in.yr$max - in.yr$min) #Maximum daily temperature range
    tab[j,"SIGMA_MN"] <- var(in.yr$mean) #variance of daily means
    tab[j,"SIGMA_MIN"] <- var(in.yr$min) #variance of daily mins
    tab[j,"SIGMA_MAX"] <- var(in.yr$max) #variance of daily maxes
    tab[j,"CV_MN"] <- sd(in.yr$mean)/mean(in.yr$mean) #CV of daily means
    tab[j,"CV_MIN"] <- sd(in.yr$min)/mean(in.yr$min) #CV of daily mins
    tab[j,"CV_MAX"] <- sd(in.yr$max)/mean(in.yr$max) #CV of daily maxes
    
    tab[j,"RNG"] <- tab[j,"MxDMT"] - tab[j,"NDNT"] #Overall max minus min
    tab[j,"D_RNG"] <- mean(in.yr$max-in.yr$min)    #Mean of daily ranges 
    
    #Frequency
#     w_dm <- as.xts(daily_max(sumsub))
#     week_max <- split.xts(w_dm, f="weeks")
#     HOTWEEK <- as.numeric(which.max(apply.weekly(w_dm, FUN = mean)))
#     WEEK_14 <- as.numeric(sum(week_max[[HOTWEEK]]>14)) #Cumulative days >14C in hottest week
#     WEEK_18 <- as.numeric(sum(week_max[[HOTWEEK]]>18)) #Cumulative days >18C in hottest week
#     WEEK_22 <- as.numeric(sum(week_max[[HOTWEEK]]>22)) #Cumulative days >22C in hottest week
    tab[j,"SUM_13"] <- as.numeric(sum(in.yr$max > 13)) #Cumulative days >13C in summer
    tab[j,"SUM_18"] <- as.numeric(sum(in.yr$max > 18)) #Cumulative days >18C in summer
    tab[j,"SUM_20"] <- as.numeric(sum(in.yr$max > 20)) #Cumulative days >20C in summer
    
    rle13 <- as.data.frame(unclass(rle(in.yr$max>13))) #dataframe with lengths of events above 13C
    rle18 <- as.data.frame(unclass(rle(in.yr$max>18))) #dataframe with lengths of events above 18C
    rle20 <- as.data.frame(unclass(rle(in.yr$max>20))) #dataframe with lengths of events above 20C
    
    tab[j,"SUM7d_13"] <- dim(rle13[rle13$values==TRUE & rle13$lengths>6,]) [1]
    tab[j,"SUM7d_18"] <- dim(rle18[rle18$values==TRUE & rle18$lengths>6,]) [1]
    tab[j,"SUM_event13"] <- dim(rle13[rle13$values==TRUE,]) [1]
    tab[j,"SUM_event18"] <- dim(rle18[rle18$values==TRUE,]) [1]
    
    tab[j,"SUM_MA7d_13"] <- sum(MA_max > 13,na.rm=TRUE)
    tab[j,"SUM_MA7d_18"] <- sum(MA_max > 18,na.rm=TRUE)
    tab[j,"SUM_MA7d_20"] <- sum(MA_max > 20,na.rm=TRUE)
    
    #Duration - THESE ARE ALL FREQUENCIES, FIGURE OUT DURATION CODE OR USE FROM ARISMENDI
#     DMOV10 <- as.numeric(sum(rollmean(daily_max(sumsub), 3) > 10)) #No. days seven day rolling mean of daily maxes > 10C
#     DMOV15 <- as.numeric(sum(rollmean(daily_max(sumsub), 3) > 15)) #No. days seven day rolling mean of daily maxes > 15C
#     DMOV20 <- as.numeric(sum(rollmean(daily_max(sumsub), 3) > 20)) #No. days seven day rolling mean of daily maxes > 20C
#     DMAX10 <- as.numeric(sum(daily_max(sumsub) > 10)) #No. days daily max exceeded 10 C
#     DMAX15 <- as.numeric(sum(daily_max(sumsub) > 15)) #No. days daily max exceeded 15 C
#     DMAX20 <- as.numeric(sum(daily_max(sumsub) > 20)) #No. days daily max exceeded 20 C
#     PG5 <- (sum(sumsub > 5) / length (sumsub))*100 #% of temps > 5C
#     PG10 <- (sum(sumsub > 10) / length (sumsub))*100 #% of temps > 10C
#     PG15 <- (sum(sumsub > 15) / length (sumsub))*100 #% of temps > 15C
#     PG20 <- (sum(sumsub > 20) / length (sumsub))*100 #% of temps > 20C

    tab[j,"DUR_mx13"] <- ifelse(dim(rle13[rle13$values==TRUE,])[1]>0,max(rle13[rle13$values==TRUE,"lengths"]),0)
    tab[j,"DUR_mn13"] <- ifelse(dim(rle13[rle13$values==TRUE,])[1]>0,mean(rle13[rle13$values==TRUE,"lengths"]),0)
    tab[j,"DUR_mx18"] <- ifelse(dim(rle18[rle18$values==TRUE,])[1]>0,max(rle18[rle18$values==TRUE,"lengths"]),0)
    tab[j,"DUR_mn18"] <- ifelse(dim(rle18[rle18$values==TRUE,])[1]>0,mean(rle18[rle18$values==TRUE,"lengths"]),0)
    tab[j,"DUR_mx20"] <- ifelse(dim(rle20[rle20$values==TRUE,])[1]>0,max(rle20[rle20$values==TRUE,"lengths"]),0)
    tab[j,"DUR_mn20"] <- ifelse(dim(rle20[rle20$values==TRUE,])[1]>0,mean(rle20[rle20$values==TRUE,"lengths"]),0)
    
    
    #Timing
    tab[j,"MxDMT_jd"] <- as.POSIXlt(in.yr[in.yr$max==tab[j,"MxDMT"],"date"])$yday[1]
    tab[j,"MA7d_DMT_jd"] <- as.POSIXlt(in.yr[which(MA_max==tab[j,"MA7d_DMT"])[1],"date"])$yday
    tab[j,"WDMT_jd"] <- as.POSIXlt(in.yr[in.yr$week==as.numeric(names(which.max(wk.mx))),"date"][4])$yday

  }

# twb<-loadWorkbook(paste(outtable,".xlsx",sep=""),create=TRUE)
# createSheet(twb,name="Results")
# writeWorksheet(twb,tab,sheet="Results")
# saveWorkbook(twb)

return(tab)
}

