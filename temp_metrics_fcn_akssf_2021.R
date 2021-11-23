
# Temp Metrics
#####################
# 10 metrics from previous Mat-Su work
# updates based on original R script: Temperature Descriptor Function - daily inputs.R
####################

#Notes:
# - convert to kelvin for correct calculation of cv
# - added sd since units are in c vs c2 for var
# - converted several metrics based on daily maximum temperatures to daily mean 
#   since we only have means for UW dataset. Also removed anything related to
#   daily ranges.

# datainput <- deshka_sites_complete_summers

tempmetrics <- function(datainput) {
  
  #look for and install needed packages
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(tidyverse, zoo, lubridate)
  
  tab <- data.frame()
  
  site_years <- datainput %>% 
    distinct(SiteID, year) %>% 
    ungroup()
  
  for(i in 1:nrow(site_years)) {
    dat <- left_join(site_years %>% slice(i), datainput) %>% 
      arrange(sampleDate) 
    
    #magnitude metrics
    ma <- dat %>%
      mutate(ma_mean = rollapply(meanDT, 7, mean, align = 'center', fill = NA)) %>% 
      summarize(MA7d_DAT = max(ma_mean, na.rm = TRUE),
                MxDAT = max(meanDT, na.rm = TRUE)) %>% 
      select(MA7d_DAT, MxDAT)
    
    monthly <- dat %>% 
      mutate(month_abb = month(sampleDate, label = TRUE, abbr = TRUE),
             month = month(sampleDate),
             month_den = days_in_month(sampleDate)) %>% 
      filter(month %in% 5:9) %>% 
      select(-month) %>% 
      group_by(month_abb, month_den) %>% 
      summarize(mon_mn = mean(meanDT, na.rm = TRUE),
                mon_ct = n()) %>%
      filter(mon_ct/month_den > 0.8) %>% 
      select(-month_den, -mon_ct) %>% 
      pivot_wider(values_from = mon_mn, names_from = month_abb)
    
    mag <- bind_cols(ma, monthly)
    
    #variability metrics
    var <- dat %>% 
      mutate(meanDT_Kelvin = meanDT + 273.15) %>% 
      summarize(SIGMA_DAT = var(meanDT),
                SD = sd(meanDT),
                CV_DAT = sd(meanDT)/mean(meanDT), #for comparison only
                CV_DAT_K = sd(meanDT_Kelvin)/mean(meanDT_Kelvin)) %>% 
      select(SIGMA_DAT, SD, CV_DAT, CV_DAT_K)
    
    #frequency metrics
    freq <- dat %>% 
      summarize(SUM_13 = sum(meanDT > 13),
                SUM_18 = sum(meanDT > 18),
                SUM_20 = sum(meanDT > 20))
    
    #duration metrics
    rle13 <- data.frame(unclass(rle(dat$meanDT > 13)))
    dur13 <- ifelse(nrow(rle13[rle13$values==TRUE,]) > 0, mean(rle13[rle13$values==TRUE,"lengths"]), 0)
    rle18 <- data.frame(unclass(rle(dat$meanDT > 18)))
    dur18 <- ifelse(nrow(rle18[rle18$values==TRUE,]) > 0, mean(rle18[rle18$values==TRUE,"lengths"]), 0)
    rle20 <- data.frame(unclass(rle(dat$meanDT > 20)))
    dur20 <- ifelse(nrow(rle20[rle20$values==TRUE,]) > 0, mean(rle20[rle20$values==TRUE,"lengths"]), 0)
    dur <- data.frame(DUR_mn13 = dur13, DUR_mn18 = dur18, DUR_mn20 = dur20)
    
    #timing metrics
    timmax <- dat %>% 
      filter(meanDT == mag$MxDAT)%>% 
      slice(1) %>% 
      pull(sampleDate) %>% 
      yday()
    timmwmt <- dat %>%
      mutate(ma_mean = rollapply(meanDT, 7, mean, align = 'center', fill = NA)) %>% 
      filter(ma_mean == mag$MA7d_DAT) %>% 
      slice(1) %>% 
      pull(sampleDate) %>% 
      yday()
    tim <- data.frame(MxDAT_jd = timmax, MA7d_DAT_jd = timmwmt)
    
    newrow <- bind_cols(site_years %>% slice(i), mag, var, freq, dur, tim)
    tab <- bind_rows(tab, newrow) 
  }
  return(tab)
}