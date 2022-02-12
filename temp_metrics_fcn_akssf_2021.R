
# Temp Metrics
#####################
# 40 metrics from previous Mat-Su work
# updates based on original R script: Temperature Descriptor Function - daily inputs.R
####################

#Notes:
# - convert to kelvin for correct calculation of cv
# - added sd since units are in c vs c2 for var


# datainput <- sites_complete_summers


tempmetrics <- function(datainput) {
  
  #look for and install needed packages
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(tidyverse, zoo, lubridate)
  
  tab <- data.frame()
  
  site_years <- datainput %>% 
    distinct(SiteID, Year) %>% 
    ungroup()
  
  for(i in 1:nrow(site_years)) {
    dat <- left_join(site_years %>% slice(i), datainput) %>% 
      arrange(sampleDate) 
    
    #magnitude metrics
    mag_mets <- dat %>%
      mutate(ma_mean = rollapply(meanDT, 7, mean, align = 'center', fill = NA),
             ma_max = rollapply(maxDT, 7, mean, align = 'center', fill = NA)) %>% 
      summarize(MA7d_DAT = max(ma_mean, na.rm = TRUE),
                MA7d_DMT = max(ma_max, na.rm = TRUE),
                MxDAT = max(meanDT, na.rm = TRUE),
                MnDAT = mean(meanDT, na.rm = TRUE),
                MxDMT = max(maxDT, na.rm = TRUE),
                MnDMT = mean(maxDT, na.rm = TRUE),
                MnDNT = mean(minDT, na.rm = TRUE),
                NDNT = min(minDT, na.rm = TRUE)) %>% 
      mutate(MxDIFF = MxDMT - NDNT)
    
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

    mag <- bind_cols(mag_mets, monthly)
    
    #variability metrics
    var <- dat %>% 
      mutate(meanDT_Kelvin = meanDT + 273.15,
             drange = maxDT - minDT) %>% 
      summarize(SIGMA_DAT = var(meanDT),
                SIGMA_DMT = var(maxDT),
                SIGMA_DNT = var(minDT),
                SD = sd(meanDT),
                CV_DAT = sd(meanDT)/mean(meanDT), #for comparison only
                CV_DAT_K = sd(meanDT_Kelvin)/mean(meanDT_Kelvin),
                RANGE_MAX = max(drange, na.rm = TRUE),
                RANGE_MN = mean(drange, na.rm = TRUE)) 
    

    #frequency metrics
    freq <- dat %>% 
      summarize(SUM_13_DMT = sum(maxDT > 13),
                SUM_18_DMT = sum(maxDT > 18),
                SUM_20_DMT = sum(maxDT > 20))
    
    #duration metrics
    #these are mean event length and max event length greater than thresholds
    #added also original duration metrics: no. of events longer than 6 days and number of events total
    dur13 <- rle(dat$maxDT > 13) %>%
      unclass() %>% 
      as_tibble() %>% 
      summarize(dur13 = case_when(sum(values == TRUE) > 0 ~ mean(lengths[values == TRUE]),
                                  TRUE ~ 0),
                dur13mx = case_when(sum(values == TRUE) > 0 ~ as.numeric(max(lengths[values == TRUE])),
                                    TRUE ~ 0),
                dur13gt6 = sum(values[values == TRUE & lengths > 6]),
                dur13ct = sum(values[values == TRUE]))
    dur18 <- rle(dat$maxDT > 18) %>%
      unclass() %>% 
      as_tibble() %>% 
      summarize(dur18 = case_when(sum(values == TRUE) > 0 ~ mean(lengths[values == TRUE]),
                                  TRUE ~ 0),
                dur18mx = case_when(sum(values == TRUE) > 0 ~ as.numeric(max(lengths[values == TRUE])),
                                    TRUE ~ 0),
                dur18gt6 = sum(values[values == TRUE & lengths > 6]),
                dur18ct = sum(values[values == TRUE]))
    dur20 <- rle(dat$maxDT > 20) %>%
      unclass() %>% 
      as_tibble() %>% 
      summarize(dur20 = case_when(sum(values == TRUE) > 0 ~ mean(lengths[values == TRUE]),
                                  TRUE ~ 0),
                dur20mx = case_when(sum(values == TRUE) > 0 ~ as.numeric(max(lengths[values == TRUE])),
                                    TRUE ~ 0),
                dur20gt6 = sum(values[values == TRUE & lengths > 6]),
                dur20ct = sum(values[values == TRUE]))
    
    dur <- cbind(dur13, dur18, dur20)
    
    #timing metrics
    timmax <- dat %>% 
      filter(maxDT == mag$MxDMT)%>% 
      slice(1) %>% #first day when the daily maximum occurs 
      pull(sampleDate) %>% 
      yday()
    timmwmt <- dat %>%
      mutate(ma_max = rollapply(maxDT, 7, mean, align = 'center', fill = NA)) %>% 
      filter(ma_max == mag$MA7d_DMT) %>% 
      slice(1) %>% 
      pull(sampleDate) %>% 
      yday()
    tim <- data.frame(MxDMT_jd = timmax, MA7d_DMT_jd = timmwmt)
    
    newrow <- bind_cols(site_years %>% slice(i), mag, var, freq, dur, tim)
    tab <- bind_rows(tab, newrow) 
  }
  return(tab)
}