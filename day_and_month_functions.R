
# 2/11/21, Rebecca Shaftel

# Function to assign the mode of the frequency of temperature measurements for each day. This value can be used to determine
# if there are adequate daily measurements to calculate daily summaries. 

# Requirements:
# Input data frame needs variables named SiteID (character), sampleDate (date), and sampleTime (hms).
# Output data frame will have new variables called time_diff (difference in minutes between measurement and most recent measurement)
# and mode_diff (most common time_diff for a site and sample date).

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

test <- c(NA, 1, 1, NA)

Mode(test)


temp_msmt_freq <- function(data.in) {
  dat.out <- data.in %>% 
    group_by(SiteID, sampleDate) %>% 
    arrange(SiteID, sampleDate, sampleTime) %>% 
    mutate(time_diff = as.numeric(sampleTime - lag(sampleTime), units = 'mins'),
           mode_diff = Mode(time_diff)) 
  return(dat.out)
}

daily_screen <- function(data.in) {
  dat.out <- data.in %>% 
    filter(UseData == 1, !is.na(mode_diff)) %>% 
    group_by(SiteID, sampleDate, mode_diff) %>%
    summarize(meanDT = mean(Temperature, na.rm = TRUE),
              minDT = min(Temperature, na.rm = TRUE),
              maxDT = max(Temperature, na.rm = TRUE),
              msmtCt = n()) %>% 
    filter(msmtCt > (0.9 * 1440/mode_diff))
}


numberOfDays <- function(date) {
  m <- format(date, format="%m")
  while (format(date, format="%m") == m) {
    date <- date + 1
  }
  return(as.integer(format(date - 1, format="%d")))
}