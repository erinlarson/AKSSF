
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


# save_daily_file <- function(data.in, acronym) {
#   daily_fields <- c("SiteID", "SampleDate", "meanDT", "minDT", "maxDT")
#   date.now <- as.Date(Sys.time())
#   data.in %>% 
#     select(one_of(daily_fields)) %>% 
#     write.csv(paste0("data_preparation/final_data/", acronym, "Daily_Data", as.character(date.now),".csv"), row.names = F)
# }
# 
# save_aktemp_file <- function(data.in, acronym) {
#   aktemp_fields <- c("SiteID", "sampleDate", "sampleTime", "Temperature", "UseData")
#   date.now <- as.Date(Sys.time())
#   data.in %>% 
#     select(one_of(aktemp_fields)) %>% 
#     write.csv(paste0("data_preparation/final_data/", acronym, "AKTemp_Data", as.character(date.now),".csv"), row.names = F)
# }
# 
# save_metadata_file <- function(data.in, acronym) {
#   date.now <- as.Date(Sys.time())
#   data.out <- paste0("data_preparation/final_data/", acronym, "Metadata", as.character(date.now),".csv")
#   akoats_fields <- c("SiteID", "seq_id", "Agency_ID", "SourceName", "Contact_person", 
#                      "Contact_email", "Contact_telephone", "Latitude", "Longitude", 
#                      "Sensor_Placement", "Waterbody_name", "Waterbody_type", "Sensor_accuracy", 
#                      "Sensor_QAQC")
#   data.in %>% 
#     select(one_of(akoats_fields)) %>% 
#     write.csv(data.out, row.names = F)
#   return(data.out)
# }

# Save csv locally and google drive copy of Metadata
save_metadata_files <- function(data.in, acronym) {
  # Date of file creation
  date.now <- as.Date(Sys.time())
  # Address to folder on google drive
  drive_path <- drive_get("https://drive.google.com/drive/u/0/folders/1_qtmORSAow1fIxvh116ZP7oKd_PC36L0")
  # Use same name for csv and googlesheet copy
  data.name <- paste0 ( acronym, "Metadata", as.character(date.now))
  # Location to save local csv
  csvname <- paste0 (data.name,".csv")
  # Path of csv - need to use a input for drive_upload
  csvpath <- paste0 ( rprojroot::find_rstudio_root_file(),"/data_preparation/final_data/", csvname)
  akoats_fields <- c ("SiteID", "seq_id", "Agency_ID", "SourceName", "Contact_person", 
                     "Contact_email", "Contact_telephone", "Latitude", "Longitude", 
                     "Sensor_Placement", "Waterbody_name", "Waterbody_type", "Sensor_accuracy", 
                     "Sensor_QAQC")
  data.in %>% 
    select(one_of(akoats_fields)) %>% 
    write.csv(csvpath, row.names = F)
  # Upload to google drive
  drive_upload(csvpath, path = as_id(drive_path), csvname ,overwrite = TRUE)
  print(paste0( "Metadata saved to ", csvpath))
  return(csvpath)
  
}

# Save csv locally and google drive copy of Daily Summaries
save_daily_files <- function(data.in, acronym) {
  # Date of file creation
  date.now <- as.Date(Sys.time())
  # Address to folder on google drive
  drive_path <- drive_get("https://drive.google.com/drive/u/0/folders/1_qtmORSAow1fIxvh116ZP7oKd_PC36L0")
  # Use same name for csv and googlesheet copy
  data.name <- paste0 ( acronym, "Daily_Data", as.character(date.now))
  # Location to save local csv
  csvname <- paste0 (data.name,".csv")
  # Path of csv - need to use a input for drive_upload
  csvpath <- paste0 ( rprojroot::find_rstudio_root_file(),"/data_preparation/final_data/", csvname)
  # Fields to save
  daily_fields <- c("SiteID", "sampleDate", "meanDT", "minDT", "maxDT")
  
  data.in %>% 
    select(one_of(daily_fields)) %>% 
    write.csv(csvpath, row.names = F)
  # Upload to google drive
  drive_upload(csvpath, path = as_id(drive_path), csvname,overwrite = TRUE)
  print(paste0( "Daily summaries saved to ", csvpath))
  return(csvpath)
  
}

# Save copy for AKTEMP
save_aktemp_files <- function(data.in, acronym) {
  # Date of file creation
  date.now <- as.Date(Sys.time())
  # Address to folder on google drive
  drive_path <- drive_get("https://drive.google.com/drive/u/0/folders/1_qtmORSAow1fIxvh116ZP7oKd_PC36L0")
  # Use same name for csv and googlesheet copy
  data.name <- paste0 ( acronym, "AKTEMP_Data", as.character(date.now))
  # Location to save local csv
  csvname <- paste0 (data.name,".csv")
  # Path of csv - need to use a input for drive_upload
  csvpath <- paste0 ( rprojroot::find_rstudio_root_file(),"/data_preparation/final_data/", csvname)
  
  aktemp_fields <- c("SiteID", "sampleDate", "sampleTime", "Temperature", "UseData")
  date.now <- as.Date(Sys.time())
  data.in %>% 
    select(one_of(aktemp_fields)) %>% 
    write.csv(csvpath, row.names = F)
  # Upload to google drive
  drive_upload(csvpath, path = as_id(drive_path), csvname ,overwrite = TRUE)
  print(paste0( "AKTEMP copy saved to ", csvpath))
  return(csvpath)
}

# Read in csv and add column with filename
read_csv_and_name <- function(csv_file_path) {
  sheet_name <- str_match(csv_file_path, "\\/\\s*(.*?)\\s*\\.csv")[2]
  dat <- read_csv(csv_file_path) %>% 
    mutate(file_name = sheet_name)
}

# Function to pull logger serial# from temp column
log_sn_fun <- function(string){ 
  
  str_extract(string, "\\-*\\d+\\.*\\d*")
} 


