# AKSSF
AKSSF project and workflow

Project to describe and map thermal sensitivities for Southcentral Alaska - Cook Inlet, Prince William Sound, and Copper River.

(Note: supplemental funding from FWS to include Kodiak temperature data in this project. And, second AKSSF award to do a similar analysis for Bristol Bay.)

1. Data Availability

* Deshka temperature model output - CIK and FWS
* Kenai temperature model output - Ben Meyers, CIK and KWF
* APU data for Anchor R from John Hagan - cleaned for anchor temperature model
* KBRR data from Lower Kenai from Steve Baird
* CIK data for Anchor from Sue Mauger for temperature modeling project
* CIK data on KNB
* USFS data for streams in Chugach National Forest from Luca Adelfio
* ACCS temperature data from Little Susitna watershed - summer 2016 and 2019-2020
* USGS sites should be selected from dataRetrieval library in R
* FWS data for Kodiak Refuge and other sites from Meg Perdue - this is post data archive on KNB
* FWS data on KNB from Meg Perdue and Don Rivard
* ADFG data on KNB for Kodiak Island - but maybe only a couple sites have continuous data
* NPS data on KNB from Trey Simmons
* ARRI data from Mat-Su - dailies are an output from thermal regimes project

All Bristol Bay data are being cleaned in a separate [gihub repo](https://github.com/rsshaftel/SWSHP-Bristol-Bay-Thermal-Diversity).



2. Data cleaning

All data need to be formatted consistently so that we can combine them for the data analysis step. All datasets should be saved as separate data files and sites files for import to AKTEMP later. 

Data file

* SiteID
* sampleDate
* sampleTime
* Temeperature
* useData

Sites file

* SiteID
* AKOATS_ID
* latitude
* longitude
* Source_Name
* Contact_Name

3. Air temperature extraction

The sites file will be used in ArcGIS to extract catchments for each site. In R, we can calculate the 3-day moving average of air temperatures from DAYMET data for the appropriate months and years of data. We aren't predicting with air temperatures so it only needs to match the empirical data. There is code in the temperature modeling repos to get this started.

4. Thermal sensitivity analysis

Per Daniel, Tim Cline has the correct code for running the DFA. He did a slightly different procedure than what Peter Lisi did in the original paper for Bristol Bay. See Cline's snow paper and DFA.


