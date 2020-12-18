# AKSSF
AKSSF project and workflow

Project to describe and map thermal sensitivities for Southcentral Alaska - Cook Inlet, Prince William Sound, and Copper River.

(Note: supplemental funding from FWS to include Kodiak temperature data in this project. And, second AKSSF award to do a similar analysis for Bristol Bay.)

1. Data Availability

Metadata and data ready:

* Deshka temperature model output - CIK and FWS data from 2017-2019
* Kenai temperature model output - Ben Meyers, CIK and KWF
* Anchor temperature model output - CIK, APU, and KBRR
* USFS data for streams in Chugach National Forest from Luca Adelfio

Metadata ready, but not data:

* CIK data on KNB - this will include all 48 sites in network, but only through 2017.
* FWS data for Kodiak Refuge, OSM sites on Kodiak (2) and Copper River (2), and WRB sites on Kodiak (2) -- from Meg Perdue, received Oct 2020.
* NPS data on KNB from Trey Simmons
* USGS sites selected in GIS from akoats working. Get data from dataRetrieval library in R
* ACCS temperature data from Little Susitna watershed - summer 2016 and 2019-2020. Dustin is doing data prep in a separate [github repo](https://github.com/rsshaftel/Little_Susitna). Grabbed 2016 sites from akoats metadata and combined with 2019 sites from logger database, some might be the same, will need to check in GIS.

Low priority datasets (wait to see if we decide to only look at sensitivity for sites with a minimum number of years of data before we pursue additional datasets):

* ADFG data on KNB for Kodiak Island - but maybe only a couple sites have continuous data
* Tribal data in AKOATS for Kodiak Island - needs to be requested
* Mat-Su thermal regimes project - daily data were saved as a final product. This would include CIK data already archived on KNB, USGS data that I can grab separately, and ARRI data from Mat-Su. I could just re-import the raw ARRI data, I believe they cleaned it for me.
* AEA data for Susitna-Watana dam. I looked at this data for the thermal regimes project and it wasn't clean. They had temperature arrays on the Susitna River with some pretty anomalous readings between loggers. We could decide to clean this dataset for this project, but it would probably be a week of someone's time. 


Bristol Bay:

All Bristol Bay data are being cleaned in a separate [github repo](https://github.com/rsshaftel/SWSHP-Bristol-Bay-Thermal-Diversity).
(Note that the only data we are missing in Bristol Bay are Dan Young's sites, otherwise everything should be good.)



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


