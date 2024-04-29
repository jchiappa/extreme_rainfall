# extreme_rainfall
  <p>
Detecting 12-hour extreme rainfall events from Stage IV data. Used to generate the database analyzed in Chiappa et al. (2024).

[insert reference]

  </p>
</div>

## Prerequisites

* At least 1 TB of storage space recommended
* Python 3 installed with dependencies listed in code


## Part I: Data Retrieval and Prep

### Stage IV

1. Download Stage IV data from EOL archive via FTP order https://doi.org/10.5065/D6PG1QDD
2. Extract all files that are not in .grb2 format (.Z or .gz format) and separate 01h, 06h, and 24h files into separate directories (conus only, when applicable)
3. Run `01_stageiv_renaming.py` -- Renames files to keep consistent naming convention and add grb2 extension where necessary (Note: Files before 20 July 2020 are technically GRIB1 format, but adding the grb2 extension does not impact the readability of the file.)
4. Optional: Run `02_check_missing.py` and see exported csv files for times of missing Stage IV data files

### NOAA Atlas 14

1. Download `NOAA_Atlas_14_CONUS.nc` from https://hdsc.nws.noaa.gov/pub/hdsc/data/tx/
2. Place in a designated folder

### States
1. If not downloaded, get the ArcGIS states shapefiles (by jeff.elrod) from https://www.arcgis.com/home/item.html?id=66d95ce91afb48829a9fe9c3145c755d
2. Extract and place folder in main data directory if not already there.

### Tropical Cyclones
1. Download the latest IBTrACS data in csv format (`ibtracs.since1980.list.v04r00.csv`) from https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r00/access/csv/
2. Place in main data directory


## Part II: Preliminary ERE Detection Algorithm

Run `03_ere_preliminary_algorithm.py` after changing paths and variables in the code (will take several hours to run -- around 1 hour per year of data)


## Part III: Quality Control
1. Run `04_qc1_params.py` to compute relevant quality control parameters for each preliminary event. Takes several minutes.
2. Run `05_qc1.py` to create and apply QC conditions. Will be applied to entire dataset regardless of resuming.
3. Run `06_qc2_combine_duplicates.py` to combine events with shared exceedance points occurring within 12 hours of each other. Also exports arrays with locations of all exceedance points for each event, updated 12-hr accumulation map arrays for combined events, and renumbers the event IDs. May take a few hours.


## Part IV: Additional Attributes and Auto-Classification
1. Run `07_ari_thresholds.py` to pull higher-end ARI thresholds from the NOAA Atlas 14 dataset for each event at the point of max exceedance. Takes several minutes.
2. Run `08_exceedance_array_exports.py` to export map arrays of 12-hr exceedance above the 10-yr ARI threshold for each event. Takes several minutes.
3. Run `09_exceedance_attrs.py` to calculate attributes on the exceedance arrays. Takes several minutes.
4. Run `10_pfs.py` to obtain the 1 mm precipitation feature statistics for each hour of each event. This will take several hours to run.
5. Run `11_classification.py` to apply conditions for classifying events as tropical cyclone, isolated, MCS/synoptic, and nocturnal. Also adds attributes related to the 1 mm/hr precipitation features and diurnal cycle. (Sunset/sunrise calculations use the code from https://michelanders.blogspot.com/2010/12/calulating-sunrise-and-sunset-in-python.html saved as "sun.py" to the code directory.)
6. Run `12_final_dataset.py` to organize the data table and save as a csv to be used for analysis.
7. Optional: Run `13_event_maps.py` to export image files with each event's 12-hr and peak 1-hr accumulation maps. Takes several hours, and requires a lot of RAM. Memory error likely if not splitting up. Run <2000 at a time and restart kernel before resuming.
 
Recommend backing up all exported csv files in case of accidental replacement/altering!


## Contact

Jason Chiappa - jchiappa@ou.edu

Please contact me if you have questions or would like any code used to perform analysis on the database.

Project Link: [https://github.com/jchiappa/extreme_rainfall](https://github.com/jchiappa/extreme_rainfall)

<p align="right">(<a href="#readme-top">back to top</a>)</p>
