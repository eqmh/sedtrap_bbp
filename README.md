# sedtrap_bbp

This repository contains:

-    particulate flux and backscattering data from sediment trap deployments in:
    - the Gulf of Mexico from December 2022 through May 2025
    - the Santa Barbara Basin from October 2016 through November of 2017
    - Bbp data contained in tsv files 'bbp_SN_YYYYMMDD.tsv' where SN is sensor serial number (1314 or 1315) and YYYYMMDD is deployment date.

-    code for data processing, compiling, and visualization and analysis:
    -   Use 'bbp_data_processing.m' for processing .RAW files downloaded from [ECO Triplet](https://www.seabird.com/eco-triplet-w/product?id=60762467721) sensors (532 and 700 nm bands)
        - This script generates a .tsv table with 'datetime', and Bbp at 532 and 700 nm
        - Requires the 'betasw_ZHH2009.m' function to calculate backscattering of seawater (Bw).
        - Inputs: 
          - filename
          - startDate and endDate
          - calibration coefficients
          - physical parameters: seawater temperature and salinity
          - TSV file name to be used, e.g. 'bbp_1315_20170510.tsv'
    -   Use 'data_compiler.m' to compile .tsv table and plot Bbp time series data
          - select directory: data_gom or data_sbb
          - TSV file name to be used for compiled table, e.g. 'eco_bbp_time_series.tsv'
    -   Use 'bbp_data_explorer.mlx' to:
        - Loads Bbp data
        - Plots Bbp data
        - Filters Bbp data
        - Smooths data
        - Loads and filters sediment trap data
        - Plots sediment trap data
        - Overlays Bbp curves onto flux plots
        - Calculates mean Bbp value for each sediment trap cup cycle
        - Creates least-square regression plots showing flux vs Bbp
    
