## Component Specification

The package `RSView` analyzes RSV G gene sequences in conjunction with information on the burden of disease (as extrapolated from childhood death rates due to acute respiratory tract infection) in order to better understand global RSV circulation dynamics and the role of RSV subtype on disease severity. In order to carry out the analysis, this package must be able to:
    
    1. Download RSV G gene sequences and appropriate metadata
    2. Tabulate this data and organize by genotype, collection date, and country
    3. Plot this RSV genotype data by location and year (on an interactive map)
    4. Compare RSV genotype data with disease severity metrics.  


## Software Components

- Argument Parser: parsearguments.py will handle user input

- Sequence Downloader: with seq_download.py, the user can download the RSV genotype and subtype data from GenBank and process it into a usable format

- Genotype Assigner: genotype.py will assign genotypes to RSV G sequences. This program relies on some sequences already being genotyped. It will not add sequences to a genotype that does not already have a high quality (< 60 gaps in the alignment) sequence in the input data.

- Health Data Processor: health_download.py will download and process the data on deaths resulting from acute respiratory infection into a usable dataframe.

- RSV Mapper: map_rsv.py will map the global distribution of Respiratory Syncytial Virus (RSV) by collection date, location,
and viral subtype or genotype.

- Health Data Plotter: plot_rsv.py will plot health metrics from the health data set either using summary data for each country or yearly data for a specified country.

- Health:Subtype Plotter: plot_correlation.py will integrate the health and RSV datasets and plot health metrics as a function of the relative prevalence of subtypes A and B in that country.

## Interactions to Accomplish Use Cases

### Use Case 1: Analyze global distribution of RSV genotypes or subtypes

    1. Argument Parser will handle user input for the following steps:
    2. Sequence Downloader will take arguments from Argument Parser to down load sequence data from GenBank
    3. Genotype Assigner will take arguments from the Argument Parser and data downloaded by the Sequence Downloader to assign genotypes to sequences and add them to the dataframe.
    4. RSV Mapper will then take the downloaded sequences from Sequence Downloader and genotypes from Genotype Assigner and plot them on a world map

### Use Case 2 and 3: Analyze health impact of acute respiratory infections around the world

    1. Argument Parser will handle user input for the following steps:
    2. Health Data Processor will download the health metrics dataset and process it into a usable format
    3. Health Data Processor will generate interactive graphs for analyzing health data on a global scale or on a country-specific, yearly basis

### Use Case 4: Analyze the correlation between RSV subtype prevalence and health impact

    1. Argument Parser will handle user input for the following steps:
    2. Sequence Downloader will take arguments from Argument Parser to down load sequence data from GenBank
    3. Genotype Assigner will take arguments from the Argument Parser and data downloaded by the Sequence Downloader to assign genotypes to sequences and add them to the dataframe.
    4. Health Data Processor will download the health metrics dataset and process it into a usable format
    5. Health:Subtype Plotter will integrate the health and sequence datasets and generate interactive graphs for analysis of the relationship between RSV subtype prevalence and health impact in different countries, both on a summary level and as these values change from year to year.

## Preliminary Plans

    1. Download RSV G gene sequences and appropriate metadata
    1b. Tabulate this data and organize by genotype, collection date, and country
    
    2. Download data on childhood deaths due to pneumonia
    2b. Tabulate this data and align with RSV G gene sequence dataset via aggregation and/or relabeling, as necessary
    
    3. Plot the RSV genotype data by location and year
    3b. Add interactive components to map to enable display of subsets of data (i.e. year, subtype, genotype)
    
    4. Add disease severity metrics to the interactive map display
    4b. Compare RSV genotype data with disease severity metrics.  

