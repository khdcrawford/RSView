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

A clinician interested in determining which RSV subtypes are most prevalent in their region would see that displayed on the RSView map. They could then select that subtype via the Subtype Selector to view where else that subtype has been recorded. To compare the prevalence of the subtype in different regions, they could use the Visualization Manager to change the dot size to represent disease prevalence. If they were in a region where RSV circulates more in the winter months, they could select that in the Year Selector. They could then use the Data Comparer to view the childhood pneumonia death analyses for this subtype in their region compared to other regions.


## Preliminary Plans

    1. Download RSV G gene sequences and appropriate metadata
    1b. Tabulate this data and organize by genotype, collection date, and country
    
    2. Download data on childhood deaths due to pneumonia
    2b. Tabulate this data and align with RSV G gene sequence dataset via aggregation and/or relabeling, as necessary
    
    3. Plot the RSV genotype data by location and year
    3b. Add interactive components to map to enable display of subsets of data (i.e. year, subtype, genotype)
    
    4. Add disease severity metrics to the interactive map display
    4b. Compare RSV genotype data with disease severity metrics.  

