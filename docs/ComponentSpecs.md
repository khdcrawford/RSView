## Component Specification

The package `RSView` analyzes RSV G gene sequences in conjunction with information on the burden of disease (as extrapolated from childhood death rates due to acute respiratory tract infection) in order to better understand global RSV circulation dynamics and the role of RSV subtype on disease severity. In order to carry out the analysis, this package must be able to:
    
    1. Download RSV G gene sequences and appropriate metadata
    2. Tabulate this data and organize by genotype, collection date, and country
    3. Plot this RSV genotype data by location and year (on an interactive map)
    4. Compare RSV genotype data with disease severity metrics.  


## Software Components

- Visualization Manager: Via dropdown menus, the User can specify options for how the data is displayed (i.e. color dots by year or by subtype; dot size representing prevalence or disease severity)

- Subtype Selector: The User provides input by selecting which subtype or genotype they wish to display from a dropdown menu. The output is a map with that subset of data.

- Year Selector: The User provides input via a) a slider to select which years to include and/or b) a checkbox to select which seasons to include (i.e. winter, summer). The output is a map with that subset of data.

- Data Comparer: This component takes the subset of data currently displayed on the map and generates a graph displaying the number of childhood deaths due to pneumonia for that region and time period 


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

