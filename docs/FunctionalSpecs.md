## Functional Specifications
 
###Background

###User Profile

RSView is an analysis of the global distribution of RSV subtypes intended for laboratory scientists and clinicians. In
order to access the RSView, users will navigate to the Github repository where it is stored
(https://github.com/khdusenbury/viralseq_mapping). Though some of the intended users will have a programming background,
this is not necessary in order for the user to understand or interact with the analysis. Users that are fluent in Python
will be able to view the code used to generate the interactive subtype distribution map, plots of correlations between
RSV subtypes and childhood deaths, and other analyses.

###Data Sources

RSView utilizes publicly available viral sequence data and metadata from
[GenBank](https://www.ncbi.nlm.nih.gov/genbank/). This data is downloaded in .fasta format. Specifically, the genetic
sequences of tens of thousands of RSV viruses are used to infer genotype and subtype based on the viral G protein.
Metadata that associated with these sequences is used to locate the virus temporally and geographically. The data on
childhood deaths due to pneumonia is taken from [UNICEF](https://data.unicef.org/topic/child-health/pneumonia/). This
data is available in tabular format.

###Use Cases

RSView is intended to be used by clinicians and laboratory scientists.

We expect clinicians will be interested in using RSView to understand how certain genotypes and subtypes of RSV may
correlate with disease severity and deaths due to other respiratory infections, such as pneumonia. Clinicians may be
interested in comparing these health metrics in their geographic region to those in disparate regions where the same
subtypes of RSV have circulated. To do this, clinicians would use RSView's interactive map to determine the most
prevalent subtypes of RSV in their region and what other locations have seen this subtype. Clinicians would also use the
RSV subtype/childhood pneumonia death analyses to compare health outcomes in their country to these other regions.

Laboratory scientists studying RSV or similar viruses may be interested in RSView as a method for understanding how
viral diversity is distributed globally and how it fluctuates over time. Scientists studying RSV may primarily use
RSView's interactive map as a preliminary tool for understanding and visualizing evolution of RSV. Scientists studying
related viruses may be interested in adapting RSView to visualize the global and temporal distribution of their virus of
interest.     