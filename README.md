# Plankton community composition and size structure in the North Sea in autumn and winter from 2013-2019

Script and data material related to :  
**Plankton community composition and size structure in the North Sea in autumn and winter from 2013-2019** by Gregor Börner, Romain Frelat, Anna Akimova, Cindy van Damme, Myron A. Peck, and Marta Moyano.



The data sets are stored in one Rdata files NBSS_Borner_2023.Rdata. It contains two data.frames:  
- tab: abundance of plankton per station  
- env: environmental and biological conditions per stations  
*(do we provide the data preprocessing scripts and the raw data?)*


The scripts are organised as follow:  
- 01_Cluster_Analysis.R : average abundance per station and per year, and run hierarchical clustering (Fig. 3)   
- 02_RDA_Dynamics.R : run RDA for either Buchan/Banks in autumn or Downs in winter (Fig 4 and 5)   
- 03_NASS_Correlation.R : calculate Normalized Abundance Size Spectra and correlate with evnironmental variables (Fig 6 and 7)    





[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)



