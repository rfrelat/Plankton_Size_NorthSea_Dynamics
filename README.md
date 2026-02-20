# Plankton community composition and size structure in the North Sea

Script and data material related to :  

Börner G, Frelat R, Akimova A, van Damme C, Peck MA, Moyano M (2025) **Autumn and winter plankton composition and size structure in the North Sea**. *Marine Ecology Progress Series*, 753:1-18 [DOI 10.3354/meps14767](https://doi.org/10.3354/meps14767)  


The data sets are stored in one `Rdata` files [`NBSS_Borner_2023.Rdata`](https://github.com/rfrelat/Plankton_Size_NorthSea_Dynamics/raw/main/NBSS_Borner_2023.Rdata). It contains two data.frames:  
- tab: abundance of plankton per station  
- env: environmental and biological conditions per stations  


The scripts are organized as follow:  
- 01_Cluster_Analysis.R : average abundance per station and per year, and run hierarchical clustering (Fig. 3)   
- 02_RDA_Dynamics.R : run RDA for either Buchan/Banks in autumn or Downs in winter (Fig 4 and 5)   
- 03_NASS_Correlation.R : calculate Normalized Abundance Size Spectra and correlate with environmental variables (Fig 6 and 7)    





[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13616727.svg)](https://doi.org/10.5281/zenodo.13616727)


