# Sentinel-2-bidirectional-hot-spot-scattering-angle-
### Code for calculation of Sentinel-2 swath edge scattering angles and investigation of Sentinel-2 bidirectional hot-spot sensing condition 

### Contact information:

Dr. Zhongbin Li (zhongbin.li@sdstate.edu), GSCE, SDSU, US

Dr. Hankui K. Zhang (hankui.zhang@sdstate.edu), GSCE, SDSU, US

Dr. David P. Roy (david.roy@sdstate.edu), GSCE, SDSU, US

### Usage:

#### 4 parameter model 

- implement "Descending.west.sa.4.parameter.model.r" to calculate the Sentinel-2 scattering angle on descending orbit western edges using 4 parameter Sentinel-2 specific model; 
- implement "Descending.east.sa.4.parameter.model.r" to calculate the Sentinel-2 scattering angle on descending orbit eastern edges using 4 parameter Sentinel-2 specific model;
- implement "Ascending.west.sa.4.parameter.model.r" to calculate the Sentinel-2 scattering angle on ascending orbit western edges using 4 parameter Sentinel-2 specific model;
- implement "Ascending.east.sa.4.parameter.model.r" to calculate the Sentinel-2 scattering angle on ascending orbit eastern edges using 4 parameter Sentinel-2 specific model;

#### 6 parameter model 

- implement "COVE.SA.6.parameter.model.r" to calculate the Sentinel-2 scattering angle on descending and ascending orbit western and eastern edges using 6 parameter model with time and location information provided by COVE KML files;

##### NOTE some variables in "COVE.SA.6.parameter.model.r": 

##### Descending orbit or ascending orbit 
- index <- ad=="ascending" ### for COVE Sentinel-2 ascending orbit granules 
- index <- ad=="descending"  ### for COVE Sentinel-2 descending orbit granules 

##### Western edges 
- West_lon[i] <- as.numeric(temp[3]) ### for COVE Sentinel-2  ascending orbit western edges    
- West_lat[i] <- as.numeric(temp[4]) ### for COVE Sentinel-2  ascending orbit western edges 

- West_lon[i] <- as.numeric(temp[7]) ### for COVE Sentinel-2  descending orbit western edges 
- West_lat[i] <- as.numeric(temp[8]) ### for COVE Sentinel-2  descending orbit western edges 

##### Visual azimuth angle 
- vaa <- j + 3*pi/2   ### for descending eastern edges 
- vaa <- j + pi/2     ### for descending western edges
- vaa <- pi/2 - j     ### for ascending western edges  
- vaa <- 3 * pi/2 - j ### for ascending eastern edges   

References:

Li, Z., Zhang, H.K., & Roy, D.P., (2018). Investigation of Seninel-2 bidirectional hot-spot sensing conditions. 
IEEE Transactions on Geoscience and Remote Sensing (under review).  

