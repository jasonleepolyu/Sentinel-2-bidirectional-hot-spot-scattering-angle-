############################################################################################################ 
### This code is used to calculate the 4-parameter model scattering angles on COVE descending orbit western edges 
### Zhongbin.li@sdstate.edu
### hankui.zhang@sdstate.edu
### david.roy@sdstate.edu 
### This code is implemented in the following submitted paper:
### Investigation of Sentinel-2 bidirectional hot-spot sensing conditions, IEEE Transactions on Geoscience and Remote Sensing, 2018. 
##########################################################################################################################################
source("Descending.west.sa.4.parameter.model.r")

rm(list = ls())
if (!exists("cove_data"))  cove_data <- read.table(file="./COVE.data.2016.v3.txt", header=TRUE, sep="\t")
names(cove_data)

 ad <- cove_data$Adescending
 index <- ad=="descending" 
 print(paste("length:", sum(index)))
  
 doy <- cove_data$DOY[index];
 hour <- cove_data$UTC[index];
 year <- cove_data$Year[index];

 coordinates <- cove_data$Coordinates[index];
 coor <- as.character(coordinates)
 
 West_lon <-  vector(,length(coor))
 West_lat <-  vector(,length(coor))
 East_lon <-  vector(,length(coor))
 East_lat <-  vector(,length(coor))
 Cent_lon <-  vector(,length(coor))
 Cent_lat <-  vector(,length(coor))
 

 for (i in 1:length(coor)) {
 
	 temp <- unlist(strsplit(coor[i], ","))
	 
	 East_lon[i] <- as.numeric(temp[1])
	 East_lat[i] <- as.numeric(temp[2])
	 
	 West_lon[i] <- as.numeric(temp[7])    ### need to change for descending and ascending 
	 West_lat[i] <- as.numeric(temp[8])
	 
	 Cent_lat[i] <- (West_lat[i] + East_lat[i])/2
	  
 }
 
 ### inclination  
 inc <- 98.62 
 
 index <- Cent_lat < (180-inc) & Cent_lat > -1*(180-inc) & East_lat< (180-inc) & East_lat> -1*(180-inc) & West_lat< (180-inc) & West_lat> -1*(180-inc);
 print(paste("length:", sum(index)))
 
 coordinates <- coordinates[index]
 Clat <- Cent_lat[index]
 doy <- doy[index];
 hour <- hour[index];
 year <- year[index];
 
 
 # View azimuth angle for descending 
 j <-  vector(,length(Clat))
 
 for (i in 1:length(Clat)) {
	 if ( abs(cos(inc*pi/180)/cos(Clat[i]*pi/180)) >= 1 ) { 
		j[i] <- pi/2 
	 } else {
		j[i] <- -1 * (asin( cos(inc*pi/180)/cos(Clat[i]*pi/180) ))
	 }   
    ### in book: Orbit and Ground Track of a Satellite, Equation (5.30)
 }
 
 ### Visual azimuth angle calculation 
 ### For Western side of COVE data
 vaa <- j + pi/2
 va <- vaa * 180/pi
 
 temp <- tan( Clat * pi/180 )/tan( inc * pi/180 )
 for (i in 1:length(Clat)) {
	if( abs(temp[[i]]) > 1.0 ) {
			temp[[i]] <- round(temp[[i]],0);
		}
 }
 model_cent_lt <- 10.5 - ( asin( temp )*180/pi )/15 ;  ### Equation (3) in Ignatov et al (2004) paper: Equator crossing times for NOAA, ERS and EOS sun-synchronous satellites   
 

 half_orbit <- 145;
 earthR <- 6371.007181;
 
 ############## ITERATION COMPUTATION
 # Cent_lat <- West_lat  ## Using West_lat as the initial values of Cent_lat
 # for (i in 1:1000) { 
	# Cent_lat <- West_lat - 2*asin(cos((180-98.62)/180*pi)*half_orbit/2/earthR/cos(Cent_lat/180*pi))*180/pi;
	# Cent_lat <- Cent_lat
 # }
  
 modeled_west_lat <- 2*asin( cos((180-va)*pi/180)*half_orbit/2/earthR )*180/pi + Clat;
 
 small_r <- earthR*sin((90-abs(modeled_west_lat))*pi/180);

 delta_lon <- 2*asin( sin((180-va)*pi/180)*half_orbit/2/small_r )*180/pi;
 
 model_west_lt <- model_cent_lt - delta_lon/15;
 
 #### The beginning of the Blanco-Muriel astronomical model    
 pii <- 3.14159265358979323846
 twopi <- 2 * pii
 rad <- pii/180
 dEarthMeanRadius <- 6371.007181
 dAstronomicalUnit <- 149597890 
 
 iDay_gmt = doy;
 
 ### Calculate Year, Month, Day from DOY. Successful !!
 month <-  vector(,length(iDay_gmt))
 day <-  vector(,length(iDay_gmt))
 MonthDays <- c(31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366);
 
 for (i in 1:length(iDay_gmt)) {
 
	 m <- 0
	 
	 while ( MonthDays[m+1] < iDay_gmt[[i]] ) {
		m <- (m + 1) %% 12
	 }
	 
	 if (m >= 1) {
		d <- iDay_gmt[[i]] - MonthDays[m]
	 } else {
		d <- iDay_gmt[[i]]
	 }
	 
	 month[[i]] <- (m + 1) 
	 day[[i]] <- d
 
 }
 
 ### Julian Day
 
 liAux1 <- (month-14)/12;
 liAux2 <- (1461*(year + 4800 + liAux1))/4 + (367*(month - 2-12*liAux1))/12- (3*((year + 4900 + liAux1)/100))/4 + day - 32075;
 dJulianDate <- as.double(liAux2)  ### - 0.5  + hour/24.0;   
 ## Calculate difference between current Julian Day and JD 2451545.0 
 dElapsedJulianDays <- dJulianDate-2451545.0;
 
 ## Calculate ecliptic coordinates (ecliptic longitude and obliquity of the ecliptic in radians but without limiting the angle to be less than 2*Pi 
 
 dOmega <- 2.1429-0.0010394594 * dElapsedJulianDays;
 dMeanLongitude <- 4.8950630 + 0.017202791698 * dElapsedJulianDays; ## Radians
 dMeanAnomaly <- 6.2400600 + 0.0172019699 * dElapsedJulianDays;
 dEclipticLongitude <- dMeanLongitude + 0.03341607 * sin( dMeanAnomaly ) + 0.00034894 * sin( 2 * dMeanAnomaly ) - 0.0001134 - 0.0000203 * sin(dOmega);
 dEclipticObliquity <- 0.4090928 - 6.2140e-9 * dElapsedJulianDays + 0.0000396 * cos(dOmega);

 ## Calculate celestial coordinates ( right ascension and declination ) in radians 
 ## but without limiting the angle to be less than 2*Pi (i.e., the result may be 
 ## greater than 2*Pi)
 
 dSin_EclipticLongitude <- sin( dEclipticLongitude );
 dY <- cos( dEclipticObliquity ) * dSin_EclipticLongitude;
 dX <- cos( dEclipticLongitude );
 dRightAscension <- atan2( dY,dX );
 
 for (i in 1:length(dRightAscension)) {
	 if( dRightAscension[[i]] < 0.0 ) {
		dRightAscension[[i]] <- dRightAscension[[i]] + twopi;
	 }
 }
 
 dDeclination <- asin( sin( dEclipticObliquity )*dSin_EclipticLongitude );
 
 dLocalMeanSiderealTime <- (6.6974243242 + 0.0657098283 * dElapsedJulianDays + model_west_lt ) * 15 * rad;  ###    West_lt    hour + West_lon/15    hour + longitude/15
 
 dHourAngle <- dLocalMeanSiderealTime - dRightAscension;
 dLatitudeInRadians <- modeled_west_lat * rad;
 dCos_Latitude <- cos( dLatitudeInRadians );
 dSin_Latitude <- sin( dLatitudeInRadians );
 dCos_HourAngle <- cos( dHourAngle );
 dZenithAngle <- (acos( dCos_Latitude*dCos_HourAngle*cos(dDeclination) + sin( dDeclination )*dSin_Latitude));
 
 dY <- -sin( dHourAngle );
 dX <- tan( dDeclination )*dCos_Latitude - dSin_Latitude*dCos_HourAngle;
 dAzimuth <- atan2( dY, dX );
 
 for (i in 1:length(dAzimuth)) {
	 if ( dAzimuth[[i]] < 0.0 ) {
		dAzimuth[[i]] <- dAzimuth[[i]] + twopi;
	 }
 }
 
 dAzimuth <- dAzimuth/rad;
 
 ### Parallax Correction
 dParallax <- (dEarthMeanRadius/dAstronomicalUnit)*sin(dZenithAngle);
 
 dZenithAngle <- (dZenithAngle + dParallax)/rad;
 
 # View zenith angle
 vza <- 11.93 * pi/180 
 
 raa <- dAzimuth*pi/180 - vaa
 
 sza <- dZenithAngle *pi/180
  
 SA <- acos(-cos(sza)*cos(vza)-cos(raa)*sin(sza)*sin(vza)) * 180/pi
  
 # write.table(SA, "./COVE.data.2016.SA.reference.descending.Western.txt", sep="\t", row.names = FALSE, col.names = FALSE)
      
 
