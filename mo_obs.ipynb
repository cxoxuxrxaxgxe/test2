#################################################################################
#---  ReadMetOfficeObservations.R
#---  Loads, analyzes and cleans the MetOffice observations dataset
#---  EU WS ICM - Section "Physical Observations"
#---  Contact : GRM CAT Risks & Reinsurance Team; tristan.perotin@axa.com
#---  April 2020
#################################################################################

# Data manipulation libraries
require(data.table)
require(dplyr)
require(ncdf4)            
# Plotting librairies
require(highcharter)
require(leaflet)
require(htmlwidgets)
require(colorRamps)
require(RColorBrewer)
# Spatial libraries
require(raster)
# Statistics libraries
require(fitdistrplus)
require(evd)              # For gev distributions
require(crch)             # For truncated distributions

########################################
# 1/7. Open and read the observations ("Description of the dataset")

# Setting input and output directories
Dir_In <- "C:/Users/a-diouf/OneDrive - AXA/EU-WS-ICM/1-Wind Observations/MetOffice Observations/axa_obs_253/"
Dir_Out <- "C:/Users/a-diouf/OneDrive - AXA/EU-WS-ICM/1-Wind Observations/Figures/"


# Creating a data table of our files by looping over them
MetOfficeObs <- data.table()
Filenames <- list.files(path = Dir_In,pattern = ".nc$",recursive = T)
for (filename in Filenames) { 
  
  # Retrieving the file, its variables number and converting it to a table by adding variables column by column with cbind 
  file <- nc_open(paste0(Dir_In,filename))
  NVars <- length(file$var)  
  NC_Table <- data.table()   
  for (nvar in 2:NVars) {    # Exclusion of the first variable: the storm's time window "time_bounds"
    NC_Table <- cbind(NC_Table, ncvar_get(file, varid = file$var[[nvar]]$name)) 
  }
  
  # Setting variables names, converting times to UTC posixct and max_wind_gust to km/h
  setnames(NC_Table, sapply(2:NVars,function(nvar){file$var[[nvar]]$name}))
  NC_Table[,time := time*3600+as.POSIXct(substr(file$var[[2]]$units,13,31),tz="UTC")]
  NC_Table[,max_wind_gust := max_wind_gust*3.6]
  NC_Table[,stormname := filename]
  NC_Table[,datebeg := as.POSIXct(substr(file$var[[2]]$units,13,31),tz="UTC")+ncvar_get(file, varid = file$var[[1]]$name)[1]*3600]
  NC_Table[,dateend := as.POSIXct(substr(file$var[[2]]$units,13,31),tz="UTC")+ncvar_get(file, varid = file$var[[1]]$name)[2]*3600]
  MetOfficeObs <- bind_rows(MetOfficeObs,NC_Table)
}

# Removing duplicates and 0 gust observations (impossible value)
MetOfficeObs <- MetOfficeObs[!duplicated(MetOfficeObs[,.(time,platform_id)])] 
MetOfficeObs <- MetOfficeObs[max_wind_gust>0]                                 
MetOfficeObs
length(unique(MetOfficeObs$platform_id)) # Question pr Hugo : à enlever ?
length(unique(MetOfficeObs$stormname))

########################################
# 2/7. Plot the observation network ("Description of the dataset")

# Selection of most frequent coordinates for stations whose coordinates vary slightly
mostfreq_extract  <- function (x,...) {as.numeric(names(sort(table(x),decreasing=TRUE)[1]))}

# Description of coordinates and available observations for each platforms
Platforms         <- MetOfficeObs[, .(NObs = .N, longitude=mostfreq_extract(longitude), latitude=mostfreq_extract(latitude), 
                                      altitude=mostfreq_extract(altitude), StationStart=min(time), StationEnd=max(time)), by=platform_id]

setorder(Platforms, NObs)                      
coordinates(Platforms) <- ~longitude+latitude # convert to SpatialPointsDataFrame
pal   <- colorNumeric(colorRampPalette(matlab.like(256))(256), c(0,220),na.color = "transparent") # color palette
MAPll <- leaflet() %>% addTiles()            %>% 
  addProviderTiles( "OpenStreetMap.Mapnik" ) %>%                                                  # map background
  addLegend(pal=pal,values = c(0,220), position="bottomright") %>%                                # legend
  addCircleMarkers(data=Platforms, color="black", radius=5, weight=1, opacity=1,
                   fillColor = ~pal(NObs), fillOpacity=0.8, group="Ground Stations" ) %>%         # observations
  addLayersControl(overlayGroups=c("Ground Stations"),options=layersControlOptions(collapsed=F))  # control
saveWidget(MAPll, paste0(Dir_Out,"GroundStations.html"), selfcontained=F)                         # save leaflet as html file

########################################
# 3/7. Study the distribution of all of the observations ("Analysis of the dataset")
fitlnorm  <- fitdist(MetOfficeObs$max_wind_gust,"lnorm")
fittnorm  <- fitdist(MetOfficeObs$max_wind_gust,"tnorm", fix.arg=list(left=0), start = list(mean=50,sd=10))
fittlogis <- fitdist(MetOfficeObs$max_wind_gust,"tlogis", fix.arg=list(left=0), start = list(location=50,scale=10))
fitgamma  <- fitdist(MetOfficeObs$max_wind_gust,"gamma")
fitgev    <- fitdist(MetOfficeObs$max_wind_gust,"gev",start=list(loc=50,scale=10,shape=0))
fits      <- list(fitlnorm,fittnorm,fittlogis,fitgamma,fitgev)
GoFit     <- gofstat(fits)
GoFit     # display the goodness of fit statistics for all of the fits
# plot the top 3 fits
fits <- list(fittlogis,fitgamma,fitgev)
png(paste0(Dir_Out,"ObservationsDistribution.png"),units="in", width=5, height=5,res=300)
par(mfrow = c(2,2), mar = c(4, 4, 1, 1))
denscomp(fits, main = "", xlab = "wind gusts(km/h)", fitcol = c("seagreen","blue","red"), fitlty = c("solid"))
qqcomp(fits, addlegend = FALSE, main = "", fitpch = 16, fitcol = c("seagreen","blue","red"), line01lty = 2)
cdfcomp(fits, addlegend = FALSE, main = "", xlab = "wind gusts(km/h)", fitcol = c("seagreen","blue","red"), lines01 = TRUE, fitlty = c("solid"))
ppcomp(fits, addlegend = FALSE, main = "", fitpch = 16, fitcol = c("seagreen","blue","red"), line01lty = 2)
dev.off()

########################################
# 4/7. Study the goodness of fit for each station with a 0.05 significance KS test ("Analysis of the dataset")
Platforms              <- as.data.table(Platforms) # convert back to data table
Platforms$KSTGEV       <- "not computed"
Platforms$KSTGAM       <- "not computed"
Platforms$KSTTNM       <- "not computed"
Platforms$KSTLNM       <- "not computed"
Platforms$KSTTLG       <- "not computed"
Platforms$GEVSHAPE     <- 0
platforms <- as.integer(names(sort(table(MetOfficeObs$platform_id),decreasing=TRUE)))
for (platform in platforms) {
  PID <- which(MetOfficeObs$platform_id==platform)
  if (Platforms[platform_id==platform]$NObs>=5){
    # perform fits
    fitlnorm  <- fitdist(MetOfficeObs[PID]$max_wind_gust,"lnorm")
    fittnorm  <- fitdist(MetOfficeObs[PID]$max_wind_gust,"tnorm", fix.arg=list(left=0), start = list(mean=50,sd=10))
    fittlogis <- fitdist(MetOfficeObs[PID]$max_wind_gust,"tlogis", fix.arg=list(left=0), start = list(location=50,scale=10))
    fitgamma  <- fitdist(MetOfficeObs[PID]$max_wind_gust,"gamma")
    fitgev    <- tryCatch({fitdist(MetOfficeObs[PID]$max_wind_gust,"gev",start=list(loc=65,scale=20,shape=-0.1))},
                          error=function(error_message) {
                            return(fitdist(MetOfficeObs[PID]$max_wind_gust,"weibull",start=list(scale=50,shape=3)))}) 
    # the gev fit can produce error in rare cases caused by poor starting parameters, we force a weibull fit when that happens as most fitted gevs are of the weibull family
    # analyze fits and store results
    Platforms[platform_id==platform]$KSTGEV   <- gofstat(fitgev)$kstest[1]
    Platforms[platform_id==platform]$KSTGAM   <- gofstat(fitgamma)$kstest[1]
    Platforms[platform_id==platform]$KSTTNM   <- gofstat(fittnorm)$kstest[1]
    Platforms[platform_id==platform]$KSTLNM   <- gofstat(fitlnorm)$kstest[1]
    Platforms[platform_id==platform]$KSTTLG   <- gofstat(fittlogis)$kstest[1]
    Platforms[platform_id==platform]$GEVSHAPE <- ifelse(length(fitgev$estimate)>2,fitgev$estimate[3],0) # save the fitted gev shapes
    # show results: uncomment if you want to display results along the way
    print(Platforms[platform_id==platform])
    # print(GoFit)
    # fits <- list(fitgev,fitgamma,fittnorm,fitlnorm,fittlogis)
    # par(mfrow = c(2,2), mar = c(4, 4, 1, 1))
    # denscomp(fits, main = "", xlab = "wind gusts(km/h)", fitcol = c("seagreen","blue","red","darkorchid","dodgerblue"), fitlty = c("solid"))
    # qqcomp(fits, addlegend = FALSE, main = "", fitpch = 16, fitcol = c("seagreen","blue","red","darkorchid","dodgerblue"), line01lty = 2)
    # cdfcomp(fits, addlegend = FALSE, main = "", xlab = "wind gusts(km/h)", fitcol = c("seagreen","blue","red","darkorchid","dodgerblue"), lines01 = TRUE, fitlty = c("solid"))
    # ppcomp(fits, addlegend = FALSE, main = "", fitpch = 16, fitcol = c("seagreen","blue","red","darkorchid","dodgerblue"), line01lty = 2)
    # readline(prompt="Press [enter] to continue")
  }
}
length(which(Platforms$KSTGEV=="not rejected"))/length(which(Platforms$NObs>=30))
length(which(Platforms$KSTGAM=="not rejected"))/length(which(Platforms$NObs>=30))
length(which(Platforms$KSTTNM=="not rejected"))/length(which(Platforms$NObs>=30))
length(which(Platforms$KSTLNM=="not rejected"))/length(which(Platforms$NObs>=30))
length(which(Platforms$KSTTLG=="not rejected"))/length(which(Platforms$NObs>=30))
length(which(Platforms$GEVSHAPE  < 0))/length(which(Platforms$NObs>=5)) # Weibull types
length(which(Platforms$GEVSHAPE  > 0))/length(which(Platforms$NObs>=5)) # Frechet types
length(which(Platforms$GEVSHAPE == 0))/length(which(Platforms$NObs>=5)) # Gumbel  types

########################################
# 5/7. Detect outlier observations, inspect extreme stations  ("Analysis of the dataset")
Platforms$TLGLocation <- 0
Platforms$TLGScale    <- 0
Platforms$VOutlier    <- 0
Platforms$NOutlier    <- 0
for (platform in rev(Platforms[NObs>=5]$platform_id)) {                 
  PID <- which(MetOfficeObs$platform_id==platform)
  if (Platforms[platform_id==platform]$NObs>=5){
    # perform fits
    fittlogis <- fitdist(MetOfficeObs[PID]$max_wind_gust,"tlogis", fix.arg=list(left=0), start = list(location=50,scale=10))
    # analyze fit and store results
    Platforms[platform_id==platform]$TLGLocation <- fittlogis$estimate[1]
    Platforms[platform_id==platform]$TLGScale    <- fittlogis$estimate[2]
    # detect outliers
    Platforms[platform_id==platform]$VOutlier    <- quantile(fittlogis,probs=0.98^(1/Platforms[platform_id==platform]$NObs))[[1]][1,1]
    NOutlier                                     <- nrow(MetOfficeObs[PID][max_wind_gust>Platforms[platform_id==platform]$VOutlier])
    Platforms[platform_id==platform]$NOutlier    <- NOutlier
    # display results : uncomment lines as necessary
    print(paste0(platform,": ",NOutlier," outliers detected"))
    # plot(fittlogis)
    # readline(prompt="Press [enter] to continue")
    }
}

# Print an example of station with a high mode
fittlogis <- fitdist(MetOfficeObs[platform_id==59230]$max_wind_gust,"tlogis", fix.arg=list(left=0), start = list(location=50,scale=10))
png(paste0(Dir_Out,"Station59230Distribution.png"),units="in", width=5, height=5,res=300)
par(mfrow = c(2,2), mar = c(4, 4, 1, 1))
denscomp(fittlogis, main = "", xlab = "wind gusts(km/h)", fitcol = c("red"), fitlty = c("solid"))
qqcomp(fittlogis, addlegend = FALSE, main = "", fitpch = 16, fitcol = c("red"), line01lty = 2)
cdfcomp(fittlogis, addlegend = FALSE, main = "", xlab = "wind gusts(km/h)", fitcol = c("red"), lines01 = TRUE, fitlty = c("solid"))
ppcomp(fittlogis, addlegend = FALSE, main = "", fitpch = 16, fitcol = c("red"), line01lty = 2)
dev.off()
Platforms[platform_id==59230]
quantile(fittlogis,probs=(Platforms[platform_id==59230]$NObs - 0.5)/Platforms[platform_id==59230]$NObs)[[1]][1,1]

########################################
# 6/7. Manually inspect outlier observations and remove observations above a likelihood threshold ("Cleaning of the dataset")
for (platform in rev(Platforms[NObs>=5][NOutlier>0]$platform_id)) {
  PID <- which(MetOfficeObs$platform_id==platform)
  NObs <- nrow(MetOfficeObs[PID])
  Platforms[platform_id==platform]$NObs <- NObs
  if (NObs>=5){
    # perform fits
    fittlogis <- fitdist(MetOfficeObs[PID]$max_wind_gust,"tlogis", fix.arg=list(left=0), start = list(location=50,scale=10))
    print(paste0(platform,": ",Platforms[platform_id==platform]$NOutlier," outliers detected"))
    print(paste0("Outliers exceeding threshold by : ", 
                 MetOfficeObs[platform_id==platform][max_wind_gust>Platforms[platform_id==platform]$VOutlier]$max_wind_gust/
                   Platforms[platform_id==platform]$VOutlier, " factor"))
    # display results
    plot(fittlogis)
    readline(prompt="Press [enter] to continue")
  }
}

# Flag and examine suspicious measures
MetOfficeObs <- Platforms %>%
  .[i=MetOfficeObs, on=c("platform_id"="platform_id")]
hist(MetOfficeObs[NObs>=5][max_wind_gust>VOutlier]$max_wind_gust/MetOfficeObs[NObs>5][max_wind_gust>VOutlier]$VOutlier)
MetOfficeObs[NObs>=5][max_wind_gust/VOutlier>1.1]
# plot one example of a station with an outlier
fittlogis <- fitdist(MetOfficeObs[platform_id==21929]$max_wind_gust,"tlogis", fix.arg=list(left=0), start = list(location=50,scale=10))
png(paste0(Dir_Out,"Station21929Distribution.png"),units="in", width=5, height=5,res=300)
par(mfrow = c(2,2), mar = c(4, 4, 1, 1))
denscomp(fittlogis, main = "", xlab = "wind gusts(km/h)", fitcol = c("red"), fitlty = c("solid"))
qqcomp(fittlogis, addlegend = FALSE, main = "", fitpch = 16, fitcol = c("red"), line01lty = 2)
cdfcomp(fittlogis, addlegend = FALSE, main = "", xlab = "wind gusts(km/h)", fitcol = c("red"), lines01 = TRUE, fitlty = c("solid"))
ppcomp(fittlogis, addlegend = FALSE, main = "", fitpch = 16, fitcol = c("red"), line01lty = 2)
dev.off()
StormsWithOutliers <- table(MetOfficeObs[NObs>=5][max_wind_gust/VOutlier>1.1]$stormname)
StormsWithOutliers
StormsWithOutliers[StormsWithOutliers>=3]

# look at the 5 storms with the most outliers
SP_MetOfficeObs <- MetOfficeObs
coordinates(SP_MetOfficeObs) <- ~longitude+latitude
pal <- colorNumeric(colorRampPalette(matlab.like(256))(256), c(50,275),na.color = "transparent")
MAPll = leaflet() %>% addTiles() %>% 
  addProviderTiles( "OpenStreetMap.Mapnik" ) %>%
  addLegend(pal=pal,values = c(75,275), position="bottomright") %>%
  addCircleMarkers(data=SP_MetOfficeObs[which(SP_MetOfficeObs$stormname %in% names(StormsWithOutliers[StormsWithOutliers>=3])[1]),], color="black", radius=5, weight=1, opacity=1, fillColor = ~pal(max_wind_gust), fillOpacity=0.8, group="Unnamed storm, 03-02-1999" ) %>%
  addCircleMarkers(data=SP_MetOfficeObs[which(SP_MetOfficeObs$stormname %in% names(StormsWithOutliers[StormsWithOutliers>=3])[2]),], color="black", radius=5, weight=1, opacity=1, fillColor = ~pal(max_wind_gust), fillOpacity=0.8, group="Anatol, 02-12_1999" ) %>%
  addCircleMarkers(data=SP_MetOfficeObs[which(SP_MetOfficeObs$stormname %in% names(StormsWithOutliers[StormsWithOutliers>=3])[3]),], color="black", radius=5, weight=1, opacity=1, fillColor = ~pal(max_wind_gust), fillOpacity=0.8, group="Xynthia, 28-03-2010" ) %>%
  addCircleMarkers(data=SP_MetOfficeObs[which(SP_MetOfficeObs$stormname %in% names(StormsWithOutliers[StormsWithOutliers>=3])[4]),], color="black", radius=5, weight=1, opacity=1, fillColor = ~pal(max_wind_gust), fillOpacity=0.8, group="Friedhelm-Joachim, 14-12-2011" ) %>%
  addCircleMarkers(data=SP_MetOfficeObs[which(SP_MetOfficeObs$stormname %in% names(StormsWithOutliers[StormsWithOutliers>=3])[5]),], color="black", radius=5, weight=1, opacity=1, fillColor = ~pal(max_wind_gust), fillOpacity=0.8, group="Dirk, 25-12-2013" ) %>%
  addLayersControl(overlayGroups=c("Unnamed storm, 03-02-1999","Anatol, 02-12_1999","Xynthia, 28-03-2010","Friedhelm-Joachim, 14-12-2011","Dirk, 25-12-2013"),options=layersControlOptions(collapsed=F)) %>%
  hideGroup(c("Unnamed storm, 03-02-1999","Anatol, 02-12_1999","Xynthia, 28-03-2010","Friedhelm-Joachim, 14-12-2011","Dirk, 25-12-2013"))
saveWidget(MAPll, paste0(Dir_Out,"OutliersObsFPs.html"), selfcontained=T)
MetOfficeObs[stormname %in% names(StormsWithOutliers[StormsWithOutliers>=3])][NObs>=5][max_wind_gust/VOutlier>1.1] # show coordinates of the flagged outliers for these 5 storms

# remove outliers and summarize the remaining observations
# MetOfficeObs <- MetOfficeObs[NObs>=5][max_wind_gust/VOutlier<=1.1] # another methodology might be considered, including nearest neighbour based outlier detection techniques
MetOfficeObs
length(unique(MetOfficeObs$platform_id))
length(unique(MetOfficeObs$stormname))
table(MetOfficeObs$source) # less than 20% of the observed gusts are estimated from sustained winds 
# save cleaned observations
fwrite(MetOfficeObs[,.(max_wind_gust,time,platform_id,longitude,latitude,altitude,source,stormname)],paste0(Dir_Out,"../MetOfficeObservations.csv"))

########################################
# 7/7. Illustration of climate variability
LongTermPlatforms <- Platforms[StationStart<as.POSIXct("1990/01/01")][StationEnd>as.POSIXct("2015/01/01")][NObs >= 70] # more than 2 storms / year observed
MetOfficeObs_Long <- MetOfficeObs[platform_id %in% LongTermPlatforms$platform_id][time>=as.POSIXct("1990/01/01")][time<as.POSIXct("2015/01/01")] %>% setorder(-max_wind_gust)
MetOfficeObs_Long[, rank := rank(-max_wind_gust, ties.method = "first"), by=platform_id]

plot(MetOfficeObs_Long[rank<=7]$time, MetOfficeObs_Long[rank<=7]$max_wind_gust, main= ">5 year return period extreme measures for each station", 
     ylab = "Wind gusts (km/h)", xlab = "Time", 
     col =  rgb(red = 0, green = 0, blue = 1, alpha = 0.03), pch=16) 
# 5 year return period observations for every stations
# not included in the documentation as ideally we should select a spatially representative sample of industry exposure (the number of observations varies a lot depending on the country)
# we will look at SSI variations in uncalibrated and calibrated wind footprints instead
