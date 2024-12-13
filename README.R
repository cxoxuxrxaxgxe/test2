#################################################################################
#---  StatisticalDownscaling
#---  Calibrates a reanalysis with observed data
#---  Contact : GRM CAT Risks & Reinsurance Team; tristan.perotin@axa.com
#---  April 2020
#################################################################################
# Useful Librairies:
# Data manipulation
require(data.table)
require(dplyr)
# Statistics
require(fitdistrplus)
require(evd)
require(crch)
require(mgcv)
# Geography
require(raster)
require(maps)
require(maptools)
require(geosphere)
# Plotting
require(highcharter)
require(htmlwidgets)
require(leaflet)
require(colorRamps)
require(RColorBrewer)


########################################
# 1.Input Data : observations, model values, land data (GMTED 2010 orography & Copernicus fractional land use)
Dir_Obs <- "C:/Users/t-perotin/OneDrive - AXA/CAT Modeling/European Windstorms/EU WS ICM - Hazard/2-Reanalysis Selection/"
Dir_Local <- "C:/Users/t-perotin/OneDrive - AXA/CAT Modeling/European Windstorms/EU WS ICM - Hazard/6-Statistical Downscaling/"
Dir_Sensitivity <- "C:/Users/t-perotin/OneDrive - AXA/CAT Modeling/European Windstorms/EU WS ICM - Hazard/6-Statistical Downscaling/SensitivityAnalysis/"
MetOfficeObs <- fread(paste0(Dir_Obs,"FinalMetOfficeObs.csv")) 
MetOfficeObs[,time:=as.POSIXct(substr(time,1,13),format="%Y-%m-%dT%H",tz="UTC")] 

Dir_Predictors <- "C:/Users/t-perotin/OneDrive - AXA/CAT Modeling/European Windstorms/EU WS ICM - Hazard/6-Statistical Downscaling/Predictors/"
Bare    <- raster(paste0(Dir_Predictors,"Bare.01.tif"))
Crop    <- raster(paste0(Dir_Predictors,"Crop.01.tif"))
Grass   <- raster(paste0(Dir_Predictors,"Grass.01.tif"))
Moss    <- raster(paste0(Dir_Predictors,"Moss.01.tif"))
Shrub   <- raster(paste0(Dir_Predictors,"Shrub.01.tif"))
Snow    <- raster(paste0(Dir_Predictors,"Snow.01.tif"))
Topo    <- raster(paste0(Dir_Predictors,"Topo.01.tif"))
Tree    <- raster(paste0(Dir_Predictors,"Tree.01.tif"))
Urban   <- raster(paste0(Dir_Predictors,"Urban.01.tif"))

gf        <- matrix(c(0,-1,0,-1,4,-1,0,-1,0),nrow=3,ncol=3)
Lapl      <- focal(Topo,gf,pad=TRUE,padValue=NA,na.rm=T)
Lapl2     <- raster(paste0(Dir_Predictors,"Lapl2.tif")) # already computed on the virtual machine, these are smoothed versions of the Laplacian
Lapl5     <- raster(paste0(Dir_Predictors,"Lapl2.tif"))
Lapl10    <- raster(paste0(Dir_Predictors,"Lapl10.tif"))
Lapl25    <- raster(paste0(Dir_Predictors,"Lapl25.tif"))
Lapl50    <- raster(paste0(Dir_Predictors,"Lapl50.tif"))
Lapl125   <- raster(paste0(Dir_Predictors,"Lapl125.tif"))
Roughness <- terrain(Topo,opt="roughness",unit="degrees",neighbors = 8)
TPI   <- terrain(Topo,opt="TPI",unit="degrees",neighbors = 8)
gf    <- focalWeight(Topo, 2e-2/2, "Gauss") # smoothed at 2km characteristic radius (comes from a quick performance analysis, could be refined)
TPI2  <- focal(TPI,gf,pad=TRUE,padValue=NA,na.rm=T)
gf    <- focalWeight(Topo, 10e-2/2, "Gauss") # smoothed at 2km characteristic radius (comes from a quick performance analysis, could be refined)
TPI10 <- focal(TPI,gf,pad=TRUE,padValue=NA,na.rm=T)
gf    <- focalWeight(Topo, 50e-2/2, "Gauss") # smoothed at 2km characteristic radius (comes from a quick performance analysis, could be refined)
TPI50 <- focal(TPI,gf,pad=TRUE,padValue=NA,na.rm=T)

gf       <- focalWeight(Topo, 10e-2/2, "Gauss") # smoothed at 2km characteristic radius (comes from a quick performance analysis, could be refined) 
Bare     <- focal(Bare,gf,pad=TRUE,padValue=NA,na.rm=T)
Crop     <- focal(Crop,gf,pad=TRUE,padValue=NA,na.rm=T)
Grass    <- focal(Grass,gf,pad=TRUE,padValue=NA,na.rm=T)
Moss     <- focal(Moss,gf,pad=TRUE,padValue=NA,na.rm=T)
Shrub    <- focal(Shrub,gf,pad=TRUE,padValue=NA,na.rm=T)
Snow     <- focal(Snow,gf,pad=TRUE,padValue=NA,na.rm=T)
Tree     <- focal(Tree,gf,pad=TRUE,padValue=NA,na.rm=T)
Urban    <- focal(Urban,gf,pad=TRUE,padValue=NA,na.rm=T)
Altitude <- focal(Topo,gf,pad=TRUE,padValue=NA,na.rm=T)

MuOri <- raster(paste0(Dir_Local,"MuOri.tif"))
SdOri <- raster(paste0(Dir_Local,"SdOri.tif"))


########################################
# 2.Extract local features at stations
#
# plot general model/obs bias
hc <- highchart() %>%
  hc_add_series(data=data.frame(x=quantile(MetOfficeObs$max_wind_gust,c(seq(0.01,0.98,0.02),0.985,0.99,0.995,0.998,0.999)),
                                y=quantile(MetOfficeObs$ERA5Gust,c(seq(0.01,0.98,0.02),0.985,0.99,0.995,0.998,0.999),na.rm=T),
                                name=c(seq(0.01,0.98,0.02),0.985,0.99,0.995,0.998,0.999)),
                name="Q/Q plot", type = "line") %>%
  hc_add_series(data=data.frame(x=c(0,250),y=c(0,250)),type = "line", name="X=Y") %>%
  hc_yAxis(title = list(text = "Modeled Wind Gust Quantiles (km/h)"), min=0, max=200) %>%
  hc_xAxis(title = list(text = "Observed Wind Gust Quantiles (km/h)"),min=0, max=200) %>%
  hc_title(text = paste0("MetOffice Modeled Wind Gusts to Observations Bias"))
hc

# extract local predictors
mostfreq_extract       <- function (x,...) {as.numeric(names(sort(table(x),decreasing=TRUE)[1]))} 
Platforms              <- MetOfficeObs[, .(NObs = .N, longitude= mostfreq_extract(longitude), latitude = mostfreq_extract(latitude), 
                                           altitude = mostfreq_extract(altitude), StationStart=min(time), StationEnd=max(time)), by = platform_id]
setorder(Platforms,NObs)                        # Description of coordinates and available observations for each platforms
coordinates(Platforms) <- ~longitude+latitude   # convert to SpatialPointsDataFrame
crs.epsg4326           <- CRS("+init=EPSG:4326")
crs(Platforms)         <- crs.epsg4326
Platforms$Altitude     <- extract(Altitude,Platforms)
Platforms$Topo         <- extract(Topo,Platforms)
Platforms$Bare         <- extract(Bare,Platforms)
Platforms$Crop         <- extract(Crop,Platforms)
Platforms$Grass        <- extract(Grass,Platforms)
Platforms$Moss         <- extract(Moss,Platforms)
Platforms$Shrub        <- extract(Shrub,Platforms)
Platforms$Snow         <- extract(Snow,Platforms)
Platforms$Tree         <- extract(Tree,Platforms)
Platforms$Urban        <- extract(Urban,Platforms)
Platforms$Lapl         <- extract(Lapl,Platforms)
Platforms$Lapl2        <- extract(Lapl2,Platforms)
Platforms$Lapl5        <- extract(Lapl5,Platforms)
Platforms$Lapl10       <- extract(Lapl10,Platforms)
Platforms$Lapl25       <- extract(Lapl25,Platforms)
Platforms$Lapl50       <- extract(Lapl50,Platforms)
Platforms$Lapl125      <- extract(Lapl125,Platforms)
Platforms$Roughness <- extract(Roughness,Platforms)
Platforms$TPI <- extract(TPI,Platforms)
Platforms$TPI2 <- extract(TPI2,Platforms)
Platforms$TPI10 <- extract(TPI10,Platforms)
Platforms$TPI50 <- extract(TPI50,Platforms)

# Compute observed / modeled distributions & local bias
ERA5Location           <- raster(paste0(Dir_Predictors,"ERA5_Location.tif"))
ERA5Scale              <- raster(paste0(Dir_Predictors,"ERA5_Scale.tif"))
Platforms$ERA5Location <- extract(ERA5Location,Platforms)
Platforms$ERA5Scale    <- extract(ERA5Scale,Platforms)
Platforms$ObsMu        <- 0
Platforms$ObsSd        <- 0
Platforms$ERA5Mu       <- 0
Platforms$ERA5Sd       <- 0
for (platform in Platforms$platform_id) {                 
  PID <- which(MetOfficeObs$platform_id==platform)
  if (Platforms[Platforms$platform_id==platform,]$NObs>=5){
    # perform fits
    fittlogis <- fitdist(MetOfficeObs[PID]$max_wind_gust,"tlogis", fix.arg=list(left=0), start = list(location=50,scale=10))
    # analyze fit and store results
    Platforms$ObsMu[Platforms$platform_id==platform] <- fittlogis$estimate[1]
    Platforms$ObsSd[Platforms$platform_id==platform] <- fittlogis$estimate[2]
    # plot(fittlogis)
    # readline(prompt="Press [enter] to continue")
    # perform fits
    fittlogis <- fitdist(MetOfficeObs[PID]$ERA5Gust,"tlogis", fix.arg=list(left=0), start = list(location=50,scale=10))
    # analyze fit and store results
    Platforms$ERA5Mu[Platforms$platform_id==platform] <- fittlogis$estimate[1]
    Platforms$ERA5Sd[Platforms$platform_id==platform] <- fittlogis$estimate[2]
    # plot(fittlogis)
    # readline(prompt="Press [enter] to continue")
    print(Platforms[Platforms$platform_id==platform,])
  }
}
Platforms$DMU <- Platforms$ObsMu-Platforms$ERA5Mu
Platforms$DSD <- Platforms$ObsSd-Platforms$ERA5Sd
Platforms
save(Platforms,file=paste0(Dir_Local,"Platforms_Enriched.RData"))


########################################
# 3. Fit a logistic model predicting model bias based on predictor values at station locations
#
# represent relationship between bias and Laplacian of orography (most important predictor)
smoothScatter(Platforms$Lapl2,Platforms$DMU,nrpoints=0) # logistic regression appears appropriate

logis <- function(x){1/(1+exp(-x))}
logit <- function(x){log(x/(1-x))}
fit   <- nls(DMU~a*(logis(TPI2*b+c))+d+e*Tree+f*Urban+g*Crop+i*Grass+j*logis(TPI10*h+k)+l*TPI50,data=Platforms@data,weights=NObs,start=list(a=30,d=4,b=0.1,c=-2,e=-0.2,f=-0.1,g=-0.06,h=1,i=-0.1,j=40,k=-3,l=-30),control=list(warnOnly=T))
Platforms$PredictDMU2 <- predict(fit,Platforms@data,type="response")
summary(fit)

DT    <- data.table(TPI2=TPI2[],TPI10=TPI10[],TPI50=TPI50[],
                 Crop=Crop[],Grass=Grass[],Tree=Tree[],Urban=Urban[])
remove(Altitude,Bare,Moss,Lapl,Lapl2,Lapl5,Lapl10,Lapl25,Lapl50,Lapl125,Shrub,Snow,Roughness)
gc()

DMU   <- Topo
DMU[] <- predict(fit,DT,type = "response")

plot(DMU)
plot(DMU,xlim=c(0,5),ylim=c(47,50))
plot(DMU,xlim=c(2.5,7.5),ylim=c(43.5,46.5))
plot(DMU,xlim=c(0,10),ylim=c(50,55))

fit   <- nls(DSD~a*(logis(TPI2*b+c))+d+e*logis(TPI10*f+g)+h*TPI50+k*Urban+l*Tree+m*Crop+n*Grass,data=Platforms@data,weights=NObs,start=list(a=4,d=4,b=0.1,c=-2,e=5,f=2,g=-6,h=-3,k=-0.03,l=-0.04,m=-0.01,n=-0.02),control=list(warnOnly=T))
Platforms$PredictDSD2 <- predict(fit,Platforms@data,type="response")
summary(fit)

DSD   <- Topo
DSD[] <- predict(fit,DT,type = "response")
gc()
plot(DSD)
plot(DSD,xlim=c(0,5),ylim=c(47,50))
plot(DSD,xlim=c(2.5,7.5),ylim=c(43.5,46.5))
plot(DSD,xlim=c(0,10),ylim=c(50,55))

gc()

writeRaster(DMU,paste0(Dir_Local,"DMU.tif"),overwrite=T)
writeRaster(DSD,paste0(Dir_Local,"DSD.tif"),overwrite=T)

########################################
# 4. Extrapolate global bias (& before / after distributions)
MuOri <- raster(paste0(Dir_Predictors,"ERA5_Location.tif"))
SdOri <- raster(paste0(Dir_Predictors,"ERA5_Scale.tif"))
MuOri <- projectRaster(MuOri,DMU)
SdOri <- projectRaster(SdOri,DSD)
writeRaster(MuOri,paste0(Dir_Local,"MuOri.tif"))
writeRaster(SdOri,paste0(Dir_Local,"SdOri.tif"))
MuFinal <- MuOri+DMU
MuFinal[MuFinal<15] <- 15
plot(MuFinal,zlim=c(15,115))
plot(MuOri,zlim=c(15,115))
SdFinal <- SdOri+DSD
SdFinal[SdFinal<3.5] <- 3.5
plot(SdFinal,zlim=c(3.5,63.5))
plot(SdOri,zlim=c(3.5,63.5))
writeRaster(MuFinal,paste0(Dir_Local,"MuFinal.tif"))
writeRaster(SdFinal,paste0(Dir_Local,"SdFinal.tif"))
remove(DMU,DSD,MuOri,SdOri,MuFinal,SdFinal,ERA5Location,ERA5Scale)
gc()

for (res in c(2,5,10,20,50,100)){
MuOriUp  <- aggregate(MuOri,fact=res)
SdOriUp   <- aggregate(SdOri,fact=res)
MuFinalUp <- aggregate(MuFinal,fact=res)
SdFinalUp <- aggregate(SdFinal,fact=res)
writeRaster(MuOriUp,paste0(Dir_Sensitivity,"MuOri",res,".tif"),overwrite=T)
writeRaster(SdOriUp,paste0(Dir_Sensitivity,"SdOri",res,".tif"),overwrite=T)
writeRaster(MuFinalUp,paste0(Dir_Sensitivity,"MuFinal",res,".tif"),overwrite=T)
writeRaster(SdFinalUp,paste0(Dir_Sensitivity,"SdFinal",res,".tif"),overwrite=T)}

########################################
# Validate downscaling
MuOri <- raster(paste0(Dir_Local,"MuOri.tif"))
SdOri <- raster(paste0(Dir_Local,"SdOri.tif"))
MuFinal <- raster(paste0(Dir_Local,"MuFinal.tif"))
SdFinal <- raster(paste0(Dir_Local,"SdFinal.tif"))

load(paste0(Dir_Local,"Platforms_Enriched.RData"))
Platforms$MuOri <- extract(MuOri,Platforms)
Platforms$MuFinal <- extract(MuFinal,Platforms)
Platforms$SdOri <- extract(SdOri,Platforms)
Platforms$SdFinal <- extract(SdFinal,Platforms)

mapping <- function(speed,muori,sdori,mufin,sdfin){qtlogis(
  ptlogis(speed,location=muori,scale=sdori,left=0),
  location=mufin,scale=sdfin,left=0)}

MetOfficeObs <- fread(paste0(Dir_Obs,"FinalMetOfficeObs.csv")) 
MetOfficeObs[,time:=as.POSIXct(substr(time,1,13),format="%Y-%m-%dT%H",tz="UTC")] 
MetOfficeObs <- as.data.table(Platforms)[,.(platform_id,MuOri,MuFinal,SdOri,SdFinal,ObsMu,ObsSd)][MetOfficeObs,on=.(platform_id)]

MetOfficeObs[,Gustcalib:=mapping(ERA5Gust,MuOri,SdOri,MuFinal,SdFinal)]
nrow(MetOfficeObs[!is.na(Gustcalib)])/nrow(MetOfficeObs)
cor(MetOfficeObs[!is.na(Gustcalib)]$Gustcalib,MetOfficeObs[!is.na(Gustcalib)]$max_wind_gust)
sd(MetOfficeObs[!is.na(Gustcalib)]$Gustcalib-MetOfficeObs[!is.na(Gustcalib)]$max_wind_gust)
exp(sd(log(MetOfficeObs[!is.na(Gustcalib)]$Gustcalib/MetOfficeObs[!is.na(Gustcalib)]$max_wind_gust)))-1

for (res in c(2,5,10,20,50,100)){
  
  MuOri   <- raster(paste0(Dir_Sensitivity,"MuOri",res,".tif"))
  SdOri   <- raster(paste0(Dir_Sensitivity,"SdOri",res,".tif"))
  MuFinal <- raster(paste0(Dir_Sensitivity,"MuFinal",res,".tif"))
  SdFinal <- raster(paste0(Dir_Sensitivity,"SdFinal",res,".tif"))
  
  Platforms$MuOri <- extract(MuOri,Platforms)
  Platforms$MuFinal <- extract(MuFinal,Platforms)
  Platforms$SdOri <- extract(SdOri,Platforms)
  Platforms$SdFinal <- extract(SdFinal,Platforms)
  
  MetOfficeObs <- fread(paste0(Dir_Obs,"FinalMetOfficeObs.csv")) 
  MetOfficeObs[,time:=as.POSIXct(substr(time,1,13),format="%Y-%m-%dT%H",tz="UTC")] 
  MetOfficeObs <- as.data.table(Platforms)[,.(platform_id,MuOri,MuFinal,SdOri,SdFinal,ObsMu,ObsSd)][MetOfficeObs,on=.(platform_id)]
  
  MetOfficeObs[,Gustcalib:=mapping(ERA5Gust,MuOri,SdOri,MuFinal,SdFinal)]
  print(res)
  print(cor(MetOfficeObs[!is.na(Gustcalib)]$Gustcalib,MetOfficeObs[!is.na(Gustcalib)]$max_wind_gust))
  print(sd(MetOfficeObs[!is.na(Gustcalib)]$Gustcalib-MetOfficeObs[!is.na(Gustcalib)]$max_wind_gust))
  
}

########################################
# Display downscaling
# cf 8-Plots
MuOri <- raster(paste0(Dir_Local,"MuOri.tif"))
SdOri <- raster(paste0(Dir_Local,"SdOri.tif"))
MuFinal <- raster(paste0(Dir_Local,"MuFinal.tif"))
SdFinal <- raster(paste0(Dir_Local,"SdFinal.tif"))
Platforms$MuOri <- extract(MuOri,Platforms)
Platforms$MuFinal <- extract(MuFinal,Platforms)
Platforms$SdOri <- extract(SdOri,Platforms)
Platforms$SdFinal <- extract(SdFinal,Platforms)
MetOfficeObs <- as.data.table(Platforms)[,.(platform_id,MuOri,MuFinal,SdOri,SdFinal,ObsMu,ObsSd)][MetOfficeObs[,.(platform_id,max_wind_gust,time,longitude,latitude,altitude,stormname,ERA5Gust)],on=.(platform_id)]
MetOfficeObs[,Gustcalib1:=mapping(ERA5Gust,MuOri,SdOri,MuFinal,SdFinal)]
res <- 5
MuOri   <- raster(paste0(Dir_Sensitivity,"MuOri",res,".tif"))
SdOri   <- raster(paste0(Dir_Sensitivity,"SdOri",res,".tif"))
MuFinal <- raster(paste0(Dir_Sensitivity,"MuFinal",res,".tif"))
SdFinal <- raster(paste0(Dir_Sensitivity,"SdFinal",res,".tif"))
Platforms$MuOri <- extract(MuOri,Platforms)
Platforms$MuFinal <- extract(MuFinal,Platforms)
Platforms$SdOri <- extract(SdOri,Platforms)
Platforms$SdFinal <- extract(SdFinal,Platforms)
MetOfficeObs <- as.data.table(Platforms)[,.(platform_id,MuOri,MuFinal,SdOri,SdFinal,ObsMu,ObsSd)][MetOfficeObs[,.(platform_id,max_wind_gust,time,longitude,latitude,altitude,stormname,ERA5Gust,Gustcalib1)],on=.(platform_id)]
MetOfficeObs[,Gustcalib5:=mapping(ERA5Gust,MuOri,SdOri,MuFinal,SdFinal)]

hc <- highchart() %>%
  hc_add_series(data=data.frame(x=quantile(MetOfficeObs$max_wind_gust,c(seq(0.01,0.98,0.02),0.985,0.99,0.995,0.998,0.999)),
                                y=quantile(MetOfficeObs$ERA5Gust,c(seq(0.01,0.98,0.02),0.985,0.99,0.995,0.998,0.999),na.rm=T),
                                name=c(seq(0.01,0.98,0.02),0.985,0.99,0.995,0.998,0.999)),
                name="Q/Q plot, ERA5", type = "line") %>%
  hc_add_series(data=data.frame(x=c(0,250),y=c(0,250)),type = "line", name="X=Y") %>%
  hc_add_series(data=data.frame(x=quantile(MetOfficeObs$max_wind_gust,c(seq(0.01,0.98,0.02),0.985,0.99,0.995,0.998,0.999)),
                                y=quantile(MetOfficeObs$Gustcalib5,c(seq(0.01,0.98,0.02),0.985,0.99,0.995,0.998,0.999),na.rm=T),
                                name=c(seq(0.01,0.98,0.02),0.985,0.99,0.995,0.998,0.999)),
                name="Q/Q plot, Downscaled-5km", type = "line") %>%
  hc_add_series(data=data.frame(x=quantile(MetOfficeObs$max_wind_gust,c(seq(0.01,0.98,0.02),0.985,0.99,0.995,0.998,0.999)),
                                y=quantile(MetOfficeObs$Gustcalib1,c(seq(0.01,0.98,0.02),0.985,0.99,0.995,0.998,0.999),na.rm=T),
                                name=c(seq(0.01,0.98,0.02),0.985,0.99,0.995,0.998,0.999)),
                name="Q/Q plot, Downscaled-1km", type = "line") %>%
  hc_yAxis(title = list(text = "Modeled Wind Gust Quantiles (km/h)"), min=0, max=200) %>%
  hc_xAxis(title = list(text = "Observed Wind Gust Quantiles (km/h)"),min=0, max=200) %>%
  hc_title(text = paste0("MetOffice Modeled Wind Gusts to Observations Bias"))
hc


########################################
# Calibration
Dir_In <- "C:/Users/t-perotin/OneDrive - AXA/CAT Modeling/European Windstorms/EU WS ICM - Hazard/4-Calibration & Downscaling/UncalibratedTestEvents/"
Dir_Out <- "C:/Users/t-perotin/OneDrive - AXA/CAT Modeling/European Windstorms/EU WS ICM - Hazard/4-Calibration & Downscaling/CalibratedTestEvents/"
Files <- list.files(Dir_In,pattern=".tif$")
for (file in Files){
  print(file)
  Gust <- raster(paste0(Dir_In,file))
  Gust <- projectRaster(Gust,ERA5Location)
  Gust[] <- mapping(Gust[],ERA5Location[],ERA5Scale[],CorrectedLocation[],CorrectedScale[])
  writeRaster(Gust,paste0(Dir_Out,file))
}
