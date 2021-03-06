library(gstat)
library(maptools)
library(rgdal)
library(latticeExtra)
library(e1071)
library(spacetime)
library(rgeos)
library(raster)


####### Date function ########
to_POSIX = function(date) {
  day = substring(date, 1, 2)
  month = substring(date, 3, 4)
  year = "2009"
  tstr = paste(year, "-", month, "-", day, sep = "")
  time = "00:00:00 CEST"
  return (paste(tstr, time))
}

########## Read workspaces and create CSV ##############
#setwd('Data//st_airbase_mday_04and06_below800')
wfiles = list.files(pattern="*.RData")
numfiles = length(wfiles)
start_index = nchar('"st_airbase_mday"')
cnames = t(as.matrix(colnames(st.airbase.mday)))
write.table(cnames, '..//AirData.csv', sep = ",", col.names = FALSE, row.names = FALSE)
for (i in 1:numfiles) {
  load(wfiles[i])
  date = substr(wfiles[i], start_index, unlist(gregexpr(pattern ='.RData', wfiles[i])) -1)
  date = to_POSIX(date)
  st.airbase.mday$TIME = rep(date, length(st.airbase.mday$mNo))
  write.table(st.airbase.mday, '..//AirData.csv', append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
}

########## Bounding Box ##############
abpm10 = read.csv('..//AirData.csv')
max_lat = max(abpm10$lat)
min_lat = min(abpm10$lat)
max_lon = max(abpm10$lon)
min_lon = min(abpm10$lon)


################# SPATIO-TEMPORAL ANALYSIS CODE  ############
abpm10 = read.csv('Data//AirData.csv')
coordinates(abpm10) = ~easting + northing
proj4string(abpm10) = CRS("+init=epsg:3035 +units=km")
plot(abpm10, axes=T)
bubble(abpm10, zcol="pm10.obs", scales=list(draw=T), main="Bubble Plot of ABPM10")

### HISTOGRAMS AND LOG TRANSFORMATION ###
abpm10_log = abpm10
abpm10_log$pm10.obs = log(abpm10$pm10.obs)
X11()
par(mfrow=c(1,2))
hist(abpm10$pm10.obs, freq = F, xlab = 'Data', ylab = 'Probability Density', main = 'Original Dataset')
lines(density(abpm10$pm10.obs), col= 'red')
hist(abpm10_log$pm10.obs, freq=F, xlab = 'Data', ylab = 'Probability Density', main = 'Natural Log Transformed Dataset')
lines(density(abpm10_log$pm10.obs), col= 'blue')
sk_abpm10 = skewness(abpm10$pm10.obs)
sk_abpm10_log = skewness(abpm10_log$pm10.obs)
X11()
par(mfrow=c(2,2))
boxplot(abpm10$pm10.obs, main="Original dataset")
boxplot(abpm10_log$pm10.obs, main = "Natural Log Transformed Dataset")
qqnorm(abpm10$pm10.obs, col='red')
qqline(abpm10$pm10.obs)
qqnorm(abpm10_log$pm10.obs, col='blue')
qqline(abpm10_log$pm10.obs)

######### SPATIO-TEMPORAL KRIGING ##############

#### SPATIAL PART ####
abpm10.UTM = spTransform(abpm10_log,CRS("+init=epsg:3035 +units=km")) 
abpm10_sp = SpatialPoints(abpm10.UTM@coords,CRS("+init=epsg:3035 +units=km")) 
#dupl = zerodist(abpm10_sp)
#abpm10_df = data.frame(pm10.obs=abpm10.UTM$pm10.obs[-dupl[,2]])
abpm10_df = data.frame(pm10.obs=abpm10.UTM$pm10.obs)

#### TEMPORAL PART ####

#abpm10_tm = as.POSIXct(abpm10.UTM$TIME[-dupl[,2]],tz="CET")
abpm10_tm = as.POSIXct(abpm10.UTM$TIME,tz="CET")
timeDF = STIDF(abpm10_sp, abpm10_tm, data = abpm10_df)
stplot(timeDF)

var = variogramST(pm10.obs~1, data=timeDF, tunit="days", na.omit=T) 
plot(var, map = F, xlab = "distance (km)")
plot(var,wireframe=T) 


###### SPATIOTEMPORAL MODELS #########

pars.l <- c(sill.s = 0, range.s = 500, nugget.s = 0, sill.t = 0, range.t = 1, nugget.t = 0, sill.st = 0, range.st = 100, nugget.st = 0)
pars.u <- c(sill.s = 0.9, range.s = 400, nugget.s = 0.125,sill.t = 0.9, range.t = 10, nugget.t = 0.125, sill.st = 0.6, range.st = 123, nugget.st = 0.1) 


#### SEPARABLE MODEL ####
separable = vgmST("separable", space = vgm(0.9,"Exp", 400, 0.125),time = vgm(0.9,"Exp", 10, 0.125), sill=0.6)  
separable_Vgm = fit.StVariogram(var, separable, fit.method=11,method="L-BFGS-B", stAni=5, lower=pars.l,upper=pars.u)
extractPar(separable_Vgm)
plot(var, separable_Vgm, xlab = "distance (km)", map=F)

#### PRODUCTSUM MODEL ####
prodSumModel = vgmST("productSum",space = vgm(0.9, "Exp", 400, 0.125),time = vgm(0.9, "Exp", 10, 0.125),k = 0.55) 
prodSumModel_Vgm = fit.StVariogram(var, prodSumModel,method = "L-BFGS-B",lower=pars.l)
extractPar(prodSumModel_Vgm)
plot(var, prodSumModel_Vgm, xlab = "distance (km)", map=F)

#### METRIC MODEL ####
metric = vgmST("metric", joint = vgm(0.9,"Exp", 400, 0.125), stAni=50)
metric_Vgm = fit.StVariogram(var, metric, method="L-BFGS-B",lower=pars.l)
extractPar(metric_Vgm)
plot(var, metric_Vgm, xlab = "distance (km)", map=F)

#### SUM-METRIC MODEL ####
sumMetric = vgmST("sumMetric", space = vgm(0.9,"Exp", 400, 0.125),time = vgm(0.9,"Exp", 10, 0.125), joint = vgm(0.9,"Exp", 123, nugget=0.125), stAni=50)
sumMetric_Vgm = fit.StVariogram(var, sumMetric, method="L-BFGS-B",lower=pars.l,upper=pars.u,tunit="days")
extractPar(sumMetric_Vgm)
plot(var, sumMetric_Vgm, xlab = "distance (km)", map=F)

#### COMPARE MODELS ####
attr(separable_Vgm, "MSE")
attr(prodSumModel_Vgm, "MSE")
attr(metric_Vgm, "MSE")
attr(sumMetric_Vgm, "MSE")



plot(var, wireframe=T, scales=list(arrows=F, z = list(distance = 5)), xlab = "distance (km)", main='Emperical Model')
plot(var, separable_Vgm, wireframe = T, scales=list(arrows=F, z = list(distance = 5)), main='Separable Model', xlab = "distance (km)")
plot(var, sumMetric_Vgm, wireframe=T, scales=list(arrows=F, z = list(distance = 5)), main='SumMetric Model', xlab = "distance (km)")
plot(var,metric_Vgm,wireframe=T, scales=list(arrows=F, z = list(distance = 5)), main='Metric Model', xlab = "distance (km)") 
plot(var,prodSumModel_Vgm,wireframe=T, scales=list(arrows=F, z = list(distance = 5)), main='Product Sum Model', xlab = "distance (km)") 

#### GENERATE MAPS ####
xy = expand.grid(x=seq(2500, 5500, by=100), y=seq(1500,3700, by=100))
xys = SpatialPoints(xy)
gridded(xys) = TRUE
proj4string(xys) = CRS("+init=epsg:3035 +units=km")
timeDF = spTransform(timeDF, CRS("+init=epsg:3035 +units=km"))
plot(xys, axes=T)
points(as.data.frame(abpm10_log)[,4:5], col=2, pch=19)


eu = shapefile("Data//Euro//EClip.shp")
eu.UTM = spTransform(eu, CRS("+init=epsg:3035 +units=km"))
#abpm10_log = spTransform(abpm10_log, CRS("+init=epsg:3035 +units=km"))
#plot(eu)
#plot(abpm10.UTM,add=T,col="red") 
tm.grid = seq(as.POSIXct('2009-04-15 00:00 CET'),as.POSIXct('2009-05-31 00:00 CET'), length.out=5)
#sp.grid.UTM = spsample(eu.UTM, n=2000, type="random") 
grid.ST = STF(xys, tm.grid) 
spl1 = list("sp.polygons", eu.UTM, first = FALSE) 
spl2 = list("sp.points", pch="+", abpm10_log, col=2, cex=1)
spl = list(spl1, spl2)
labels1 <- layer(sp.text(coordinates(eu.UTM), txt = eu.UTM$CNTR_CODE, pos = 1))

#### INTERPOLATE ####
predicted <- krigeST(pm10.obs~1, data=timeDF, modelList=prodSumModel_Vgm, newdata=grid.ST, computeVar = TRUE) 
predicted$var1.pred = exp(predicted$var1.pred)
#abpm10.ST <- STFDF(xys, tm.grid, data.frame(pm10.inter = predicted$var1.pred))
stplot(predicted, sp.layout = spl, scales=list(draw=TRUE), main="Kriged prediction") + labels1
varmap = predicted
varmap$var1.pred = predicted$var1.var
stplot(varmap, sp.layout = spl, scales=list(draw=TRUE), main="Kriged prediction") + labels1
