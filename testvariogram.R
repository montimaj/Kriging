library(gstat)
library(maptools)
library(rgdal)
load("AirBasePM10_070409.RData")

head(abpm10)
coordinates(abpm10) <- ~easting + northing
proj4string(abpm10) <- CRS("+init=epsg:3035 +units=km")
plot(abpm10, axes=T)
bubble(abpm10, zcol="pm10.obs", scales=list(draw=T), main="Bubble Plot of ABPM10")

############ Q1 ###################
temp = abpm10
temp$pm10.obs = log10(abpm10$pm10.obs)

par(mfrow=c(1,2))
hist(abpm10$pm10.obs, freq=F, main = 'ORIGINAL ABPM10', xlab = 'Obs')
lines(density(abpm10$pm10.obs), col='red')
hist(temp$pm10.obs, freq=F, main = 'LOG TRANSFORMED ABPM10', xlab = 'Obs')
lines(density(temp$pm10.obs), col='blue')
skewness(abpm10$pm10.obs)
skewness(temp$pm10.obs)

par(mfrow=c(1,2))
qqnorm(abpm10$pm10.obs, main='ORIGINAL ABPM10', col = 'red')
qqline(abpm10$pm10.obs)

qqnorm(temp$pm10.obs, main='LOG TRANSFORMED ABPM10', col = 'blue')
qqline(temp$pm10.obs)

plot(abpm10.ev, model=abpm10.mv)
par(mfrow=c(1,2))

boxplot(abpm10$pm10.obs, main='ORIGINAL ABPM10')
boxplot(temp$pm10.obs, main='LOG TRANSFORMED ABPM10')

########### Q2 ###################
cut_seq = seq(400, 2000, by = 50)
width_seq = seq(50, 1000, by = 50)
rows = length(cut_seq)*length(width_seq)
evarr = array(list(), c(rows, 6))
i = 1
psill = 150
range = 300
nugget = 100

for (cutoff in  cut_seq){
	for(width in width_seq) { 
		abpm10.ev = variogram(pm10.obs~1, data=abpm10, cutoff=cutoff, width=width)
		abpm10.mv = fit.variogram(abpm10.ev,vgm(psill, "Exp", range, nugget))
		pmk.exp.cv = krige.cv(pm10.obs~1, abpm10, model=abpm10.mv)
		ME.exp.cv = sum(pmk.exp.cv$residual)/length(pmk.exp.cv$residual)
		MSE.exp.cv = sum(pmk.exp.cv$residual^2)/length(pmk.exp.cv$residual)
		evarr[i,] = list(cutoff, width, ME.exp.cv, sqrt(MSE.exp.cv), abpm10.ev, abpm10.mv)
		i = i + 1
	}
}

min_rmse = evarr[[1, 4]]
for (i in 2:rows) {
	v = evarr[[i,5]]
	num_points = length(v[,"np"])
	if (num_points >= 15 && evarr[[i, 4]] < min_rmse) {
		min_variogram = evarr[i,] 			
	}
}
plot(min_variogram[[5]], model=min_variogram[[6]], type='l')

#################### Q3 ####################
cutoff = min_variogram[[1]]  #1050
width = min_variogram[[2]]   #50

models = c('Exp','Sph','Gau','Mat','Ste','Cir','Lin','Bes','Pen','Per','Wav')
results = array(list(), c(length(models),4))
i=1
for (model in models) {
	abpm10_mod.ev = variogram(pm10.obs~1, data=abpm10, cutoff=cutoff, width=width)
	abpm10_mod.mv = fit.variogram(abpm10_mod.ev,vgm(psill, model, range, nugget))
	pmk.mod.cv = krige.cv(pm10.obs~1, abpm10, model=abpm10_mod.mv)
	ME.mod.cv = sum(pmk.mod.cv$residual)/length(pmk.mod.cv$residual)
	MSE.mod.cv = sum(pmk.mod.cv$residual^2)/length(pmk.mod.cv$residual)
	results[i,] = list(ME.mod.cv, sqrt(MSE.mod.cv), abpm10_mod.ev, abpm10_mod.mv)
	i = i + 1
}

#plot(results[[3, 3]], model=results[[2, 4]], type='l')








