# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~         
#  LOAD DATA ----
#  analysis and build regression models using NDVI
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(!require(rgdal)){install.packages('rgdal')};library('rgdal')
if(!require(sp)){install.packages('sp')};library('sp')
if(!require(raster)){install.packages('raster')};library('raster')
if(!require(beepr)){install.packages('beepr')};library('beepr')


## function random_points: create random points ------
set.seed(123)

random_points = function(raster, number, xbuf=200, ybuf=200, espg=32611){
  if(!require(sp)){install.packages('sp')};library('sp')
  if(!require(raster)){install.packages('raster')};library('raster')
  
  frame = extent(raster)
  easting = runif(number, frame[1]-xbuf, frame[2]+xbuf)
  northing = runif(number, frame[3]-ybuf, frame[4]+ybuf)
  locID = 1 : number
  points = data.frame(locID, easting, northing)
  coordinates(points) = ~ easting + northing
  proj4string(points) = CRS(paste("+init=epsg:",espg, sep=""))
  print(crs(points))
  return(points)
}
set.seed(123)
ran_pt = random_points(ndvi, 13000, xbuf=50, ybuf = 0, espg=32611)
# ran_pt = read.csv('palousepoints.csv', header = T)
# ran_pt$locID = 1 : nrow(ran_pt)
# plot(ran_pt[, -3])

# CDL ------------
# cropland data layer includes: confidence, crop type, and cultivation
# original -99999 replaced by NA
# cdl_list = list.files('./GEEexpo/CroplandDataLayer/', 
#                       pattern = 'tif$', full.names = T)
# cdl_img = stack(cdl_list)
# crs(cdl_img)
# names(cdl_img)

samp_loc = ran_pt
# coordinates(samp_loc) = ~easting + northing
# samp_loc@proj4string = CRS("+init=epsg:32611")

cdl_val = raster::extract(cdl_img, samp_loc); beep(2)
cdl_val = cbind(samp_loc@data, samp_loc@coords, cdl_val) 
cdl_val[cdl_val==0]=NA
nrow(na.omit(cdl_val))

write.csv(na.omit(cdl_val)[, 1:3], 'sample_location.csv', col.names = T, row.names = F)

plot(na.omit(cdl_val)[, 2:3])
over(na.omit(cdl_val)[, 2:3])

cdl_pal = read.csv(file.choose(), header = T)
cdl_df = base::merge(cdl_df, cdl_pal, by = "Value")


# cdl_cult = raster('./GEEexpo/CDL2017cultivated32611.tif')
# crs(cdl_crop)
# cdl_conf@data@values[cdl_conf@data@values==-99999] = NA
# 
# plot(raster::trim(cdl_conf))
# cdl_conf = raster:: trim(cdl_conf)
# spplot(cdl_conf, col.region=brewer.pal(5, 'Greens'))



## NDVI RAW ---------
# ndvi = raster('ndvi20170709hiConf.tif') 
# ndvi = trim(ndvi)
# ndvi_val = raster::extract(ndvi, ran_pt)

# copy the file path of NDVI images
ndvi_raw_files = list.files(path="./NDVIraw", full.names = T)
head(ndvi_raw_files)
length(ndvi_raw_files)

# exclude out the images with too many missing values
ndvi_raw_files = ndvi_raw_files[sapply(ndvi_raw_files, file.size) > 10000000]
head(ndvi_raw_files)
length(ndvi_raw_files)


# sample ndvi values at random points locations
ndvi_raw_df = list()

for(i in 1:length(ndvi_raw_files)){
  ndvi1 = raster(ndvi_raw_files[i])
  ndvi_raw_df[[i]] = raster::extract(ndvi1, ran_pt) 
}; beep(2)
ndvi_raw_df = data.frame(Reduce(cbind, ndvi_raw_df))
ndvi_raw_df[1:5, 1:8]

if(!require(stringr)){install.packages('stringr')};library('stringr')

ndvi_raw_names = str_remove_all(ndvi_raw_files, "./NDVIraw/LC08_043027_")
ndvi_raw_names = str_remove_all(ndvi_raw_names, ".tif")
ndvi_raw_time = data.frame(date = as.Date(ndvi_raw_names, "%Y%m%d"))
ndvi_raw_time$year = format(ndvi_raw_time$date, "%Y")
ndvi_raw_time$mon = format(ndvi_raw_time$date, "%b")
ndvi_raw_time$day = format(ndvi_raw_time$date, "%d")

names(ndvi_raw_df) = ndvi_raw_names
names(ndvi_raw_df)
ndvi_raw_df = ndvi_raw_df[,1:29]

dat= cbind(ran_pt, ndvi=ndvi_val, cova_val, crop=cdl_val)
dat = na.omit(dat) 
nrow(dat)
names(dat)



## NDVI FUSED ----------



# Plotting rasters ----
if(!require(tidyverse)){install.packages('tidyverse')};library('tidyverse')
library(ggplot2)
library(RColorBrewer)

par <- old.par <- par(no.readonly = T)
plot(cova,8, col = brewer.pal(6, "blues"));points(ran_pt)

ndviCol = palette(c('#F1B555', '#99B718','#74A901', '#529400', '#207401', 
                    '#056201','#004C00', '#023B01'))

pts = dat[,c(1,2)]
coordinates(pts) = ~easting + northing
proj4string(pts) = CRS("+init=epsg:32611")

plot (ndvi,col = 'gray', box=F, axes=F, legend=NULL)
points(pts, pch=19, cex=0.6, col='#e34a33')
title(main='Random sample over AOI, \n Overlapping on Least Cloud Cover Scenes \n (Top 25)', cex=0.8)

plot(dat[,c(3, 12:37)])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~         
#                         Regression Analysis ----
#  
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(!require(caret)){install.packages('caret')};library('caret')
if(!require(leaps)){install.packages('leaps')};library('leaps')
library(MASS)
library(corrplot)
corr_max = cor(dat[,c(3, 12:37)])
sort(corr_max)
corrplot(corr_max, method = "circle", type= "lower", 
         tl.srt=45, tl.col = 'gray20', tl.cex = 1.2)

corr_max2 = cor(dat[,c(3, 15, 16,24, 29, 36)])
corrplot(corr_max2, method = "square", type= "upper", 
         tl.srt=45, tl.col = 'gray20',
         addCoef.col = "black")

corr_mx_sort = corr_max[order(corr_max[,1],decreasing = TRUE),]
row.names(corr_mx_sort)
write.csv(corr_mx_sort, './outputs/ndvi_corr.csv')

## linear regression without crop type ----
length(unique(dat$crop))
dat[!dat$crop %in% c(23, 24, 42, 52, 61), 'crop'] = 300
dat$crop = as.factor(dat$crop)
contrasts(dat$crop) =contr.treatment(6, base = 5)


m1 = lm(ndvi~ ., data = dat[,-c(1,2, 4:8)])
summary(m1)

m2 = stepAIC(m1, direction='both', trace=F)
summary(m2)


# set.seed(2018)
# train.control = trainControl(method = "cv", number = 10)
# m3 = train(formular(m2), data = dat,
#            method = "leapSeq", 
#            tuneGrid = data.frame(nvmax = 1:5),
#            trControl = train.control
# )
# m3$results
# summary(m3$finalModel)



## linear regression with crop type ----
m4 = lm(ndvi~crop, data=dat)
summary(m4)

fitcrop = aov(ndvi ~ crop, data=dat)
summary(fitcrop)
TukeyHSD(fitcrop)


m5 = lm(ndvi ~ . , data = dat[,c(-1, -2, -12)])
summary(m5)

m5.1 = lm(ndvi ~ . + .^2 , data = dat[,c(-1, -2, -12)])
summary(m5.1)

m5.2 = lm(ndvi ~ . + .^2 + .^3 , data = dat[,c(-1, -2, -12)])
summary(m5.2)

# set.seed(2018)
# train.control = trainControl(method = "cv", number = 10)
# m6 = train(formula(m5), 
#            data = dat,
#            method = "leapSeq", 
#            tuneGrid = data.frame(nvmax = 1:10),
#            trControl = train.control
# )
# m6$results
# summary(m6$finalModel)

m7 = stepAIC(m5, direction = 'both', trace=F)
summary(m7)

m7.1 = stepAIC(m5.1, direction = 'both', trace=F)
summary(m7.1)

m7.2 = stepAIC(m5.2, direction = 'both', trace=F)
summary(m7.2)


anova(m7.1, m7.2)

set.seed(2018)
train.control = trainControl(method = "cv", number = 10)
m8.1 = train(formula(m7.1), data = dat,
             method = "leapSeq", # linear reg with stepwise selection
             tuneGrid = data.frame(nvmax = 1:20),
             trControl = train.control
)
m8.1$results
summary(m8.1$finalModel)

trellis.par.set(caretTheme())
ggplot(m8.1)

m8.2 = train(formula(m7.2), data = dat,
             method = "leapSeq", # linear reg with stepwise selection
             tuneGrid = data.frame(nvmax = 1:20),
             trControl = train.control
)
m8.2$results
summary(m8.2$finalModel)

trellis.par.set(caretTheme())
plot(m8.1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#           with historical ndvi----
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

summary(model0 <- lm(ndvi~., data=dat[, -c(1, 2, 11,12:40)]))
summary(model1 <-stepAIC(model0, direction = 'both', trace = F))

summary(model2 <- lm(ndvi~., data=dat[, -c(1, 2, 12:40)]))
summary(model3 <-stepAIC(model2, direction = 'both', trace = F))

summary(model4 <- lm(ndvi~., data=dat[, -c(1, 2)]))
summary(model5 <- stepAIC(model4, direction = 'both', trace = F))

summary(model6 <- lm(ndvi~ crop, data = dat))

summary(model7 <- lm(ndvi~., data = dat[,-c(1,2, 4:11)]))
summary(stepAIC(model7, direction = 'both', trace = F))

summary(model4 <- lm(ndvi ~ asp + elev + twi + crop + 
                       `20140717` + `20160417` + `20160503` + `20140802` + `20150704`, data = dat[, -c(1, 2)]))

summary(model6 <- stepAIC(model4))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~         
#                         Spatial Autocorrelation ----
#  
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
m7re = data.frame(dat[, 1:2], m7$residuals)
names(m7re)[3] = 're'
m7re_sp = m7re
coordinates(m7re_sp) = ~easting+northing
proj4string(m7re_sp) = CRS("+init=epsg:32611")


spplot(m7re_sp, "re", do.log = F,
       key.space=list(x=0,y=1,corner=c(0,1)),
       scales=list(draw=F), cuts = 5,
       col.regions= palette(c('#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6')),
       legendEntries = c("[-0.37,  -0.23]", "[-0.23,  -0.10]", 
                         "[-0.10,   0.03]", "[ 0.03,   0.17]", 
                         "[ 0.17,   0.30]"),
       par.settings=list(panel.background=list(col="grey50")),
)

length(dist(m7re[,1:2])) # number of pairs
dismax=max(dist(m7re[,1:2])) # max distance
min(dist(m7re[,1:2])) # min distance

if(!require(geoR)){install.packages('geoR')}; library(geoR)

# variogram
# cloud
re_vg_cloud = variog(as.geodata((m7re)), option = 'cloud',
                     max.dist  = dismax, estimator.type = "modulus")
plot(re_vg_cloud, pch=19, cex=0.1, 
     col='lightblue', main='Unbinned (averaged) variogram')


# bin
re_vg = variog(as.geodata((m7re)), option = 'bin',
               max.dist  = dismax, estimator.type = "modulus",
               uvec = 20)
plot(re_vg, pch=19, cex=0.1, 
     col='lightblue', main='Binned variogram')
max(m7re$re)

# Spatial ACF test by MORAN's I
if(!require(ape)){install.packages('ape')}; library(ape)
dist=as.matrix(dist(as.matrix(m7re[,1:2])))
dist.inv = log(1/dist)
diag(dist.inv) = 0
dist.inv[1:5, 1:5]

Moran.I(m7re$re, dist.inv, alternative = 'greater')

summary(as.vector(1/log(dist)))

ndvi_sp = dat[,1:3]
coordinates(ndvi_sp) = ~easting+northing
proj4string(ndvi_sp) = CRS("+init=epsg:32611")
ndvisvg = variogram(ndvi~1, data=ndvi_sp,
                    cutoff = dismax,
                    cressie=T)
plot(ndvisvg, pch=19, cex=1.2,
     main='varigram of NDVI')

## ..vrg signif ----         
vrgEnv = function(data, xcol=NA, ycol=NA, zcol=NA, nbin=20, max.dist=NA, dist.rate=1,nsim=99){
  # vrgEnv is for ploting variogram and variogram envelope
  # ~~~~~~ convert data type if necessary ~~~~~~~~~~~~~~~~
  if(!require(geoR)){install.packages('geoR')}; library(geoR)
  if(!class(data)=='geodata'){
    geodata = as.geodata(data, coords.col = c(xcol, ycol), data.col=zcol)
  }else{
    geodata = data
  }
  if(is.na(max.dist)){
    max.dist = max(dist(geodata$coords))
  }
  vg = variog(geodata, max.dist = dist.rate*max.dist, uvec=nbin)
  env = variog.mc.env(geodata, obj.var=vg, nsim=nsim)
  plot(vg, envelope=env, pch=19, cex=0.3,
       main=paste("Envelops for Empirical Variograms", nsim))
}

# png('./outputs/cor_elev.png')
vrgEnv(dat[,c(1, 2, 6)], xcol=1, ycol=2, zcol=3, nbin=100, nsim=99)
# dev.off()




library(ggplot2)

theme=theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.background = element_rect(fill = 'grey95'),
            axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),
            legend.position='none')

p1 = ggplot(dat, aes(easting, northing))+
  geom_point(aes(color=ndvi), size=3)+
  scale_color_gradient2(low = "#E53935", mid='#78909C',high = '#00CC66',  midpoint = 0.55)+
  theme+
  ggtitle('NDVI')


p2 = ggplot(dat, aes(easting, northing))+
  geom_point(aes(color=ppt), size=3)+
  scale_color_gradient2(low = "#E53935", mid='#78909C',high = '#00CC66',  midpoint = mean(dat$ppt))+
  theme+
  ggtitle('PPT')

p3 = ggplot(dat, aes(easting, northing))+
  geom_point(aes(color=elev), size=3)+
  scale_color_gradient2(low = "#E53935", mid='#78909C',high = '#00CC66',  midpoint = mean(dat$elev))+
  theme+
  ggtitle('ELEV')

p4 = ggplot(dat, aes(easting, northing))+
  geom_point(aes(color=tem), size=3)+
  scale_color_gradient2(low = "#E53935", mid='#78909C',high = '#00CC66',  midpoint = mean(dat$tem))+
  theme+
  ggtitle('TEMP')

library('gridExtra')
png('./outputs/trend.png', width = 800, height = 300)
grid.arrange(p1, p3, p2, p4, nrow=1)
dev.off()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#           lme ----
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(!require(ncf)){install.packages('ncf')};library(ncf)
if(!require(nlme)){install.packages('nlme')}; library(nlme)
library(lattice)

sac_ndvi = spline.correlog(x=dat$easting, y=dat$northing, z=dat$ndvi, resamp = 100)
png("./outputs/SpatCorr_ndvi.png")
plot(sac_ndvi, main="NDVI Spatial Autocorrelation")
dev.off()

sac_res = spline.correlog(x=dat$easting, y=dat$northing, z=resid(m7), resamp = 100)
png("./outputs/SpatCorr_resi.png")
plot(sac_res, main="Residual Spatial Autocorrelation")
dev.off()

dat$grp = rep(1, dim(dat)[1])
mm1 = lme(formula(m7), data=dat, 
          random = ~1|grp,
          corr = corSpatial(form=~easting+northing, 
                            type='exponential',
                            nugget=F),
          method='ML')

mm2 = lme(formula(m7), data=dat, 
          random = ~1|grp,
          corr = corSpatial(form=~easting+northing, 
                            type='exponential',
                            nugget=T),
          method='ML')


# sink('./outputs/sum_mixed1.txt')
# print(summary(mm1))
# sink()
# 
# 
# sink('./outputs/sum_ols7.txt')
# print(summary(m7))
# print(paste("AIC",AIC(m7), sep = " "))
# sink()





