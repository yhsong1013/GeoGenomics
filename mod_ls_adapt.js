/**
 * @yuanhong.song@wsu.edu (Yuanhong Song)
 * adapted from M. He et al, 2018, Remote Sensing
 *
 * linear interpolation of landsat using modis data
 * for temporal resolution
 */

/*****************************************************
*
*  Define Area of Interest and Study Time Frame
*
* ****************************************************/

// load AOI (palouse) boundary which saved as GEE assets
var aoi = ee.FeatureCollection("users/yhsong/palouse");

var landsat = 8;

//set time period
var compositestart = ee.Date('2013-04-01');
var compositeend   = ee.Date('2013-10-01');


/** **************************************************
*
*        Set up general work frame
*
* ****************************************************/

// set up AOI boundary for plot -----
var aoiBoundary = ee.Image().toByte().paint(aoi, 3, 5);
Map.centerObject (aoi,9);
Map.addLayer (aoiBoundary, {
    palette: '000000,FF0000,00FF00,0000FF',
    max: 3,
    opacity: 0.5
});

// set up flow control for landsat collection ----
if (landsat ==7) {
  var lsCollection = 'LANDSAT/LE07/C01/T1_SR';
  var nir = 'B4';
  var red = 'B3';
} else if (landsat == 8 ) {
  var lsCollection = 'LANDSAT/LC08/C01/T1_SR';
  var nir = 'B5';
  var red = 'B4' ;
} else {
  print ('only support Landsat7 and Landsat8');
}

// convert time format to milliseconds -----
var startMillis = compositestart.millis();
var endMillis   = compositeend.millis();
var time = 'system:time_start';

// Set up raster display parameters ----
var lsVis = {min: 0,max: 5000, gamma: 1.4,
              bands: ['B4', 'B3', 'B2'],
};
var ndviVis = {min: 0, max: 1,
    palette: ['FF0000', 'FFC300', '6BA44E', '057811', '065292']
};

// use Cropland Data Layer to to exclude non-cultivar pixels
var cdl = ee.Image('USDA/NASS/CDL/2017').clip(aoi);
function maskCDL(image){
      var mask = cdl.select('cultivated').eq(2)
                    .and(cdl.select('confidence').gt(86));
      return image.updateMask(mask);
}

/** **************************************************
*
*        Landsat (Landsat 8 Surface Reflectance)
*
* ****************************************************/

/**
 * CLOUD FILTER --------------
 * 'maskCloudLC' is a function to remvoe cloud and cloud shadow
 * cited from M. He et al, 2018, Remote Sensing
 */
var maskCloudLC = function(img){
   var clear = img.select('pixel_qa').bitwiseAnd(2).neq(0);
   return img.updateMask(clear);
};

// get Landsat data
var bandNames = ee.List (['B2', 'B3', 'B4', 'B5', 'pixel_qa']);
var ls = ee.ImageCollection(lsCollection)
          .select(bandNames)
          .filterDate(startMillis,endMillis)
          .filterBounds(aoi)  // restric the volume of data to process
          .map(function(img) {return img.clip(aoi)})
          .map(maskCloudLC);
print(ls, 'landsat');
Map.addLayer (ls.mean(), lsVis, 'landsat collection');


var getLsNDVI = function(image){
  var doy=ee.Date(image.get('system:time_start')).getRelative('day','year');
  return image.normalizedDifference([nir,red]).rename('lsNDVI')
              .copyProperties(image, [time]);
};

// lsNDVI should be a list of NDVI images from landsat, each with one band
var lsNDVI = ls.map(getLsNDVI);
print(lsNDVI, 'ls ndvi');
Map.addLayer (lsNDVI.mean(), ndviVis, 'landsat ndvi mean');

/** **************************************************
*
*        MODIS (MOD09Q1 product)
*
* ****************************************************/

/**
 * get modis QA bands --------------
 * the function compute the bits needed to be extracted
 * returns a single band image of the extracted QA bits
 * then gives the extracted band a new namge
 */
var getQABits = function(image, start, end, newName) {
    var pattern = 0;
    for (var i = start; i <= end; i++) {
        pattern += Math.pow(2, i);
       }
    return image.select([0], [newName])
                .bitwiseAnd(pattern)
                .rightShift(start);
};

/**
 * MODIS cloud mask function --------------
 * based on the MODIS QA bands and the getQABits function
 * the function returns an image masking out cloudy areas.
 */
var maskCloudMod = function(image) {
  var QA = image.select('QA');
  var internalCloud = getQABits(QA, 0, 1, 'MOD_LAND_QA');   // Get the MOD_LAND_QA bits
  return image.mask(internalCloud.eq(0));
};

//must use clip to get images for MT since MOD09Q1 is global mosaic already.
var mod = ee.ImageCollection('MODIS/MOD09Q1')
            .filterDate(startMillis, endMillis)
            .map(function(img){return img.clip(aoi)})
            .map(maskCloudMod);

print(mod, 'modis surf ref no cloud');

var getModNDVI = function(image) {
  return image.normalizedDifference(['sur_refl_b02', 'sur_refl_b01'])
              .rename('modNDVI')
              .resample('bilinear')
              .copyProperties(image, [time]);
      //      .copyProperties(image, ['id']);
}

/**
 * modNDVI should be a list of NDVI images from MODIS
 * each image is one band only
 * there should be more MODIS NDVI than Landsat NDVI
 */
var modNDVI = mod.map(getModNDVI);
print(modNDVI, 'modis ndvi');
Map.addLayer (modNDVI.mean(), ndviVis, 'modis ndvi mean');


/** **************************************************
*
*        JOIN
*
* ****************************************************/

//set time window, time difference for the join
var tempwin = 8;  // 8 days
var tempdif = tempwin * 24 * 60 * 60 * 1000;

// define the join
var join = ee.Join.saveAll({
  'matchesKey': 'ndvi', 'ordering': time, 'measureKey': 'timeDiff'
});

var maxDiff = ee.Filter.maxDifference({
  'difference': tempdif, 'leftField': time, 'rightField': time
});

// Define the prior filter ------
var prior = ee.Filter.greaterThan({'leftField': time, 'rightField': time});
var filter1 = ee.Filter.and(maxDiff, prior);
//var filter1 = max.add(ee.Filter.greaterThan({'leftField':time,'rightField':time}));

// define the post filter ----------
var post = ee.Filter.lessThan({
  'leftField': time, 'rightField': time});
var filter2 = ee.Filter.and(maxDiff,post);
//var filter2=max.add(ee.Filter.LessThan({'leftField':time,'rightField':time}));

var filter = ee.Filter.or(filter1, filter2);

// apply join
var lsmodis=ee.ImageCollection(join.apply(lsNDVI, modNDVI, filter));
print(lsmodis, 'landsat modis ndvi join')

//create images with three bands(modNDVI,constant value of 1, lsNDVI)
var ndvi=lsmodis.map(function(img){
  var maxmodis=ee.ImageCollection.fromImages(img.get('ndvi')).select('modNDVI').max();
  return maxmodis.addBands(ee.Image(1))
                 .addBands(img.select('lsNDVI'))
                 .select([0,1,2],['maxmodis','const1','lsNDVI'])
            //     .copyProperties(img,['doy'])
                 .copyProperties(img,[time]);
});
print (ndvi, 'landsat modis ndvi in one collection')

/** **************************************************
*
*        Linear Inperpolation
*
* ****************************************************/

//linear regression between MODIS and Landsat
//Landsat = slope * MODIS + intercept
var coeffs = ndvi.reduce(ee.Reducer.linearRegression(2, 1)).select('coefficients');
var slope = coeffs.arrayGet([0,0]);
var intercept = coeffs.arrayGet([1,0]);

//using slope and intercept to linear interpolate the images
var fitted=ndvi.map(function(img){
  // using the average value of original lsNDVI and fitted NDVI, kind of smooth the data,
  // remove short variations do not compute the fitted NDVI for gaps
  var fit1=img.select(2).add(intercept.add(img.select(0).multiply(slope))).multiply(0.5);
  //filling the gaps using unmask
  var fit2=fit1.unmask(intercept.add(img.select(0).multiply(slope))).rename('fittedNDVI');
  return fit2.copyProperties(img,[time]);
       //      .copyProperties(img,['doy']);
});

//compute 30m NDVI using MODIS and slope,intercept derived above.
var fittedNDVI=modNDVI.map(function(img){
  var result=intercept.add(img.select('modNDVI').multiply(slope)).rename('30mNDVI');
  return result.copyProperties(img,[time])
               .copyProperties(img,['system:index']);
});
print(fittedNDVI, 'fitted NDVI');
Map.addLayer(fittedNDVI.mean(), ndviVis, 'fitted NDVI');
