/**
 * @yuanhong.song@wsu.edu (Yuanhong Song)
 * adapted from M. He et al, 2018, Remote Sensing
 *
 * linear interpolation of landsat using modis data
 * for temporal resolution
 *
 * General outlite:
 * I. set up working frame
 *    aoi: area of interest
 *    landsat:  landsat series, input should be 7 or 8
 *    year: year of processing
 *    cdlString: cropland data layer, need to modify the year
 *    confBoudn: confidence level low bound of the cropland data
 * II. Interpolation
 *    1. process Landsat
 *    2. process MODIS
 *    3. join Landsat and MODIS series
 *    4. Interpolate
 * III. post-processing
 *    1. mask to cultivated, high CDL confidence araa
 *    2. reproject the data
 *    3. export in batch
 * IV. display data for visual check
 */

/*****************************************************
*
*  Define Area of Interest and Study Time Frame
*
* ****************************************************/

// load AOI (palouse) boundary which saved as GEE assets
var aoi = ee.FeatureCollection("users/yhsong/palouse");
var landsat = 8;
var year = 2015;
var cdlString = 'USDA/NASS/CDL/2015'
var confBound = 85
var myCRS = 'EPSG:32611'

//set time period
var compositestart = ee.Date('2015-01-01');
var compositeend   = ee.Date('2015-12-31');
var cultMask = ee.Image(cdlString).select('cultivated').eq(2);
var confMask = ee.Image(cdlString).select('confidence').gt(confBound);

/** **************************************************
*
*        Set up general work frame
*
* ****************************************************/

// set up AOI boundary for plot -----
var aoiBoundary = ee.Image().toByte().paint(aoi, 3, 5);

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

/*****************************************************
*
*        Landsat (Landsat 8 Surface Reflectance)
*
* ****************************************************/

// CLOUD FILTER --------------
var maskCloudLC = function(img){
   var clear = img.select('pixel_qa').bitwiseAnd(2).neq(0);
   return img.updateMask(clear);
};

// get Landsat data --------------
var bandNames = ee.List (['B2', 'B3', 'B4', 'B5', 'pixel_qa']);
var ls = ee.ImageCollection(lsCollection)
          .select(bandNames)
          .filterDate(startMillis,endMillis)
          .filterBounds(aoi)  // restric the volume of data to process
          .map(function(img) {return img.clip(aoi)})
          .map(maskCloudLC);
print(ls, 'landsat');


// get Landsat NDVI ------------
var getLsNDVI = function(image){
  var doy=ee.Date(image.get('system:time_start')).getRelative('day','year');
  return image.normalizedDifference([nir,red]).rename('lsNDVI')
              .copyProperties(image, [time]);
};

// lsNDVI should be a list of NDVI images from landsat, each with one band
var lsNDVI = ls.map(getLsNDVI);
print(lsNDVI, 'ls ndvi');

/** **************************************************
*
*        MODIS (MOD09Q1 product)
*
* ****************************************************/

// get modis QA bands --------------
var getQABits = function(image, start, end, newName) {
    var pattern = 0;
    for (var i = start; i <= end; i++) {
        pattern += Math.pow(2, i);
       }
    return image.select([0], [newName])
                .bitwiseAnd(pattern)
                .rightShift(start);
};

// MODIS cloud mask function --------------
var maskCloudMod = function(image) {
  var QA = image.select('QA');
  var internalCloud = getQABits(QA, 0, 1, 'MOD_LAND_QA');   // Get the MOD_LAND_QA bits
  return image.mask(internalCloud.eq(0));
};

// get MODIS::MOD09Q1
var mod = ee.ImageCollection('MODIS/MOD09Q1')
            .filterDate(startMillis, endMillis)
            .map(function(img){return img.clip(aoi)})
            .map(maskCloudMod);

print(mod, 'modis surf ref no cloud');

// get MODIS NDVI ----
var getModNDVI = function(image) {
  return image.normalizedDifference(['sur_refl_b02', 'sur_refl_b01'])
              .rename('modNDVI')
              .resample('bilinear')
              .copyProperties(image, [time]);
      //      .copyProperties(image, ['id']);
}

var modNDVI = mod.map(getModNDVI);
print(modNDVI, 'modis ndvi');

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
var joinedNDVI=lsmodis.map(function(img){
  var maxmodis=ee.ImageCollection.fromImages(img.get('ndvi')).select('modNDVI').max();
  return maxmodis.addBands(ee.Image(1))
                 .addBands(img.select('lsNDVI'))
                 .select([0,1,2],['maxmodis','const1','lsNDVI'])
            //     .copyProperties(img,['doy'])
                 .copyProperties(img,[time]);
});
print (joinedNDVI, 'landsat modis ndvi in one collection')

/** **************************************************
*
*        Linear Inperpolation
*
* ****************************************************/

// Landsat = slope * MODIS + intercept
var coeffs = joinedNDVI.reduce(ee.Reducer.linearRegression(2, 1)).select('coefficients');
var slope = coeffs.arrayGet([0,0]);
var intercept = coeffs.arrayGet([1,0]);

// Landsat_hat = slope_hat * MODIS + intercept_hat
var fitted = joinedNDVI.map(function(img){
  var fit1 = img.select(2).add(intercept.add(img.select(0).multiply(slope))).multiply(0.5);
  var fit2 = fit1.unmask(intercept.add(img.select(0).multiply(slope))).rename('fittedNDVI');
  return fit2.copyProperties(img,[time]);
       //      .copyProperties(img,['doy']);
});

var fitNDVI = modNDVI.map(function(img){
  var   result = intercept.add(img.select('modNDVI').multiply(slope)).rename('30mNDVI');
  return result.copyProperties(img,[time])
               .copyProperties(img,['system:index']);
});
print(fitNDVI, 'fitted NDVI');


/*****************************************************
*
*        Filter
*
* ****************************************************/
var cdlFill = function(image){
  var image1 = image.updateMask(cultMask)
  var image2 = image1.updateMask(confMask)
  return image2
};

var filterNDVI = fitNDVI.map(cdlFill);


/** **************************************************
*
*        Project and Export
*
* ****************************************************/

// reproject -------------

var filterNDVI32611 = filterNDVI.map( function(image){
        return image.reproject(myCRS, null, 30);
});

print (filterNDVI32611);

// Export --------------
// https://github.com/fitoprincipe/geetools-code-editor
var col2drive = function(col, folder, region, scale) {
    var n = col.size().getInfo();
    var colList = col.toList(n);
    for (var i = 0; i < n; i++) {
      var img = ee.Image(colList.get(i));
      var id = img.id().getInfo();
      Export.image.toDrive({
        image:img,
        description: id,
        folder: folder,
        fileNamePrefix: id,
        region: region,
        scale: scale});
    }};

col2drive(filterNDVI32611, 'GEE/ndvi', aoi, 30)


/*****************************************************
*
*       Visual Check
*
******************************************************/

// Set up raster display parameters ----
var lsVis = {min: 0,max: 5000, gamma: 1.4,
              bands: ['B4', 'B3', 'B2'],
};
var ndviVis = {min: 0, max: 1,
    palette: ['FF0000', 'FFC300', '6BA44E', '057811', '065292']
};
Map.centerObject(aoi, 9)
Map.addLayer (ls.mean(), lsVis, 'landsat collection');
Map.addLayer (lsNDVI.mean(), ndviVis, 'landsat ndvi mean');
Map.addLayer (modNDVI.mean(), ndviVis, 'modis ndvi mean');
Map.addLayer(fitNDVI.mean(), ndviVis, 'fitted NDVI');
