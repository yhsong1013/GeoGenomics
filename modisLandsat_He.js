//import the fusion table and focus on the state of Montana
// the fustion table id could be found by import fusion table
// masked by seven crop types in MT.

var states=ee.FeatureCollection("ft:1FTW4Q7Es7jbQO7TiOlXylLBJh7YInwEYHtV4ro4","geometry");
var MTstate=states.filterMetadata("id","equals","MT");

var MTstateBoundary = ee.Image().toByte().paint(MTstate, 3, 5);        // Outline using color 3, width 5.

Map.addLayer(MTstateBoundary, {
    palette: '000000,FF0000,00FF00,0000FF',
    max: 3,
    opacity: 0.5
});

//set time period
var compositestart=ee.Date('2008-04-01');
var compositeend=ee.Date('2008-10-01');

//set bands to retrieve NDVI and pixel_qa to select clear pixels;
var bandNames=ee.List(['B3','B4','pixel_qa']);
var time='system:time_start';

var tempwin=8;  //set time window

//get image from imagecollection by start and end date;
//convert to milliseconds;
var startMillis=compositestart.millis();
var endMillis=compositeend.millis();

//get Landsat 5 or 7 SR data
var ls=ee.ImageCollection('LANDSAT/LT05/C01/T1_SR')
//var ls=ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
          .select(bandNames)
          .filterDate(startMillis,endMillis)
          .filterBounds(MTstate);

//create function to select only clear pixels in the images
var maskClouds=function(img){
  var clear=img.select('pixel_qa').bitwiseAnd(2).neq(0);
  return img.updateMask(clear);
};


//apply the above functions to image collection
var lscloudmask=ls.map(maskClouds);
var lsndvi=lscloudmask.map(function(img){
  var doy=ee.Date(img.get('system:time_start')).getRelative('day','year');
  return img.normalizedDifference(['B4','B3']).rename('lsNDVI')
            .copyProperties(img,[time]);
      //      .set({'doy':doy});
});

lsndvi=lsndvi.map(function(img){return img.clip(MTstate)});

//MODIS surface reflectance
//must use clip to get images for MT since MOD09Q1 is global mosaic already.
var modsr=ee.ImageCollection('MODIS/MOD09Q1')
            .filterDate(startMillis,endMillis)
            .map(function(img){return img.clip(MTstate)});
//print(modsr);

var getQABits = function(image, start, end, newName) {
    // Compute the bits we need to extract.
    var pattern = 0;
    for (var i = start; i <= end; i++) {
       pattern += Math.pow(2, i);
    }
    // Return a single band image of the extracted QA bits, giving the band
    // a new name.
    return image.select([0], [newName])
                  .bitwiseAnd(pattern)
                  .rightShift(start);
};

// A function to mask out cloudy pixels.
var maskClouds = function(image) {
  // Select the QA band.
  var QA = image.select('QA');
  // Get the MOD_LAND_QA bits
  var internalCloud = getQABits(QA, 0, 1, 'MOD_LAND_QA');
  // Return an image masking out cloudy areas.
  return image.mask(internalCloud.eq(0));
};

var modcloudmask=modsr.map(maskClouds);
var modndvi=modcloudmask.map(function(img){
  return img.normalizedDifference(['sur_refl_b02','sur_refl_b01']).rename('modNDVI')
            .resample('bilinear')
            .copyProperties(img,[time]);
      //      .copyProperties(img,['id']);
});

var tempdif=tempwin * 24 * 60 * 60 * 1000;

//define the join
var join=ee.Join.saveAll({'matchesKey':'ndvi','ordering':time,'measureKey':'timeDiff'});
var maxDiff=ee.Filter.maxDifference({
  'difference':tempdif,
  'leftField':time,
  'rightField':time
});
// Define the prior filter
var prior=ee.Filter.greaterThan({'leftField':time,'rightField':time});
var filter1=ee.Filter.and(maxDiff, prior);
//var filter1=max.add(ee.Filter.greaterThan({'leftField':time,'rightField':time}));
//define the post filter
var post=ee.Filter.lessThan({'leftField':time,'rightField':time});
var filter2=ee.Filter.and(maxDiff,post);
//var filter2=max.add(ee.Filter.LessThan({'leftField':time,'rightField':time}));
var filter=ee.Filter.or(filter1,filter2);
//apply join
var lsmodis=ee.ImageCollection(join.apply(lsndvi,modndvi,filter));
//print(lsmodis)

//create images with three bands(modNDVI,constant value of 1,lsNDVI)
var ndvi=lsmodis.map(function(img){
  var maxmodis=ee.ImageCollection.fromImages(img.get('ndvi'))
                                  .select('modNDVI')
                                  .max();
  return maxmodis.addBands(ee.Image(1))
                 .addBands(img.select('lsNDVI'))
                 .select([0,1,2],['maxmodis','const1','lsNDVI'])
            //     .copyProperties(img,['doy'])
                 .copyProperties(img,[time]);
});

//linear regression between MODIS and Landsat
//Landsat=slope*MODIS+intercept
var coeffs=ndvi.reduce(ee.Reducer.linearRegression(2,1)).select('coefficients');
var slope=coeffs.arrayGet([0,0]);
var intercept=coeffs.arrayGet([1,0]);

//using slope and intercept to linear interpolate the images
var fitted=ndvi.map(function(img){
  //using the average value of original lsNDVI and fitted NDVI, kind of smooth the data, remove short variations
  //do not compute the fitted NDVI for gaps
  var fit1=img.select(2).add(intercept.add(img.select(0).multiply(slope))).multiply(0.5);
  //filling the gaps using unmask
  var fit2=fit1.unmask(intercept.add(img.select(0).multiply(slope))).rename('fittedNDVI');
  return fit2.copyProperties(img,[time]);
       //      .copyProperties(img,['doy']);
});

//compute 30m NDVI using MODIS and slope,intercept derived above.
var fittedNDVI=modndvi.map(function(img){
  var result=intercept.add(img.select('modNDVI').multiply(slope)).rename('30mNDVI');
  return result.copyProperties(img,[time])
               .copyProperties(img,['system:index']);
});
print(fittedNDVI);
