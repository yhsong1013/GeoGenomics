var palettes = require('users/gena/packages:colorbrewer').Palettes;
var utils = require('users/gena/packages:utils')
var text = require('users/gena/packages:text')

////////////////////////////////////////////////////////////////////////////////////////////
////////Code completed by Erin Trochim, UAF postdoc, for Tessa Hasbrouck, UAF, M.S. Thesis
////////////////////////////////////////////////////////////////////////////////////////////

//Creates start day to be Aug 20 and end day to be Oct 1
var startDay = 232; //use noaa day of year calendar
var endDay = 274;

//Brings in Global Surface Water data
var gsw = ee.Image('JRC/GSW1_0/GlobalSurfaceWater');
var recurrence = gsw.select('recurrence');

//Masks water bodies using GlobalSurface Water data
var maskWater = function(image) {
  var water = recurrence.gte(90).unmask(0).eq(0);
  return image.updateMask(water);
};

//Creates polygon for entire area around Nulato
var sites = new ee.FeatureCollection([
    ee.Feature(
        ee.Geometry.Polygon(
            [[-159.158935546875,65.02485746814963],[-159.19189453125,64.31039624941911],
[-156.4013671875,64.31039624941911],[-156.434326171875,65.06657315085853],[-159.158935546875,65.02485746814963]]),
        {name: 'StudyArea', fill: 2})
    ]);

//Brings in MODIS 8-Day composite data from 2001-2016
var modis_filtered = ee.ImageCollection('MODIS/MOD09Q1')
    .filterDate('2000-01-01', '2017-12-31')
   .filter(ee.Filter.dayOfYear(startDay, endDay))
   .map(maskWater)
   .filterBounds(sites);

//Calculates NDVI by week
var keepProperties = ['system:time_start', 'system:time_end'];
var getndvi = function(image) {
  return image.normalizedDifference(['sur_refl_b02', 'sur_refl_b01']).copyProperties(image,keepProperties)
};

//Creates an image collection from NDVI values
var ndvi = ee.ImageCollection(modis_filtered.map(getndvi));
print(ndvi);

//Sets the start and end year
var startyear = 2000;
var endyear = 2017;

//Create visualization palettes of NDVI
print(palettes); //to see different palette options

var bounds = sites.geometry().bounds()

//make the data 8-bit which is necessary for making a video
var ndvi_video =  ndvi.map(function(image){
  var start = ee.Date(image.get('system:time_start'))
  var end = ee.Date(image.get('system:time_end'))
  var label = start.format('YYYY-MM-dd').cat(' - ').cat(end.format('YYYY-MM-dd'))

  return image.visualize({
    forceRgbOutput: true,
    palette: palettes.BrBG[9], //original palette was "000000", "fdbb84"
    min: 0,
    max: 1
  }).clip(bounds).set({label: label});
});


// annotate
var annotations = [
  {
    position: 'left', offset: '1%', margin: '1%', property: 'label', scale: Map.getScale() * 2
  }
]

ndvi_video = ndvi_video.map(function(image) {
  return text.annotateImage(image, {}, bounds, annotations)
})



// add a few layers to map
var animation = require('users/gena/packages:animation')
animation.animate(ndvi_video, {maxFrames: 5})

Map.centerObject(bounds, 9)

//Export NDVI from whole study area to video
Export.video.toDrive({
  collection: ndvi_video,
  description: "NDVItimelapse",    // Filename, no spaces allowed
  framesPerSecond: 1,             // I.e., 1 year / second
  region: bounds,
  scale: 250,                     // Scale in m
  });
