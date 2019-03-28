/**
 * Making CDL time lapse video
 * @ yuanhong.song@wsu.edu (Yuanhong Song)
 *
 * Adapted from E. Trochim:
 * https://code.earthengine.google.com/6270df443326ec0d90a18838bd91c5a5
 *
 * v3. explicit CDL corlor palette, 255 values
 */

// DEFINE AREA OF INTEREST ----------------
var polygon = ee.Geometry.Rectangle([-117.8101, 46.6028, -116.8488, 47.4231]);
var aoi = ee.FeatureCollection(polygon);
var bounds = aoi.geometry().bounds();

// WATERBODY MASK -------------------------
// Global Surface Water data for constant water coverage area
var gsw = ee.Image('JRC/GSW1_0/GlobalSurfaceWater');
var recurrence = gsw.select('recurrence');

var maskWater = function(image) {
  var water = recurrence.gte(90).unmask(0).eq(0);
  return image.updateMask(water);
};


// get cropland collection -----------
var cdl = ee.ImageCollection('USDA/NASS/CDL').filterBounds(aoi)
              .filterDate('2008', '2018')
              .select('cropland')
              .map(maskWater)
              .filterBounds(aoi)
              // return image.copyProperties(image, ['system:time_start'])});
print(cdl);

// CREATING COLOR PALETTE ---------------
// the following is a explicit color palette containing 255 colors
// extracted from GEE cropland stantdard palette, and NA color values are
// filled with white ("ffffff").
// The palette is necessary to turn the image collection to 8-bit video correctly
// Haven't found a more efficient way yet.
var cdl_palette = [
  "ffff00", "ffd300", "ff2626", "00a8e2", "ff9e0a", "267000", "ffff00",
  "ffffff", "ffffff", "ffffff", "70a500", "00af49", "dda50a", "dda50a",
  "7cd3ff", "ffffff", "ffffff", "ffffff", "ffffff", "ffffff", "ffffff",
  "e2007c", "896054", "d8b56b", "a57000", "d69ebc", "707000", "aa007c",
  "a05989", "700049", "d69ebc", "d1ff00", "7c99ff", "d6d600", "d1ff00",
  "00af49", "ffa5e2", "a5f28c", "00af49", "d69ebc", "ffffff", "a800e2",
  "a50000", "702600", "00af49", "af7cff", "702600", "ff6666", "ff6666",
  "ffcc66", "ff6666", "00af49", "00ddaf", "54ff00", "f2a377", "ff6666",
  "00af49", "7cd3ff", "e8bfff", "afffdd", "00af49", "bfbf77", "ffffff",
  "93cc93", "c6d69e", "ccbfa3", "ff00ff", "ff8eaa", "ba004f", "704489",
  "007777", "af9970", "ffff7c", "ffffff", "b5705b", "00a582", "e8d6af",
  "af9970", "ffffff", "ffffff", "ffffff", "f2f2f2", "999999", "4970a3",
  "ffffff", "ffffff", "ffffff", "7cafaf", "e8ffbf", "ffffff", "ffffff",
  "ffffff", "00ffff", "ffffff", "ffffff", "ffffff", "ffffff", "ffffff",
  "ffffff", "ffffff", "ffffff", "ffffff", "ffffff", "ffffff", "ffffff",
  "ffffff", "ffffff", "ffffff", "ffffff", "ffffff", "ffffff", "4970a3",
  "d3e2f9", "ffffff", "ffffff", "ffffff", "ffffff", "ffffff", "ffffff",
  "ffffff", "ffffff", "999999", "999999", "999999", "999999", "ffffff",
  "ffffff", "ffffff", "ffffff", "ffffff", "ffffff", "ccbfa3", "ffffff",
  "ffffff", "ffffff", "ffffff", "ffffff", "ffffff", "ffffff", "ffffff",
  "ffffff", "93cc93", "93cc93", "93cc93", "ffffff", "ffffff", "ffffff",
  "ffffff", "ffffff", "ffffff", "ffffff", "ffffff", "c6d69e", "ffffff",
  "ffffff", "ffffff", "ffffff", "ffffff", "ffffff", "ffffff", "ffffff",
  "ffffff", "ffffff", "ffffff", "ffffff", "ffffff", "ffffff", "ffffff",
  "ffffff", "ffffff", "ffffff", "ffffff", "ffffff", "ffffff", "ffffff",
  "ffffff", "e8ffbf", "ffffff", "ffffff", "ffffff", "ffffff", "ffffff",
  "ffffff", "ffffff", "ffffff", "ffffff", "ffffff", "ffffff", "ffffff",
  "ffffff", "7cafaf", "ffffff", "ffffff", "ffffff", "ffffff", "7cafaf",
  "ffffff", "ffffff", "ffffff", "ffffff", "ffffff", "ffffff", "ffffff",
  "ffffff", "00ff8c", "d69ebc", "ff6666", "ff6666", "ff6666", "ff6666",
  "ff8eaa", "334933", "e27026", "ff6666", "ff6666", "ffffff", "ff6666",
  "af9970", "ff8eaa", "ff6666", "ff8eaa", "ff6666", "ff6666", "ff8eaa",
  "00af49", "ffd300", "ffd300", "ff6666", "ffffff", "ff6666", "896054",
  "ff6666", "ff2626", "e2007c", "ff9e0a", "ff9e0a", "a57000", "ffd300",
  "a57000", "267000", "267000", "ffd300", "000099", "ff6666", "ff6666",
  "ff6666", "ff6666", "ff6666", "ff6666", "ff6666", "ff6666", "ffffff",
  "ffffff", "ffffff", "267000"
  ];
print(cdl_palette);

// extract time footprint in image -----------
var cdlVideo = cdl.map(function(image){
  var label = ee.Date(image.get('system:time_start')).format('YYYY')

  return image.visualize({
              forceRgbOutput: true,
              palette: cdl_palette,
              min: 0,
              max: 254
    }).clip(aoi).set({label: label});
});


// add time footprint as annotation in each image -----------
var annotations = [{
   position: 'left', offset: '1%', margin: '1%',
   property: 'label', scale: Map.getScale() * 2
  }];

// turn images to vedio ----------------------
// meanwhile add annotation to each image
var text = require('users/gena/packages:text')
var cropVideo = cdlVideo.map(function(image) {
  return text.annotateImage(image, {}, bounds, annotations)
})

// show instant video in GEE map (optional) --------------
Map.centerObject(aoi, 9);
var animation = require('users/gena/packages:animation')
animation.animate(cropVideo, {maxFrames: 50})

// Export video to drive
Export.video.toDrive({
  folder:'GEE',
  collection: cropVideo,
  description: "cdlTimeLapse_255colors",
  region: bounds,
  framesPerSecond: 2,
  scale:30
  });
