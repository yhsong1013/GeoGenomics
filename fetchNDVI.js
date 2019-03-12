/*
Fetching NDVI at AOI  Y.song 11/2018
*/


// ----------------------------------------------------------------------------
//        Define the feature class boundaries
// ----------------------------------------------------------------------------

// 'palouse' is the user define Palouse boundary
var palouse = ee.Geometry.Polygon(
        [
          [
            [-117.84445422977853, 46.729012061434204],
            [-117.39676135868478, 46.73277721796108],
            [-117.00949329227853, 46.50829288685803],
            [-116.83371204227853, 46.66496418338115],
            [-117.00949329227853, 46.794864354166634],
            [-117.08090442509103, 47.428527677686745],
            [-117.36929553837228, 47.430385733400996]
          ]
        ]);

Map.addLayer(palouse, 9);

// ----------------------------------------------------------------------------
//        Cloud Filter
// ----------------------------------------------------------------------------

// 'maskL8sr' is a function to remvoe cloud and cloud shadow
function maskL8sr(image) {
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
  var cloudShadowBitMask = (1 << 3);
  var cloudsBitMask = (1 << 5);
  // Get the pixel QA band.
  var qa = image.select('pixel_qa');
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                 .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return image.updateMask(mask);
}

// ----------------------------------------------------------------------------
//        Getting CDL for data masking
// ----------------------------------------------------------------------------

// load cropland data layer, use it to exclude non-cultivar pixels
// the study area will be where landsat content is, so will first mask with the landsat boundary

function maskcdl(image){
      var mask = cdl.select('cultivated').eq(2)
                    .and(cdl.select('confidence').gt(86));
      return image.updateMask(mask);
}

// ----------------------------------------------------------------------------------------
//        Function for NDVI
// ----------------------------------------------------------------------------------------

var getNDVI = function(image){
  return image.addBands(image.normalizedDifference(['B5','B4'])
                        .rename('ndvi'));
};





// ----------------------------------------------------------------------------------------
//        Collection to Drive
// ----------------------------------------------------------------------------------------

// function for saving the image collection to Google Drive, modified from:
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
        scale: scale})
    }
  }


// ----------------------------------------------------------------------------------------
//        Getting data for one scene
// ----------------------------------------------------------------------------------------

var cdl = ee.Image('USDA/NASS/CDL/2017').clip(palouse);

var startDate = '2018-01-01';
var endDate = '2018-12-31';

var LC = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
            .filterDate(startDate, endDate)
            .filterBounds(ee.Geometry.Point(-117.2218, 47.0351))
            .map(maskL8sr)
            .map(maskcdl)
            .map(getNDVI)
            .select('ndvi');

col2drive(LC, 'LC2018', palouse, 30)


var ndviVis = {min: -1, max: 1,
    palette: ['FFFFFF', 'CE7E45', 'DF923D', 'F1B555', 'FCD163', '99B718',
               '74A901', '66A000', '529400', '3E8601', '207401', '056201',
               '004C00', '023B01', '012E01', '011D01', '011301']};


var coldisp = function(col, vis) {
    var n = col.size().getInfo();
    var colList = col.toList(n);

    for (var i = 0; i < n; i++) {
      var img = ee.Image(colList.get(i));
      var id = img.id().getInfo();

      Map.addLayer(img, vis, id);
    }
  }

Map.centerObject(palouse,9);
//coldisp(LC, ndviVis);

//var lsatVis = {bands: ['B4', 'B3', 'B2'],min: 0,max: 5000,gamma: 1.4,};
//var trueC = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
//            .filterDate(startDate, endDate)
//            .filterBounds(ee.Geometry.Point(-117.2218, 47.0351))
//            .map(maskcdl);
//coldisp(trueC, lsatVis)

//var trueCnoC = trueC.map(maskL8sr);
//coldisp(trueCnoC, lsatVis)
