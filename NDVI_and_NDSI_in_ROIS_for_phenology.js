// Link https://code.earthengine.google.com/edacb1f4b6f43d2ac8888aed0705d24f 

// var roi = DOCwatershed;
var rois = ROIS;

Map.addLayer(rois);
// ~~~~~~~~~~~~~~~~~~~~~~~~
//    MODIS snow cover
// ~~~~~~~~~~~~~~~~~~~~~~~~
var start = ('2000-01-01')
var end = ('2019-08-31')
// Mask 
var getQABits = function(image, start, end, newName) {
    // Compute the bits we need to extract
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
// ~~~~~~~~ QA NDSI~~~~~~~~ 
var filterBadObsNDSI = function(image){
    var QA = getQABits(image.select('NDSI_Snow_Cover_Basic_QA'),0,15,'atm_corr');
    var maskedImage = ee.Image(0).where(
        QA.eq(0)                           // Only corrected atmosphere
        .or(QA.eq(1))                 // No clouds
        .or(QA.eq(2))                    // Clear or unset 
    ,1);                                      //Put a 1 on good pixels
    return image.mask(maskedImage);            //Return only good pixels
};

var MaskMODISSnow = function(image) {
  var c = filterBadObsNDSI(image);
  return image.updateMask(c);
};

// ~~~~~~~~ QA NDVI ~~~~~~~~ 
var filterBadObsNDVI = function(image){
    var AtmQA = getQABits(image.select('QA'),12,12,'atm_corr');
    var ModlandQA = getQABits(image.select('QA'),0,1,'cloud');
    var dataQA = getQABits(image.select('QA'),4,11,'data_quality');
    var maskedImage = ee.Image(0).where(
        AtmQA.eq(1)                           // Only corrected atmosphere
        .and(ModlandQA.eq(0))                 // No clouds
        .and(dataQA.eq(0))                    // Clear or unset 
    ,1);                                      //Put a 1 on good pixels
    return image.mask(maskedImage);            //Return only good pixels
};

var MaskNDVI = function(image) {
  var c = filterBadObsNDVI(image);
  return image.updateMask(c);
};

// LST 
var filterBadObslst = function(image){
    var LSTerrorQA = getQABits(image.select('QC_Day'),6,7,'LST_error');
    var EmmisiQA = getQABits(image.select('QC_Day'),4,5,'Emmisivity_error');
    var MandatoryQA = getQABits(image.select('QC_Day'),0,1,'Mandatory_QA');
    var dataQA = getQABits(image.select('QC_Day'),2,3,'data_qualityQA');
    var maskedImage = ee.Image(0).where(
        EmmisiQA.eq(0).or(EmmisiQA.eq(1))       // Only corrected atmosphere
        .and(LSTerrorQA.eq(0)).or(LSTerrorQA.eq(1))//.or(LSTerrorQA.eq(2))    
        .and(MandatoryQA.eq(0))                 // No clouds
        .and(dataQA.eq(0))                    // Clear or unset 
    ,1);                                      //Put a 1 on good pixels
    return image.mask(maskedImage);            //Return only good pixels
};

var MaskLST = function(image) {
  var b = filterBadObslst(image);
  return image.updateMask(b);
};


// ~~~~~~~~ NDSI ~~~~~~~
// var NDSICollection = ee.ImageCollection('MODIS/006/MOD10A1')
//                 .filterBounds(roi)
//                 .filterDate('2000-02-26','2019-12-31')
//                 .map(MaskMODISSnow)
//                 .select('NDSI');
var NDSICollection = ee.ImageCollection('MODIS/006/MOD10A1')
  .map(MaskMODISSnow)
  .filterDate('2000-02-26','2019-12-31')
  .map(function(image) {
    var ndsi = image.select('NDSI').multiply(0.0001).rename('NDSI')
      .reproject('EPSG:25833', null, 250);
  return image.addBands(ndsi).set('system:time_start', image.get('system:time_start'));
}); 

// ~~~~~~~~ NDVI ~~~~~~~
var NDVICollection = ee.ImageCollection('MODIS/006/MOD09Q1')
  .map(MaskNDVI)
  .filterDate('2000-02-26','2019-12-31')
  .map(function(image) {
    var ndvi = image.normalizedDifference(['sur_refl_b02', 'sur_refl_b01']).rename('NDVI')
      .reproject('EPSG:25833', null, 250);
    var date = image.date().format('YYYY-MM-dd')
  return image.addBands(ndvi).set('date', date);
}); 

// ~~~~~~~~ LST ~~~~~~~
var LSTCollection = ee.ImageCollection('MODIS/006/MOD11A2')
  // .map(createTimeBand)
  .map(MaskLST)
  .filterDate('2000-02-25','2019-12-31')
  .map(function(image) {
    var lst = image.select('LST_Day_1km')
    .reproject('EPSG:25833', null,250)
      .multiply(0.02).subtract(273.15).rename('LST');
  return ee.Image(1).addBands(lst).set('system:time_start', image.get('system:time_start'));  
  }); 



// Export results with all rois 
var NDVI_ROIS = NDVICollection.select('NDVI').map(function(image) {
  return image.reduceRegions({
    collection: rois,
    reducer: ee.Reducer.mean(),
    scale: 250
  }).map(function(f) {
    // Add a date property to each output feature.
    return f.set('date', image.date().format("YYYY-MM-dd"))
  })
})

var export_ndvi = ee.FeatureCollection(NDVI_ROIS).flatten()

Export.table.toDrive({
  collection: export_ndvi,
  description: 'export-ndvi-modis',
  fileNamePrefix: 'NDVI_00_18_in_rois',
  folder: 'phenology',
  fileFormat: 'CSV'
});

// NDSI
// Export results with all rois 
var NDSI_ROIS = NDSICollection.select('NDSI').map(function(image) {
  return image.reduceRegions({
    collection: rois,
    reducer: ee.Reducer.mean(),
    scale: 250
  }).map(function(f) {
    // Add a date property to each output feature.
    return f.set('date', image.date().format("YYYY-MM-dd"))
  })
})

var export_ndsi = ee.FeatureCollection(NDSI_ROIS).flatten()


// Export.table.toDrive({
//   collection: export_ndsi,
//   description: 'export-ndsi-modis',
//   fileNamePrefix: 'NDSI_00_18_in_rois',
//   folder: 'phenology',
//   fileFormat: 'CSV'
// });


// LST export 
var LST_ROIS = LSTCollection.select('LST').map(function(image) {
  return image.reduceRegions({
    collection: rois,
    reducer: ee.Reducer.mean(),
    scale: 250
  }).map(function(f) {
    // Add a date property to each output feature.
    return f.set('date', image.date().format("YYYY-MM-dd"))
  })
})

var export_lst = ee.FeatureCollection(LST_ROIS).flatten()

Export.table.toDrive({
  collection: export_lst,
  description: 'export-lst-modis',
  fileNamePrefix: 'LST_00_18_in_rois',
  folder: 'phenology',
  fileFormat: 'CSV'
});


// NDSI in watershed
// var NDSI_watershed = ui.Chart.image.seriesByRegion({
//   imageCollection: NDSICollection.select('NDSI'),
//   regions: watershed,
//   reducer: ee.Reducer.mean(),
//   scale: 250
// }).setOptions({
//       title: 'NDVI time series at ROIs',
//     });

// print(NDSI_watershed)


var NDVI_watershed = ui.Chart.image.seriesByRegion({
  imageCollection: NDVICollection.select('NDVI'),
  regions: watershed,
  xProperty: 'date',
  reducer: ee.Reducer.mean(),
  scale: 250
}).setOptions({
      title: 'NDVI time series at ROIs',
    });

print(NDVI_watershed)
