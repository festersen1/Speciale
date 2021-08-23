Map.setCenter(19.125,68.220703125, 10);
var geometry = ee.Geometry.Rectangle([20.8212890625,68.748046875, 17.3408203125,67.552734375]);

// ROIS = https://code.earthengine.google.com/?asset=users/festersen/rois_for_phenology
// DOC watershed = https://code.earthengine.google.com/?asset=users/festersen/watershed

// ==========================================================================================
//                                  Quality assesment 
// ==========================================================================================
var start = ee.Date('2000-01-01');
var end = ee.Date('2018-12-31');

// Mask: 
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
// ~~~~~~~~ NDVI QA ~~~~~~~~ 
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

// ~~~~~~~~ LST QA ~~~~~~~~
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

// ---------------------------------- Apply mask and ready collections ----------------------------------
var createTimeBand = function(image) {
  return image.addBands(image.metadata('system:time_start').subtract(start.millis())); //Adds band with date 
};

// ~~~~~~~~ NDVI ~~~~~~~
var NDVICollection = ee.ImageCollection('MODIS/006/MOD09Q1')
  .map(createTimeBand)
  .filterDate('2000-01-01','2018-12-31')
  .map(MaskNDVI)
  .map(function(image) {
    var ndvi = image.normalizedDifference(['sur_refl_b02', 'sur_refl_b01']).rename('NDVI')
      .reproject('EPSG:25833', null, 250);
  return ee.Image(1).addBands(ndvi).set('system:time_start', image.get('system:time_start'));  
}); 


//~~~~~~~~ LST ~~~~~~~
var LSTCollection = ee.ImageCollection('MODIS/006/MOD11A2')
  .map(createTimeBand)
  .map(MaskLST)
  .map(function(image) {
    var lst = image.select('LST_Day_1km')
    .reproject('EPSG:25833', null,250)
      .multiply(0.02).subtract(273.15).rename('LST');
  return ee.Image(1).addBands(lst).set('system:time_start', image.get('system:time_start'));  
  }); 

// ==========================================================================================
//                            De-seasoning 
// ==========================================================================================
// Calculate mean for months, for all 19 years (2000-2018):
var months = ee.List.sequence(1, 12);
var years = ee.List.sequence(2000, 2018);

//~~~~~~~~ Aggregation to monthly data ~~~~~~~~
//NDVI
var byMonthYearNDVI = ee.ImageCollection.fromImages(
  years.map(function(y) {
    return months.map(function (m) {
      return NDVICollection
        .select('NDVI')
        .filter(ee.Filter.calendarRange(y, y, 'year'))
        .filter(ee.Filter.calendarRange(m, m, 'month'))
        .mean()
        .set('month', m).set('year', y)
  });
}).flatten());

// LST
var byMonthYearLST = ee.ImageCollection.fromImages(
  years.map(function(y) {
    return months.map(function (m) {
      return LSTCollection
        .select('LST')
        .filter(ee.Filter.calendarRange(y, y, 'year'))
        .filter(ee.Filter.calendarRange(m, m, 'month'))
        .mean()
        .set('month', m).set('year', y);
  });
}).flatten());

//~~~~~~~~ Average annual values ~~~~~~~~
// Calculate mean for year
var byYearNDVI = ee.ImageCollection.fromImages(
      years.map(function (y) {
        return NDVICollection
        .filter(ee.Filter.calendarRange(y, y, 'year'))
        .select('NDVI')
        .mean()
        .set('year', y);
}));

var byYearLST = ee.ImageCollection.fromImages(
      years.map(function (y) {
        return LSTCollection
        .filter(ee.Filter.calendarRange(y, y, 'year'))
        .select('LST')
        .mean()
        .set('year', y);
}));

// Calculate mean for every month, every year:
var byMonthNDVI = ee.ImageCollection.fromImages(
      months.map(function (m) {
        return NDVICollection
        .filter(ee.Filter.calendarRange(m, m, 'month'))
        .select('NDVI')
        .mean()
        .set('month', m);
}));

var byMonthLST = ee.ImageCollection.fromImages(
     months.map(function (m) {
        return LSTCollection
        .filter(ee.Filter.calendarRange(m, m, 'month'))
        .select('LST')
        .mean()
        .set('month', m);
}));

// Filter images by geometry to ensure collections are inside study area:
var monthlyNDVI = byMonthYearNDVI.filterBounds(geometry);
var annualAverageNDVI = byYearNDVI.filterBounds(geometry); 
var monthlyAverageNDVI = byMonthNDVI.filterBounds(geometry);

var monthlyLST = byMonthYearLST.filterBounds(geometry); 
var annualAverageLST = byYearLST.filterBounds(geometry);
var monthlyAverageLST = byMonthLST.filterBounds(geometry);

//  ------------- De-season by anomalies  -------------
var anomaliesNDVI = monthlyNDVI.filterMetadata('system:index',"not_equals", 0).select('NDVI')
  .map(function(image) {
  var x = image.get('month');
  var y = monthlyAverageNDVI.select('NDVI').filterMetadata('month',"equals", x).first();
  return ee.Image(y.select('NDVI').subtract(image.select('NDVI'))).rename('deseasoned_ndvi_an').set('month',
  image.get('month')).set('year', image.get('year'));
});

var anomaliesLST = monthlyLST.filterMetadata('system:index',"not_equals", 0).select('LST')
  .map(function(image) {
  var x = image.get('month');
  var y = monthlyAverageLST.select('LST').filterMetadata('month',"equals", x).first();
  return ee.Image(y.select('LST').subtract(image.select('LST'))).rename('deseasoned_lst_an').set('month', 
  image.get('month')).set('year', image.get('year'));
});

//  -------------Proportions (monthly values divided by annual mean) -------------
//~~~~~~~~ NDVI ~~~~~~~
var proportionsNDVI = monthlyNDVI
  .map(function(image) {
  // Date - used for aggregating to seasonal index
    var dateString = ee.String('1').cat('-').cat(image.get('month')).cat('-').cat(image.get('year'));
    var replaceComma = dateString.replace('\\.0','').replace('\\.0','');
    var date = ee.Date.parse('d-MM-YYYY', replaceComma);
    var yearOffset = ee.Image(date.difference(start, 'year')).rename('yearOffset');
  
  // Calculation of proportions
  var x = image.get('year');
  var y = annualAverageNDVI.filterMetadata('year',"equals", x).first();
  var proportions = image.divide(y);
  return ee.Image(proportions.toDouble().rename('proportions')).addBands(yearOffset).toDouble().set('date',date)
  .set('system:time_start',date.millis());
});

//~~~~~~~~ LST ~~~~~~~
var proportionsLST = monthlyLST
  .map(function(image) {
  // Date - used for aggregating to seasonal index
    var dateString = ee.String('1').cat('-').cat(image.get('month')).cat('-').cat(image.get('year'));
    var replaceComma = dateString.replace('\\.0','').replace('\\.0','');
    var date = ee.Date.parse('d-MM-YYYY', replaceComma);
    var yearOffset = ee.Image(date.difference(start, 'year')).rename('yearOffset');
  
  // Calculation of proportions
  var x = image.get('year');
  var y = annualAverageLST.filterMetadata('year',"equals", x).first();
  var proportions = image.divide(y);
  return ee.Image(proportions.toDouble().rename('proportions')).addBands(yearOffset).toDouble().set('date',date)
  .set('system:time_start',date.millis());
});


// ------------- Calculate seasonal index -------------
// Proportions are averaged for every month:
//~~~~~~~~ NDVI ~~~~~~~
var SiNDVI = ee.ImageCollection.fromImages(
      months.map(function (m) {
        return proportionsNDVI
        .filter(ee.Filter.calendarRange(m, m, 'month'))
        .select('proportions').mean()
        .set('month', m)
        .rename('monthly_proportions');
}));

//~~~~~~~~ LST ~~~~~~~
var SiLST = ee.ImageCollection.fromImages(
      months.map(function (m) {
        return proportionsLST
        .filter(ee.Filter.calendarRange(m, m, 'month'))
        .select('proportions').mean()
        .set('month', m)
        .rename('monthly_proportions');
}));

// ------------- Correct time serie values with the seasonal index -------------
// All values are multiplied by the seasonal index 

//~~~~~~~~ NDVI ~~~~~~~
var deseasonedNDVI = monthlyNDVI
  .map(function(image) {
  var x = image.get('month');
  var y = SiNDVI.filterMetadata('month',"equals", x).first();
  var diff = image.divide(y).rename('deseasoned_ndvi');
  return ee.Image(diff).set('month', image.get('month')).set('year', image.get('year'));
});

//~~~~~~~~ LST ~~~~~~~
var deseasonedLST = monthlyLST
  .map(function(image) {
  var x = image.get('month');
  var y = SiLST.filterMetadata('month',"equals", x).first();
  var diff = image.divide(y).rename('deseasoned_lst');
  return ee.Image(diff).set('month', image.get('month')).set('year', image.get('year'));
});

// ---------------------------------- Create date proporties ----------------------------------
// Anomali method
//~~~~~~~~ NDVI ~~~~~~~
var withDateNDVI_an = anomaliesNDVI
  .map(function(image) {
    var dateString = ee.String('1').cat('-').cat(image.get('month')).cat('-').cat(image.get('year'));
    var replaceComma = dateString.replace('\\.0','').replace('\\.0','');
    var date = ee.Date.parse('d-MM-YYYY', replaceComma);
    var yearOffset = ee.Image(date.difference(start, 'year')).rename('yearOffset');
    return ee.Image(image.addBands(yearOffset).toDouble().set('date',date).set('system:time_start',date.millis()));
});

//~~~~~~~~ LST ~~~~~~~
var withDateLST_an = anomaliesLST
  .map(function(image) {
    var dateString = ee.String('1').cat('-').cat(image.get('month')).cat('-').cat(image.get('year'));
    var replaceComma = dateString.replace('\\.0','').replace('\\.0','');
    var date = ee.Date.parse('d-MM-YYYY', replaceComma);
    var yearOffset = ee.Image(date.difference(start, 'year')).rename('yearOffset');
    return ee.Image(image.addBands(yearOffset).toDouble().set('date',date).set('system:time_start',date.millis()));
});

// Index method
//~~~~~~~~ NDVI ~~~~~~~
var withDateNDVI = deseasonedNDVI
  .map(function(image) {
    var dateString = ee.String('1').cat('-').cat(image.get('month')).cat('-').cat(image.get('year'));
    var replaceComma = dateString.replace('\\.0','').replace('\\.0','');
    var date = ee.Date.parse('d-MM-YYYY', replaceComma);
    var yearOffset = ee.Image(date.difference(start, 'year')).rename('yearOffset');
    return ee.Image(image.addBands(yearOffset).toDouble().set('date',date).set('system:time_start',date.millis()));
});

//~~~~~~~~ LST ~~~~~~~
var withDateLST = deseasonedLST
  .map(function(image) {
    var dateString = ee.String('1').cat('-').cat(image.get('month')).cat('-').cat(image.get('year'));
    var replaceComma = dateString.replace('\\.0','').replace('\\.0','');
    var date = ee.Date.parse('d-MM-YYYY', replaceComma);
    var yearOffset = ee.Image(date.difference(start, 'year')).rename('yearOffset');
    return ee.Image(image.addBands(yearOffset).toDouble().set('date',date).set('system:time_start',date.millis()));
});

// ==========================================================================================
//                    Analysis of data 
// ==========================================================================================
// --------------------------------------------------------------------
//                       Plot charts
// --------------------------------------------------------------------
var NDVI_both_deseasoned = withDateNDVI.combine(withDateNDVI_an,true)

var deseaonedNDVIinROIS = ui.Chart.image.seriesByRegion({
  imageCollection: NDVI_both_deseasoned.select('deseasoned_ndvi','deseasoned_ndvi_an'),
  xProperty: 'date',
  regions: ROIS,
  reducer: ee.Reducer.mean(),
  scale: 250
}).setOptions({
      title: 'De-seasoned NDVI time series at ROIs',
    });
  
  
var deseaonedLSTinROIS = ui.Chart.image.seriesByRegion({
  imageCollection: withDateLST.select('deseasoned_lst'),
  xProperty: 'date',
  regions: ROIS,
  reducer: ee.Reducer.mean(),
  scale: 1000
}).setOptions({
      title: 'De-seasoned LST time series at ROIs',
    });
    
// print(deseaonedNDVIinROIS)


// --------------------------------------------------------------------
//                       Trend analysis 
// --------------------------------------------------------------------
//  Linear trends
// ~~~~~~~~  Index  ~~~~~~~~ 
// Calculate linear trend (linear regression with date as independent value)
var linTrendNDVI = withDateNDVI.select(['yearOffset', 'deseasoned_ndvi'])
  .reduce(ee.Reducer.linearFit());

// print('Linear trend NDVI',linTrendNDVI)
var linTrendLST = withDateLST.select(['yearOffset', 'deseasoned_lst'])
  .reduce(ee.Reducer.linearFit());

//Merge image collections to make correlations 
var NDVI_and_LST = withDateNDVI
.combine(withDateLST)
.map(function(image) { 
  return image.clip(geometry); //Clip to decrease processing time 
});

var linTrendBoth = NDVI_and_LST.select(['LST', 'NDVI'])
  .reduce(ee.Reducer.linearFit());
  
//  ~~~~~~~~ Anomalies ~~~~~~~~ 
// Calculate linear trend (linear regression with date as independent value)
var linTrendNDVI_an = withDateNDVI_an.select(['yearOffset', 'deseasoned_ndvi_an'])
  .reduce(ee.Reducer.linearFit());

// print('Linear trend NDVI',linTrendNDVI)
var linTrendLST_an = withDateLST_an.select(['yearOffset', 'deseasoned_lst_an'])
  .reduce(ee.Reducer.linearFit());

//Merge image collections to make correlations 
var NDVI_and_LST_an = withDateNDVI_an
.combine(withDateLST_an)
.map(function(image) { 
  return image.clip(geometry); //Clip to decrease processing time 
});

var linTrendBoth_an = NDVI_and_LST_an.select(['LST', 'NDVI'])
  .reduce(ee.Reducer.linearFit());

// Select scale (slope of linear trend)
// Index
var ScaleNDVI = linTrendNDVI.select('scale');
var ScaleLST = linTrendLST.select('scale');
var ScaleBoth = linTrendBoth.select('scale');

// Anomalies 
var ScaleNDVI_an = linTrendNDVI_an.select('scale');
var ScaleLST_an = linTrendLST_an.select('scale');
var ScaleBoth_an = linTrendBoth_an.select('scale');


// ---------------------------------- Calculate significance and apply masks ----------------------------------
//  ~~~~~~~~ t-test  ~~~~~~~~
// FormaTrend test - caluclates t-test on residuals 
// Index
var formatrendNDVI = withDateNDVI.select('deseasoned_ndvi').formaTrend();
var formatrendLST = withDateLST.select('deseasoned_lst').formaTrend();


// Anomalies
var formatrendNDVI_an = withDateNDVI_an.select('deseasoned_ndvi_an').formaTrend();
var formatrendLST_an = withDateLST_an.select('deseasoned_lst_an').formaTrend();


//  ~~~~~~~~ Pearson test  ~~~~~~~~
// Index
var PearsonNDVI = withDateNDVI.select('yearOffset','deseasoned_ndvi').reduce(ee.Reducer.pearsonsCorrelation());
var PearsonLST = withDateLST.select('yearOffset','deseasoned_lst').reduce(ee.Reducer.pearsonsCorrelation());

// Anomali
var PearsonNDVI_an = withDateNDVI_an.select('yearOffset','deseasoned_ndvi_an').reduce(ee.Reducer.pearsonsCorrelation());
var PearsonLST_an = withDateLST_an.select('yearOffset','deseasoned_lst_an').reduce(ee.Reducer.pearsonsCorrelation());


// ~~~~~~~~ Masks  ~~~~~~~~
// Variation mask
// Aggregate to annual data and find mean of all 19 years
var VarMaskNDVI = ee.ImageCollection.fromImages(
      years.map(function (y) {
        return NDVICollection
        .filter(ee.Filter.calendarRange(y, y, 'year'))
        .select(1).reduce(ee.Reducer.stdDev())
        .set('year', y);
}));

var stdMaskNDVI = VarMaskNDVI.mean();

var VarMaskLST = ee.ImageCollection.fromImages(
      years.map(function (y) {
        return LSTCollection
        .filter(ee.Filter.calendarRange(y, y, 'year'))
        .select(1).reduce(ee.Reducer.stdDev())
        .set('year', y);
}));

var stdMaskLST = VarMaskLST.mean();

// ---------- Significance t-test mask ---------- 
// Index
var t_testNDVI = formatrendNDVI.select('long-tstat')
var t_testLST = formatrendLST.select('long-tstat')

var formaMaskNDVI = t_testNDVI.gt(1.97).or(t_testNDVI.lt(-1.97))
.and(stdMaskNDVI.select('NDVI_stdDev').gt(0.1));

var formaMaskLST = t_testLST.gt(1.97).or(t_testLST.lt(-1.97))
.and(stdMaskNDVI.select('NDVI_stdDev').gt(0.1));

// Anomalies
var t_testNDVI_an = formatrendNDVI_an.select('long-tstat')
var t_testLST_an = formatrendLST_an.select('long-tstat')

var formaMaskNDVI_an = t_testNDVI_an.gt(1.97).or(t_testNDVI_an.lt(-1.97))
.and(stdMaskNDVI.select('NDVI_stdDev').gt(0.1));

var formaMaskLST_an = t_testLST_an.gt(1.97).or(t_testLST_an.lt(-1.97))
.and(stdMaskNDVI.select('NDVI_stdDev').gt(0.1));

// ---------- Significance Pearson mask ---------- 
// Index
var P_NDVI = PearsonNDVI.select('p-value')
var P_LST = PearsonLST.select('p-value')

var PearsonMaskNDVI = P_NDVI.gt(-0.05).and(P_NDVI.lt(0.05))
var PearsonMaskLST = P_LST.gt(-0.05).and(P_LST.lt(0.05))


// Anomali
var P_NDVI_an = PearsonNDVI_an.select('p-value')
var P_LST_an = PearsonNDVI_an.select('p-value')


var PearsonMaskNDVI_an = P_NDVI_an.gt(-0.05).and(P_NDVI_an.lt(0.05))
var PearsonMaskLST_an = P_LST_an.gt(-0.05).and(P_LST_an.lt(0.05))


// ~~~~~~~~  Apply masks ~~~~~~~~ 
// T-test
// Index
var NDVIScaleWithMask = ScaleNDVI.updateMask(formaMaskNDVI);
var LSTScaleWithMask = ScaleLST.updateMask(formaMaskLST);

// Anomalies
var NDVIScaleWithMask_an = ScaleNDVI_an.updateMask(formaMaskNDVI_an);
var LSTScaleWithMask_an = ScaleLST_an.updateMask(formaMaskLST_an);

// Pearson correlation 
// Index
var NDVIScaleWithMask_pearson = ScaleNDVI.updateMask(PearsonMaskNDVI);
var LSTScaleWithMask_pearson  = ScaleLST.updateMask(PearsonMaskLST);

// Anomalies
var NDVIScaleWithMask_an_pearson  = ScaleNDVI_an.updateMask(PearsonMaskNDVI_an);
var LSTScaleWithMask_an_pearson  = ScaleLST_an.updateMask(PearsonMaskLST_an);



// ---------------------------------- Show linear trends  ----------------------------------
var VisParamNDVI = {
  palette: ['e33521', 'ffff00', '27a929'],
  min: -0.6,
  max: 0.7,
};

var VisParamLST = {
  palette: ['27a929', 'ffff00', 'e33521'],
  min: -100,
  max: 50,
};

Map.addLayer(NDVIScaleWithMask,VisParamNDVI,'NDVI lin trend - index method t-test',false);
Map.addLayer(NDVIScaleWithMask_pearson,VisParamNDVI,'NDVI lin trend - index method Pearson',false);

// Map.addLayer(NDVIScaleWithMask_an,VisParamNDVI,'NDVI anomali trend - t-test',false);
// Map.addLayer(NDVIScaleWithMask_an_pearson,VisParamNDVI,'NDVI anomali trend - pearson ',false);


// Map.addLayer(LSTScaleWithMask,VisParamLST,'LST lin trend - index method t-test',false);
// Map.addLayer(LSTScaleWithMask_pearson,VisParamLST,'LST lin trend - index method Pearson',false);


// Map.addLayer(LSTScaleWithMask_an,VisParamLST,'LST anomali trend - t-test',false);
// Map.addLayer(LSTScaleWithMask_an_pearson,VisParamLST,'LST anomali trend - pearson ',false);


// charts

var AnnualMeanNDVI__2 = ee.ImageCollection.fromImages(
      years.map(function (y) {
        return NDVICollection
        .filter(ee.Filter.calendarRange(y, y, 'year'))
        .select('NDVI')
        .mean()
        .set('year', y);
}));


print(ui.Chart.image.series(AnnualMeanNDVI__2.select('NDVI'), geometry4, ee.Reducer.mean(), 250,'year'));
// print(ui.Chart.image.series(byYearNDVI.select('NDVI'), geometry2, ee.Reducer.mean(), 250),'year');


// ---------- Export linear trends with masks ----------
// Combine scaled images for export 
// Without masks
var ScaleCombineNDVI = ScaleNDVI.select('scale').rename('scale_ndvi')
.addBands(ScaleNDVI_an.select('scale').rename('scale_ndvi_an'));

var ScaleCombineLST = ScaleLST.select('scale')
.rename('scale_lst').addBands(ScaleLST_an.select('scale').rename('scale_lst_an'));

// With masks 
var NDVIScaleWithMask_combined = NDVIScaleWithMask.select('scale').rename('scale_ndvi')
.addBands(NDVIScaleWithMask_an.select('scale').rename('scale_ndvi_an'));

var LSTScaleWithMask_combined = LSTScaleWithMask.select('scale')
.rename('scale_lst').addBands(LSTScaleWithMask_an.select('scale').rename('scale_lst_an'));

var LSTScaleWithMask_combined_pearson = LSTScaleWithMask_pearson.select('scale')
.rename('scale_lst').addBands(LSTScaleWithMask_an_pearson.select('scale').rename('scale_lst_an'));


// Export images of masked correlations
Export.image.toDrive({  
          image: NDVIScaleWithMask,
          folder: 'Final data GEE', 
          description:'NDVIScaleWithMask',
          fileNamePrefix:'NDVI_modis_00_12_deseasoned_index',  
          scale: 250, 
          region: geometry3 
    }); 
  
  
// Export.image.toDrive({  
//           image: LSTScaleWithMask_combined_pearson,
//           folder: 'Final data GEE', 
//           description:'LSTScaleWithMask_combined_pearson',
//           fileNamePrefix:'LST_modis_00_19_deseasoned_trends_pearson',  
//           scale: 1000, 
//           region: geometry3 
//     }); 
  

// Print band names to ensure right bands are used 
// print(NDVIScaleWithMask_combined.bandNames())

// =========================================================================================================
// Export of data in ROIS. Series of de-seasoned NDVI and LST as well as non-de-seasoned data is exported 
// for further analysis 
// =========================================================================================================

// ---------- Add date proporties to monthly values ----------
var withDateNDVI_no_deseason = monthlyNDVI
  .map(function(image) {
    var dateString = ee.String('1').cat('-').cat(image.get('month')).cat('-').cat(image.get('year'));
    var replaceComma = dateString.replace('\\.0','').replace('\\.0','');
    var date = ee.Date.parse('d-MM-YYYY', replaceComma);
    var yearOffset = ee.Image(date.difference(start, 'year')).rename('yearOffset');
    return ee.Image(image.addBands(yearOffset).toDouble().set('date',date).set('system:time_start',date.millis()));
});

var withDateLST_no_deseason = monthlyLST
  .map(function(image) {
    var dateString = ee.String('1').cat('-').cat(image.get('month')).cat('-').cat(image.get('year'));
    var replaceComma = dateString.replace('\\.0','').replace('\\.0','');
    var date = ee.Date.parse('d-MM-YYYY', replaceComma);
    var yearOffset = ee.Image(date.difference(start, 'year')).rename('yearOffset');
    return ee.Image(image.addBands(yearOffset).toDouble().set('date',date).set('system:time_start',date.millis()));
});

// ---------- Combine different image series ----------
// -- NDVI --
var NDVI_combined1 = withDateNDVI
.combine(withDateNDVI_no_deseason,true)
.map(function(image) { 
  return image.clip(geometry); //Clip to decrease processing time 
});

var NDVI_combined = NDVI_combined1
.combine(withDateNDVI_an,true)
.map(function(image) { 
  return image.clip(geometry); //Clip to decrease processing time 
});

// -- LST --
var LST_combined1 = withDateLST
.combine(withDateLST_no_deseason,true)
.map(function(image) { 
  return image.clip(geometry); //Clip to decrease processing time 
});

var LST_combined = LST_combined1
.combine(withDateLST_an,true)
.map(function(image) { 
  return image.clip(geometry); //Clip to decrease processing time 
});


// ---------- Make dictionary with values to combine in csv ----------
// ~~~~~~~~ NDVI  ~~~~~~~~
var NDVI_in_rois = NDVI_combined.map(function(image) {
  return image.reduceRegions({
    collection: ROIS,
    reducer: ee.Reducer.mean(),
    scale: 250
  }).map(function(f) {
    // Add a date property to each output feature.
    return f.set('date', image.date().format("YYYY-MM-dd"))
  })
})

var exportNDVI = ee.FeatureCollection(NDVI_in_rois).flatten()

var orig_0 = exportNDVI.filterMetadata('ORIG_FID',"equals",0).reduceColumns(ee.Reducer.toList(2), ['yearOffset', 'deseasoned_ndvi']).values().get(0)
var orig_1 = exportNDVI.filterMetadata('ORIG_FID',"equals",1).reduceColumns(ee.Reducer.toList(2), ['yearOffset', 'deseasoned_ndvi']).values().get(0)
var orig_2 = exportNDVI.filterMetadata('ORIG_FID',"equals",2).reduceColumns(ee.Reducer.toList(2), ['yearOffset', 'deseasoned_ndvi']).values().get(0)
var orig_3 = exportNDVI.filterMetadata('ORIG_FID',"equals",3).reduceColumns(ee.Reducer.toList(2), ['yearOffset', 'deseasoned_ndvi']).values().get(0)
var orig_4 = exportNDVI.filterMetadata('ORIG_FID',"equals",4).reduceColumns(ee.Reducer.toList(2), ['yearOffset', 'deseasoned_ndvi']).values().get(0)
var orig_5 = exportNDVI.filterMetadata('ORIG_FID',"equals",5).reduceColumns(ee.Reducer.toList(2), ['yearOffset', 'deseasoned_ndvi']).values().get(0)
var orig_6 = exportNDVI.filterMetadata('ORIG_FID',"equals",6).reduceColumns(ee.Reducer.toList(2), ['yearOffset', 'deseasoned_ndvi']).values().get(0)
var orig_7 = exportNDVI.filterMetadata('ORIG_FID',"equals",7).reduceColumns(ee.Reducer.toList(2), ['yearOffset', 'deseasoned_ndvi']).values().get(0)
var orig_8 = exportNDVI.filterMetadata('ORIG_FID',"equals",8).reduceColumns(ee.Reducer.toList(2), ['yearOffset', 'deseasoned_ndvi']).values().get(0)
var orig_9 = exportNDVI.filterMetadata('ORIG_FID',"equals",9).reduceColumns(ee.Reducer.toList(2), ['yearOffset', 'deseasoned_ndvi']).values().get(0)


var t0 = ee.List(orig_0).reduce(ee.Reducer.pearsonsCorrelation());
var t1 = ee.List(orig_1).reduce(ee.Reducer.pearsonsCorrelation());
var t2 = ee.List(orig_2).reduce(ee.Reducer.pearsonsCorrelation());
var t3 = ee.List(orig_3).reduce(ee.Reducer.pearsonsCorrelation());
var t4 = ee.List(orig_4).reduce(ee.Reducer.pearsonsCorrelation());
var t5 = ee.List(orig_5).reduce(ee.Reducer.pearsonsCorrelation());
var t6 = ee.List(orig_6).reduce(ee.Reducer.pearsonsCorrelation());
var t7 = ee.List(orig_7).reduce(ee.Reducer.pearsonsCorrelation());
var t8 = ee.List(orig_8).reduce(ee.Reducer.pearsonsCorrelation());
var t9 = ee.List(orig_9).reduce(ee.Reducer.pearsonsCorrelation());


var dataTable = {
  correlations: [
        {OB1: t0},
        {WH1:t1},
        {WH2:t2},
        {H1:t3},
        {OB2:t4},
        {H2:t5},
        {WH3:t6},
        {OB3:t7},
        {H3:t8},
        {National_park:t9}
        ]
};
// print(dataTable)

// ~~~~~~~~ LST  ~~~~~~~~
var LST_in_rois = LST_combined.map(function(image) {
  return image.reduceRegions({
    collection: ROIS,
    reducer: ee.Reducer.mean(),
    scale: 250
  }).map(function(f) {
    // Add a date property to each output feature.
    return f.set('date', image.date().format("YYYY-MM-dd"))
  })
})

var exportLST = ee.FeatureCollection(LST_in_rois).flatten()

// Export.table.toDrive({
//   collection: exportLST,
//   description: 'export-LST-deseaoned',
//   fileNamePrefix: 'LST_deseasoned_and_not-in-ROIS',
//   folder: 'Final data GEE',
//   fileFormat: 'CSV'
// });

// Export.table.toDrive({
//   collection: exportNDVI,
//   description: 'export-NDVI-deseaoned',
//   fileNamePrefix: 'NDVI_deseasoned_and_not-in-ROIS',
//   folder: 'Final data GEE',
//   fileFormat: 'CSV'
// });

// ---------------------------------------------
// Combination of NDVI and LST without deseason
// ---------------------------------------------
var combo = LSTCollection.select('LST')
.combine(NDVICollection.select('NDVI'),true)
.map(function(image) { 
  return image.clip(geometry); //Clip to decrease processing time 
});

// ~~~~~~~~ Create feature collection ~~~~~~~~
var NDVI_LST_no_deseason = combo.map(function(image) {
  return image.reduceRegions({
    collection: ROIS,
    reducer: ee.Reducer.mean(),
    scale: 250
  }).map(function(f) {
    // Add a date property to each output feature.
    return f.set('date', image.date().format("YYYY-MM-dd"))
  })
})

var exportBothNoDeseason = ee.FeatureCollection(NDVI_LST_no_deseason).flatten()
// Export.table.toDrive({
//   collection: exportBothNoDeseason,
//   description: 'export-NDVI-deseaoned',
//   fileNamePrefix: 'NDVI_LST_non_deseasoned-in-ROIS',
//   folder: 'Final data GEE',
//   fileFormat: 'CSV'
// });


// ---------------------------------------------
//  Correlation between NDVI and LST
// ---------------------------------------------
// Collections
var NDVICollection_corr = ee.ImageCollection('MODIS/006/MOD09Q1')
  .filterBounds(geometry2)
  .filterDate('2002-01-01', '2018-12-31')
  .map(createTimeBand)
  .map(function(image) {
    var ndvi = image.normalizedDifference(['sur_refl_b02', 'sur_refl_b01']).rename('NDVI')
      .reproject('EPSG:25833', null, 250);
  return ee.Image(1).addBands(ndvi).set('system:time_start', image.get('system:time_start'));  
}); 

var LSTCollection_corr = ee.ImageCollection('MODIS/006/MOD11A2')
  .filterDate('2002-01-01', '2018-12-31')
  .filterBounds(geometry2)
  .map(createTimeBand)
  .map(function(image) {
    var lst = image.select('LST_Day_1km')
    .reproject('EPSG:25833', null,250)
      .multiply(0.02).subtract(273.15).rename('LST');
  return ee.Image(1).addBands(lst).set('system:time_start', image.get('system:time_start'));  
  }); 

// X and Y values (list)
var yvalues = NDVICollection_corr.map(function(i) {
  return i.set(i.reduceRegion({ 
    reducer: ee.Reducer.mean(), 
    geometry: geometry, 
    scale: 1000
  }))
}).aggregate_array('NDVI')


var xvalues = LSTCollection_corr.map(function(i) {
  return i.set(i.reduceRegion({ 
    reducer: ee.Reducer.mean(), 
    geometry: geometry, 
    scale: 1000
  }))
}).aggregate_array('LST');

var chartNDVIandLST =  ui.Chart.array.values(yvalues, 0, xvalues)
    // .setSeriesNames(['LST', 'NDVI'])
    .setOptions({
      title: 'LST vs NDVI 2002-2018',
      hAxis: {'title': 'LST'},
      vAxis: {'title': 'NDVI'},
      pointSize: 3,
      trendlines: { 0: { 
        color: 'green',
        lineWidth: 3,
        opacity: 0.2,
      } }
});

// print(chartNDVIandLST);

// print(ee.List(xvalues).zip(yvalues))
var corr = ee.List(xvalues).zip(yvalues).reduce(ee.Reducer.pearsonsCorrelation())
// print('Correlation coefficient',corr)

// --------------------------------------------------
//               Annual mean of NDVI and LST 
// --------------------------------------------------
// ~~~~~~~~ Annual means  ~~~~~~~~
var AnnualMeanNDVI = ee.ImageCollection.fromImages(
      years.map(function (y) {
        return NDVICollection
        .filter(ee.Filter.calendarRange(y, y, 'year'))
        .select('NDVI')
        .mean()
        .set('year', y);
}));

var AnnualMeanLST = ee.ImageCollection.fromImages(
      years.map(function (y) {
        return LSTCollection
        .filter(ee.Filter.calendarRange(y, y, 'year'))
        .select('LST')
        .mean()
        .set('year', y);
}));

var AnnualMeans = AnnualMeanNDVI.combine(AnnualMeanLST)

// Make charts
var chartAnnualMeans = ui.Chart.image.series({
  imageCollection: AnnualMeans,
  xProperty: 'year',
  region: DOCwatershed,
  reducer: ee.Reducer.mean(),
  scale: 250
}).setOptions({title: 'NDVI and LST annual mean'});

var chartMeanLST = ui.Chart.image.series({
  imageCollection: AnnualMeanLST,
  xProperty: 'year',
  region: DOCwatershed,
  reducer: ee.Reducer.mean(),
  scale: 1000
}).setOptions({title: 'LST annual mean'});

// print(chartAnnualMeans)

// Charts with values in all ROIS 
var AnnualMeanNDVIInROIS = ui.Chart.image.seriesByRegion({
  imageCollection: AnnualMeanNDVI.select('NDVI'),
  xProperty: 'year',
  regions: ROIS,
  reducer: ee.Reducer.mean(),
  scale: 250
}).setOptions({
      title: 'NDVI time series at ROIs',
    });

// print(AnnualMeanNDVIInROIS)

    
// // ------------------------- Thumbnails -------------------------
// // Visualization parameters.
// var args = {
//   crs: 'EPSG:25833',  
//   dimensions: '300',
//   region: geometry,
//   min: -0.1,
//   max: 1, 
//   palette: ['e33521', 'ffff00', '27a929','104a11'],
//   framesPerSecond: 5,
// };

// // NDVI in several years
// var thumbCOL = NDVICollection.filterDate('2000-03-01','2006-10-01')

// // Create a video thumbnail and add it to the map.
// var thumb = ui.Thumbnail({
//   // Specifying a collection for "image" animates the sequence of images.
//   image: thumbCOL.select('NDVI'),
//   params: args,
//   style: {
//     position: 'bottom-right',
//     width: '320px'
//   }});
// // Map.add(thumb);

// // NDVI in august 
// var thumbCOL2 = NDVICollection
// .select('NDVI')
//   .filter(ee.Filter.calendarRange(2000,2018,'year'))
//   .filter(ee.Filter.calendarRange(8,8,'month'));

// var thumb2 = ui.Thumbnail({
//   // Specifying a collection for "image" animates the sequence of images.
//   image: thumbCOL2,
//   params: args,
//   style: {
//     position: 'bottom-right',
//     width: '320px'
//   }});
// // Map.add(thumb2);
