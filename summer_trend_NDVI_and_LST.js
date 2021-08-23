// link https://code.earthengine.google.com/ef39e8b838f3ad70643d3594df36e502 

Map.setCenter(18.874082062258594,68.35693288294851, 9);

Map.addLayer(rois_insitu)

var ROI =  rois_insitu

var start = ee.Date('2000-01-01');
var end = ee.Date('2019-12-31');

var months = ee.List.sequence(1, 12);
var years = ee.List.sequence(2000, 2019);

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


var createTimeBand = function(image) {
  return image.addBands(image.metadata('system:time_start').subtract(start.millis())); //Adds band with date 
};

// -------------------------------------
//            Collections
// -------------------------------------
// NDVI
var summerNDVI = ee.ImageCollection('MODIS/006/MOD09Q1')
  .filter(ee.Filter.calendarRange(2000,2019,'year'))
  .filter(ee.Filter.calendarRange(6,8,'month')) //Filtrer sommermåneder
  .map(createTimeBand)
  .filterBounds(ROI)
  .map(MaskNDVI)
  .map(function(image) {
    var date = ee.Date(image.get('system:time_start'));
    var ndvi = image.normalizedDifference(['sur_refl_b02', 'sur_refl_b01']).rename('NDVI')
      .reproject('EPSG:25833', null, 250);
    var yearOffset = ee.Image(date.difference(ee.Date(start), 'year')).rename('yearOffset');
    var yearOffset2 = ee.Number(date.difference(start, 'year'))
        var system = image.getNumber('system:time_start').toDouble()
    return ee.Image(ndvi).addBands(yearOffset).toDouble().addBands(system).set('system:time_start', image.get('system:time_start')).set('date', date).set('yearOffset',yearOffset2)
});


//LST
var summerLST = ee.ImageCollection('MODIS/006/MOD11A2')
  .filterBounds(ROI)
  .map(MaskLST)
 .filter(ee.Filter.calendarRange(2000,2019,'year'))
 .filter(ee.Filter.calendarRange(6,8,'month'))
  .map(createTimeBand)
  .select('LST_Day_1km') 
  .filterDate('2000-01-01','2019-12-31')
  .map(function(image) {
    var date = ee.Date(image.get('system:time_start'));
    var lst = image.reproject('EPSG:25833', null, 500)
      .multiply(0.02).subtract(273.15).rename('LST');
    var yearOffset = ee.Image(date.difference(ee.Date(start), 'year')).rename('yearOffset');
    var yearOffset2 = ee.Number(date.difference(start, 'year'))

    return ee.Image(lst).addBands(yearOffset).toDouble().set('system:time_start', image.get('system:time_start')).set('date', date).set('yearOffset',yearOffset2)
  });

// -------------------------------------
//            Annual summer values
// -------------------------------------
// NDVI
var annualNDVI = ee.ImageCollection.fromImages(
      years.map(function (y) {
        return summerNDVI
        .filter(ee.Filter.calendarRange(y, y, 'year'))
        .select('NDVI')
        .mean()
        .set('year', y);
}));

// LST
var annualLST = ee.ImageCollection.fromImages(
      years.map(function (y) {
        return summerLST
        .filter(ee.Filter.calendarRange(y, y, 'year'))
        .select('LST')
        .mean()
        .set('year', y);
}));

// Apply date to annual data
var withDateAnnualSummerNDVI = annualNDVI
  .map(function(image) {
    var dateString = ee.String('01-07-').cat(image.get('year'));
    var replaceComma = dateString.replace('\\.0','').replace('\\.0','');
    var date = ee.Date.parse('d-MM-YYYY', replaceComma);
    var yearOffset = ee.Image(date.difference(start, 'year')).rename('yearOffset');
    var yearOffset2 = ee.Number(date.difference(start, 'year'))
    return ee.Image(image.addBands(yearOffset).toDouble().set('date',date).set('system:time_start',date.millis()).set('yearOffset',yearOffset2));
});

var withDateAnnualSummerLST = annualLST
  .map(function(image) {
    var dateString = ee.String('01-07-').cat(image.get('year'));
    // var dateString = image.get('year').cat()
    var replaceComma = dateString.replace('\\.0','').replace('\\.0','');
    var date = ee.Date.parse('d-MM-YYYY', replaceComma);
    var yearOffset = ee.Image(date.difference(start, 'year')).rename('yearOffset');
    return ee.Image(image.addBands(yearOffset).toDouble().set('date',date).set('system:time_start',date.millis()));
});


// -------------------------------------
//            Linear fits
// -------------------------------------
// Non-aggregated
var linearFitSummerNDVI = summerNDVI.select(['yearOffset', 'NDVI'])
  .reduce(ee.Reducer.linearFit());
  
var linearFitSummerNDVI_00_11 = summerNDVI.filterDate('2000-01-01','2011-12-31').select(['yearOffset', 'NDVI'])
  .reduce(ee.Reducer.linearFit().unweighted());


var linearFitSummerNDVI_12_19 = summerNDVI.filterDate('2012-01-01','2019-12-31').select(['yearOffset', 'NDVI'])
  .reduce(ee.Reducer.linearFit());  

// reggresion
var linearFitSummerNDVI_00_11_linereg = summerNDVI.filterDate('2000-01-01','2011-12-31').select(['yearOffset', 'NDVI'])
  .reduce(ee.Reducer.linearRegression({
    numX: 1,
    numY: 1
}));

// Map.addLayer(linearFitSummerNDVI_00_11,null,'allmont')

  
  

var linearFitSummerLST = summerLST.select(['yearOffset', 'LST'])
  .reduce(ee.Reducer.linearFit());

var linearFitSummerLST_00_11 = summerLST.filterDate('2000-01-01','2011-12-31').select(['yearOffset', 'LST'])
  .reduce(ee.Reducer.linearFit());
  
var linearFitSummerLST_12_19 = summerLST.filterDate('2012-01-01','2019-12-31').select(['yearOffset', 'LST'])
  .reduce(ee.Reducer.linearFit());
print(linearFitSummerLST_12_19)

var combo_lst_ndvi = summerLST.filterDate('2000-01-01','2011-12-31').combine(summerNDVI.filterDate('2000-01-01','2011-12-31'), true);
// print(combo_lst_ndvi);

var linearFitLST_NDVI = combo_lst_ndvi.select(['LST', 'NDVI'])
  .reduce(ee.Reducer.linearFit());

Map.addLayer(linearFitLST_NDVI.select('scale'),null, 'Correlation')

// Annual summer values 
var linearFitAnnualSummerNDVI = withDateAnnualSummerNDVI.select(['yearOffset', 'NDVI'])
  .reduce(ee.Reducer.linearFit());
// Map.addLayer(linearFitAnnualSummerNDVI,null,'annual')


var linearFitAnnualSummerLST = withDateAnnualSummerLST.select(['yearOffset', 'LST'])
  .reduce(ee.Reducer.linearFit());
  
  
var combo_lst_ndvi_annual = withDateAnnualSummerLST.filterDate('2000-01-01','2011-12-31').combine(withDateAnnualSummerNDVI.filterDate('2000-01-01','2011-12-31'), true);
// print(combo_lst_ndvi);

// var linearFitLST_NDVI_annual = combo_lst_ndvi_annual.select(['LST', 'NDVI'])
//   .reduce(ee.Reducer.linearFit());
// Map.addLayer(linearFitLST_NDVI_annual.select('scale'),null,'correlation annual 00 11')

// -------------------------------------
//           T-test
// -------------------------------------
var formatrendSummerNDVI = summerNDVI.select('NDVI').formaTrend();
var formatrendSummerNDVI_00_11 = summerNDVI.filterDate('2000-01-01','2011-12-31').select('NDVI').formaTrend();
var formatrendSummerNDVI_12_19 = summerNDVI.filterDate('2012-01-01','2019-12-31').select('NDVI').formaTrend();

var formatrendSummerLST = summerLST.select('LST').formaTrend();
var formatrendSummerLST_00_11 = summerLST.filterDate('2000-01-01','2011-12-31').select('LST').formaTrend();
var formatrendSummerLST_12_19 = summerLST.filterDate('2012-01-01','2019-12-31').select('LST').formaTrend();

// print('Summer NDVI 2000-2019 no.',summerNDVI.size())
// print('Summer NDVI 2000-2011 no.',summerNDVI.filterDate('2000-01-01','2011-12-31').size())
// print('Summer NDVI 2012-2019 no.',summerNDVI.filterDate('2012-01-01','2019-12-31').size())

var formatrendAnnualSummerNDVI = withDateAnnualSummerNDVI.select('NDVI').formaTrend();
var formatrendAnnualSummerLST = withDateAnnualSummerLST.select('LST').formaTrend();

// Masks
// All months
var t_testNDVI = formatrendSummerNDVI.select('long-tstat')
var t_testNDVI_00_11 = formatrendSummerNDVI_00_11.select('long-tstat')
var t_testNDVI_12_19 = formatrendSummerNDVI_12_19.select('long-tstat')



var t_testLST = formatrendSummerLST.select('long-tstat')
var t_testLST_00_11 = formatrendSummerLST_00_11.select('long-tstat')
var t_testLST_12_19 = formatrendSummerLST_12_19.select('long-tstat')

// Criteria
var formaMaskNDVI = t_testNDVI.lt(-1.97).or(t_testNDVI.gt(1.97))
var formaMaskNDVI_00_11 = t_testNDVI_00_11.lt(-1.97).or(t_testNDVI_00_11.gt(1.97))
var formaMaskNDVI_12_19 = t_testNDVI_12_19.lt(-1.99).or(t_testNDVI_12_19.gt(1.99))


var formaMaskLST = t_testLST.lt(-1.97).or(t_testLST.gt(1.97))
var formaMaskLST_00_11 = t_testLST_00_11.lt(-1.97).or(t_testLST_00_11.gt(1.97))
var formaMaskLST_12_19 = t_testNDVI_12_19.lt(-1.99).or(t_testNDVI_12_19.gt(1.99))

// Annual summer values
var t_testNDVI_annual = formatrendAnnualSummerNDVI.select('long-tstat')
var t_testLST_annual = formatrendAnnualSummerLST.select('long-tstat')

var formaMaskNDVI_annual = t_testNDVI_annual.lt(-1.97).or(t_testNDVI_annual.gt(1.97))
var formaMaskLST_annual = t_testLST_annual.lt(-1.97).or(t_testLST_annual.gt(1.97))

// -------------------------------------
//            Pearson
// -------------------------------------
// All months
var PearsonSummerNDVI = summerNDVI.select('yearOffset','NDVI').reduce(ee.Reducer.pearsonsCorrelation());
var PearsonSummerNDVI_00_11 = summerNDVI.filterDate('2000-01-01','2011-12-31').select('yearOffset','NDVI').reduce(ee.Reducer.pearsonsCorrelation());
var PearsonSummerNDVI_12_19 = summerNDVI.filterDate('2012-01-01','2019-12-31').select('yearOffset','NDVI').reduce(ee.Reducer.pearsonsCorrelation());

var PearsonSummerLST = summerLST.select('yearOffset','LST').reduce(ee.Reducer.pearsonsCorrelation());
var PearsonSummerLST_00_11 = summerLST.filterDate('2000-01-01','2011-12-31').select('yearOffset','LST').reduce(ee.Reducer.pearsonsCorrelation());
var PearsonSummerLST_12_19 = summerLST.filterDate('2012-01-01','2019-12-31').select('yearOffset','LST').reduce(ee.Reducer.pearsonsCorrelation());


// Annual values 
var PearsonAnnualSummerNDVI = withDateAnnualSummerNDVI.select('yearOffset','NDVI').reduce(ee.Reducer.pearsonsCorrelation());
var PearsonAnnualSummerLST = withDateAnnualSummerLST.select('yearOffset','LST').reduce(ee.Reducer.pearsonsCorrelation());


// Masks
// All months
var P_NDVI = PearsonSummerNDVI.select('p-value')
var P_NDVI_00_11 = PearsonSummerNDVI_00_11.select('p-value')
var P_NDVI_12_19 = PearsonSummerNDVI_12_19.select('p-value')


var P_LST = PearsonSummerLST.select('p-value')
var P_LST_00_11 = PearsonSummerLST_00_11.select('p-value')
var P_LST_12_19 = PearsonSummerLST_12_19.select('p-value')



var PearsonMaskNDVI = P_NDVI.gt(-0.05).and(P_NDVI.lt(0.05))
var PearsonMaskNDVI_00_11 = P_NDVI_00_11.gt(-0.05).and(P_NDVI_00_11.lt(0.05))
var PearsonMaskNDVI_12_19 = P_NDVI_12_19.gt(-0.05).and(P_NDVI_12_19.lt(0.05))

var PearsonMaskLST = P_LST.gt(-0.05).and(P_LST.lt(0.05))
var PearsonMaskLST_00_11 = P_LST_00_11.gt(-0.05).and(P_LST_00_11.lt(0.05))
var PearsonMaskLST_12_19 = P_LST_12_19.gt(-0.05).and(P_LST_12_19.lt(0.05))

// Annual summer values
var Peason_NDVI_annual = PearsonAnnualSummerNDVI.select('p-value')
var PearsonLST_annual = PearsonAnnualSummerLST.select('p-value')

var PearsonMaskNDVI_annual = Peason_NDVI_annual.gt(-0.05).and(Peason_NDVI_annual.lt(0.05))
var PearsonMaskLST_annual = PearsonLST_annual.gt(-0.05).and(PearsonLST_annual.lt(0.05))

// -------------------------------
//            Apply masks
// -------------------------------
// T-test
// No aggregation
var scaleSumerNDVI = linearFitSummerNDVI.select('scale')
var scaleSumerNDVI_00_11 = linearFitSummerNDVI_00_11.select('scale')
var scaleSumerNDVI_12_19 = linearFitSummerNDVI_12_19.select('scale')

var scaleSummerLST = linearFitSummerLST.select('scale')
var scaleSummerLST_00_11 = linearFitSummerLST_00_11.select('scale')
var scaleSummerLST_12_19 = linearFitSummerLST_12_19.select('scale')


var summerNDVITrend = scaleSumerNDVI.updateMask(formaMaskNDVI);
var summerNDVITrend_00_11 = scaleSumerNDVI_00_11.updateMask(formaMaskNDVI_00_11);
// uden scale
var summerNDVITrend_00_11_med_alt = linearFitSummerNDVI_00_11.updateMask(formaMaskNDVI_00_11);
var summerNDVITrend_12_19 = scaleSumerNDVI_12_19.updateMask(formaMaskNDVI_12_19);


var summerLSTTrend = scaleSummerLST.updateMask(formaMaskLST);
var summerLSTTrend_00_11 = scaleSummerLST_00_11.updateMask(formaMaskLST_00_11);
var summerLSTTrend_12_19 = scaleSummerLST_12_19.updateMask(formaMaskLST_12_19);

// Annual summer value
var annualsummerNDVItrend = linearFitAnnualSummerNDVI.select('scale')
var annualsummerLSTtrend = linearFitAnnualSummerLST.select('scale')
var annualSummerNDVITrend = annualsummerNDVItrend.updateMask(formaMaskNDVI);
var annualSummerLSTTrend = annualsummerLSTtrend.updateMask(formaMaskLST);

// Pearson
// No aggregation
var summerNDVITrend_pearson = scaleSumerNDVI.updateMask(PearsonMaskNDVI);
var summerNDVITrend_pearson_00_11 = scaleSumerNDVI_00_11.updateMask(PearsonMaskNDVI_00_11);
var summerNDVITrend_pearson_12_19 = scaleSumerNDVI_12_19.updateMask(PearsonMaskNDVI_12_19);


var summerLSTTrend_pearson = scaleSummerLST.updateMask(PearsonMaskLST);
var summerLSTTrend_pearson_00_11 = scaleSummerLST_00_11.updateMask(PearsonMaskLST_00_11);
var summerLSTTrend_pearson_12_19 = scaleSummerLST_12_19.updateMask(PearsonMaskLST_12_19);

// Annual summer value
var annualSummerNDVITrend_pearson = annualsummerNDVItrend.updateMask(PearsonMaskNDVI_annual);
var annualSummerLSTTrend_pearson = annualsummerLSTtrend.updateMask(PearsonMaskLST_annual);


// View fits 
var palette = {
  palette: ['2166ac','67a9cf','d1e5f0','ffffbf','fddbc7','ef8a62','b2182b'],
  min: -1,
  max: 1,
}

var no1 = {palette : ['b2182b']}
var no2 = {palette : ['ffffbf']}
var no3 = {palette : ['2166ac']}

// Map.addLayer(summerNDVITrend.select('scale'),no1,'NDVI 00-19',false);
// Map.addLayer(summerNDVITrend_00_11.select('scale'),no2,'NDVI 00-11',false);
// Map.addLayer(summerNDVITrend_pearson_00_11.select('scale'),no3,'NDVI 00-11 - pearson',false);
// Map.addLayer(summerLSTTrend.select('scale'),no1,'LST 00-19',false);
// Map.addLayer(summerLSTTrend_00_11.select('scale'),no2,'LST 00-11',false);
Map.addLayer(summerLSTTrend_pearson_00_11.select('scale'),no3,'LST 00-11 - pearson',false);
Map.addLayer(summerLSTTrend_pearson_12_19.select('scale'),no3,'LST 12-19 - pearson',false);
Map.addLayer(summerLSTTrend_12_19.select('scale'),no2,'LST 12-19 - t-test');
// Map.addLayer(annualSummerLSTTrend.select('scale'),no3,'Annual LST 00-18');

print(summerNDVITrend_00_11_med_alt)
// Export images
Export.image.toDrive({
  image: summerNDVITrend_00_11_med_alt,
  description: 'summerNDVITrend_00_11_med_alt',
  fileNamePrefix: 'summerNDVITrend_00_11_med_alt',
  folder: 'summer values',
  region: geometry, 
  scale: 250
});

// Export.image.toDrive({
//   image: summerLSTTrend_12_19,
//   description: 'summerLSTTrend_12_19',
//   fileNamePrefix: 'summerLSTTrend_12_19',
//   folder: 'summer values',
//   region: geometry, 
//   scale: 1000
// });
// Export.image.toDrive({
//   image: summerLSTTrend_pearson,
//   description: 'summerLSTTrend_pearson_00_19',
//   fileNamePrefix: 'summerLSTTrend_pearson_00_19',
//   folder: 'summer values',
//   region: geometry, 
//   scale: 1000
// });


// Export.image.toDrive({
//   image: summerLSTTrend_pearson_00_11,
//   description: 'summerLSTTrend_pearson_00_11',
//   fileNamePrefix: 'summerLSTTrend_pearson_00_11',
//   folder: 'summer values',
//   region: geometry, 
//   scale: 1000
// });






// // Charts of annual summer values
// var annualNDVIChart = ui.Chart.image.seriesByRegion({
//   imageCollection: summerNDVI.select('NDVI'),
//   xProperty: 'yearOffset',
//   regions: ROI,
//   reducer: ee.Reducer.mean(),
//   scale: 250
// }).setOptions({
//       title: 'NDVI time series at ROIs',
//     });

// var annualLSTChart = ui.Chart.image.seriesByRegion({
//   imageCollection: annualLST.select('LST'),
//   xProperty: 'year',
//   regions: ROI,
//   reducer: ee.Reducer.mean(),
//   scale: 250
// }).setOptions({
//       title: 'LST time series at ROIs',
//     });
    


// // Combine tables for export
// var summerValues = summerNDVI.combine(summerLST,true)
// var summerAnnualValues = annualNDVI.combine(annualLST,true)

// // print(summerValues)
// // Export features to csv
// var export_summerValues = summerValues.select('NDVI','LST').map(function(image) {
//   return image.reduceRegions({
//     collection: ROI,
//     reducer: ee.Reducer.mean(),
//     scale: 250
//   }).map(function(f) {
//     // Add a date property to each output feature.
//     return f.set('date', image.date().format("YYYY-MM-dd"))
//   })
// })

// var summerCSV = ee.FeatureCollection(export_summerValues).flatten()
// // print(summerCSV)
// var export_annualSummerValues = summerAnnualValues.select('NDVI','LST').map(function(image) {
//   return image.reduceRegions({
//     collection: ROI,
//     reducer: ee.Reducer.mean(),
//     scale: 250
//   }).map(function(f) {
//     // Add a date property to each output feature.
//     return f.set('year', image.get('year'))
//   })
// })

// var annualSummerCSV = ee.FeatureCollection(export_annualSummerValues).flatten()


// // Export tables 
// Export.table.toDrive({
//   collection: summerCSV,
//   description: 'LST_and_NDVI_in_summer_all_months-ndvi-modis',
//   fileNamePrefix: 'LST_and_NDVI_in_summer_all_months',
//   folder: 'summer values',
//   fileFormat: 'CSV'
// });

// // Export.table.toDrive({
// //   collection: annualSummerCSV,
// //   description: 'annual_summer_values_00_19_maj_sep',
// //   fileNamePrefix: 'annual_summer_values_00_19_maj_sep',
// //   folder: 'summer values',
// //   fileFormat: 'CSV'
// // });

