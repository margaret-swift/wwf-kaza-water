/**~ * ~ * ~ * ~ * ~ * ~ * ~ * ~ * ~ * ~ * ~ * ~ * 
 * WaterDefinitionWorkflow.js
 * Created: 05 August 2025
 * Author:  Margaret Swift 
 * Contact: <mes473@cornell.edu> 
 * Website: www.maggie.earth
 * GITHub:  https://github.com/margaret-swift/wwf-kaza-water/blob/main/WaterDefinitionWorkflow.js
 * ~ * ~ * ~ * ~ * ~ * ~ * ~ * ~ * ~ * ~ * ~ * ~ * 
 * 
 * This code uses Otsu thresholding on median AWEI values, calculated from Sentinel-2
 *    MSI imagery, to generate a binary "isWater" raster where: 
 *      1 = water
 *      0 = not water
 * 
 * I've put in several "SPEED BOOKMARK" notes if you're having issues with performance or memory
 * 
 * The data are exported to a Google Drive folder. This code was originally written to interact
 *    with the Google Cloud Storage platform.  If you have questions about batch execution for 
 *    multiple AOI (I ran mine over the WWF Hydrosheds database), email me at <mes473@cornell.edu>
 * 
 * MAIN CODE STARTS ~ LINE 120
 */
 
 
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
//                                    !!! LOADING DATA FILES !!!
// -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

// LOAD PEKEL GLOBAL SURFACE WATER DATASET
// https://www.nature.com/articles/nature20584
var JRC = ee.Image("JRC/GSW1_2/GlobalSurfaceWater").select('recurrence').selfMask();

// LOAD SENTINEL DATA
var s2 = ee.ImageCollection('COPERNICUS/S2_HARMONIZED');
var s2c = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY');
var s2sr = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');

// LOAD MODIS BURNED AREA PRODUCT
var burns = ee.ImageCollection('MODIS/006/MCD64A1');

// LOAD BUILDINGS
var buildings = ee.ImageCollection('GOOGLE/Research/open-buildings-temporal/v1')
  .filterDate('2023-01-01', '2023-12-31')
  .select('building_presence')
  .reduce(ee.Reducer.median())
  .rename('isBuilding')

// LOAD WORLD SETTLEMENT FOOTPRINT
var wsf = ee.Image('DLR/WSF/WSF2015/v1').unmask()

// LOAD ROADS FROM OSM 
// you can download from here: https://www.openstreetmap.org/#map=5/38.01/-95.84
// IF YOU USE ROAD MASKING: You need to make sure you're zoomed in quite far so that the painting mask 
//   doesn't get too large (it's pixel-based)
// var roads = ee.FeatureCollection('projects/kaza-waterhole-mapping/assets/osm-files/primary_trunk_rds_5countries');



// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
//                               !!! USER DEFINED PARAMETERS !!!
// -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

//..........................................................................
// USER-DEFINED AOI AND WINDOW-CENTER POINT
var aoi = ee.Geometry.Polygon(
            [[[23.362973058178042, -18.457400995456705],
              [23.362973058178042, -18.646180731915678],
              [23.621838414623355, -18.646180731915678],
              [23.621838414623355, -18.457400995456705]]], null, false);
var poi = ee.Geometry.Point([23.497416585821693, -18.537746016891045]);


//..........................................................................
// EXPORT FILES?
//..........................................................................

// yes/no
var doExport = false;

// BATCH EXPORT OVER WWF HYDROSHEDS, YES OR NO?
var runBatches = false;

// SINGLE AOI EXPORT PARAMETERS
var xname = "waterhole_fill_AOI"

// BATCH EXPORT PARAMETERS
// var listSize = 25
// var hydro = ee.FeatureCollection("WWF/HydroSHEDS/v1/Basins/hybas_8");
// var batch_aoi = hydro.filterBounds(aoi)
// var id_name = "HYBAS_ID"
// var bname = "kaza_waterhole_rasters" // GCS bucket name
// var xname_slug = 'maxfill_h8' // export URI prefix


//..........................................................................
// USER-DEFINED WET SEASON
//..........................................................................

// We want to define maximum water fill at the wettest time of the wettest years.
// this will differ by location, so we let the user define the wet season. 
// I count backward from the end of the wet season, since wet seasons in KAZA
// span the previous November and December. You can change this code to fit
// the climate of your region of interest.

var wet_end = '-05-01' // end of wet season
var wet_length = -5 // number of months to go back from end of wet season
function defineWetSeason (y) {
  var d2 = ee.Date(ee.String(ee.Number(y)).cat(wet_end));
  var d1 = d2.advance(wet_length, 'month');
  return(ee.Dictionary({'d1':d1, 'd2':d2}));
}



// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
//                               !!! MAIN CODE EXPLORATION !!!
// -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

//..........................................................................
// PLOTTING
//..........................................................................

// SET UP MAP OPTIONS
// Map needs to be zoomed in for proper road masking and pixel-based calcs
Map.centerObject(poi, 17) 
Map.setOptions('SATELLITE');
var aweiVis = {"bands":["AWEIsh_median"],"min":-0.74,"max":0.014};

// PLOTTING AWEI
// a wet year - 2022
var dates_22 = defineWetSeason(2022)
var s2_22 = s2Prep(dates_22, aoi).select("AWEIsh_median");
Map.addLayer(s2_22, aweiVis, "AWEIsh wet year - 2022")
// a dry year - 2019
var dates_19 = defineWetSeason(2019)
var s2_19 = s2Prep(dates_19, aoi).select("AWEIsh_median");
Map.addLayer(s2_19, aweiVis, "AWEIsh dry year - 2019", false)

//..........................................................................
// CREATE MAXIMUM WATER FILL DATASET
//..........................................................................
var water_max = createMaxPotentialFill(aoi);
Map.addLayer(water_max, {'palette':'yellow'}, "Maximum Water Fill")


//..........................................................................
// AN EXAMPLE IN DRY YEAR FOR COMPARISON TO MAX WATER FILL
//..........................................................................
var debug = true; // If you run Otsu in debug mode, you get a histogram in the console
var ot = computeThresholdUsingOtsu(s2_19, aoi, debug)
var water_19 = s2_19.where(s2_19.lt(ot), 0).where(s2_19.gte(ot), 1).selfMask()
Map.addLayer(water_19, {'palette':'red'}, "Dry year (2019) water fill")

// we can zoom back out now
Map.centerObject(poi, 14) 

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
//                                        !!! EXPORT !!!
// -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

//..........................................................................
// SINGLE AOI EXPORT
//..........................................................................
if (!runBatches & doExport) {
  Export.image.toDrive({
    image: water_max,
    description: xname,
    region: aoi,
    scale: 10,
    maxPixels: 10000000000000
  }); 

//..........................................................................
// BATCH EXPORT INTO BUCKETS
//..........................................................................
/* 
 * GOOGLE CLI CODE TO MOVE FILES: 
 *    gcloud storage mv gs://kaza_waterhole_rasters/SLUG* gs://kaza_waterhole_rasters/TO_FOLDER/
 * 
 * BATCH TASK EXECUTION - STEPS
 *    https://benny.istan.to/blog/20220319-batch-task-execution-in-google-earth-engine-code-editor

 * PASTE INTO DEVELOPER CONSOLE
      function runTaskList(){
        // var tasklist = document.getElementsByClassName('task local type-EXPORT_IMAGE awaiting-user-config');
        // for (var i = 0; i < tasklist.length; i++)
        //         tasklist[i].getElementsByClassName('run-button')[0].click();
        $$('.run-button' ,$$('ee-task-pane')[0].shadowRoot).forEach(function(e) {
             e.click();
        })
      }
      runTaskList();
*/
} else if (runBatches & doExport) {
  
  // get list of IDS
  var all_ids = batch_aoi.aggregate_array(id_name).sort()
  var ids_list = ee.List.sequence(0, all_ids.size(), listsize);
  
  //  BATCH TASK EXECUTION
  // batchExport(1) // 
  // batchExport(2) // 
  // ... for as many as you have.
}

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
//                                    !!! FUNCTIONS !!!
// -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

//..........................................................................
//   PLOTTING
//..........................................................................
function outlineObject(obj, name, col, w) {
  w = w || 1; col = col || 'yellow'
  /* outlines a geometry */
  var outline = ee.Image().byte().paint({
    featureCollection:ee.FeatureCollection([obj]),color:1,width:w
  });
  Map.addLayer(outline, {palette:col}, name);
}

//..........................................................................
//   BAND ALGEBRA
//..........................................................................
function panSharpen(params) {
  /* 
    Pan-sharpen SWIR bands: https://leclab.wixsite.com/spatial/post/pansharpening-sentinel-2-imagery-in-google-earth-engine
  */
  var image = params.image;
  var geometry = params.geometry || image.geometry();
  var crs = params.crs || image.select(0).projection();
  var maxPixels = params.maxPixels;
  var bestEffort = params.bestEffort || true;
  var tileScale = params.tileScale || 1;
  image = image.clip(geometry);
  var bands10m = ['B2', 'B3', 'B4', 'B8'];
  var bands20m = ['B5', 'B6', 'B7', 'B8A', 'B11', 'B12'];
  var panchromatic = image
    .select(bands10m)
    .reduce(ee.Reducer.mean());
  var image20m = image.select(bands20m);
  var image20mResampled = image20m.resample('bilinear');
  var stats20m = image20m
    .reduceRegion({
      reducer: ee.Reducer.stdDev().combine(
        ee.Reducer.mean(), null, true
      ),
      geometry: geometry,
      scale: 20,
      crs: crs, 
      bestEffort: bestEffort, 
      maxPixels: maxPixels, 
      tileScale: tileScale
    })
    .toImage();
    
  var mean20m = stats20m
    .select('.*_mean')
    .regexpRename('(.*)_mean', '$1');
    
  var stdDev20m = stats20m
    .select('.*_stdDev')
    .regexpRename('(.*)_stdDev', '$1');
    
  var kernel = ee.Kernel.fixed({
    width: 5,
    height: 5, 
    weights: [
      [-1, -1, -1, -1, -1],
      [-1, -1, -1, -1, -1],
      [-1, -1, 24, -1, -1],
      [-1, -1, -1, -1, -1],
      [-1, -1, -1, -1, -1]
    ], 
    x: -3, 
    y: -3, 
    normalize: false
  });
  
  var highPassFilter = panchromatic
    .convolve(kernel)
    .rename('highPassFilter');

  var stdDevHighPassFilter = highPassFilter
    .reduceRegion({
      reducer: ee.Reducer.stdDev(),
      geometry: geometry,
      scale: 10,
      crs: crs, 
      bestEffort: bestEffort, 
      maxPixels: maxPixels, 
      tileScale: tileScale
    })
    .getNumber('highPassFilter');
  
  function calculateOutput(bandName) {
    bandName = ee.String(bandName);
    var W = ee.Image().expression(
      'stdDev20m / stdDevHighPassFilter * modulatingFactor', {
        stdDev20m: stdDev20m.select(bandName),
        stdDevHighPassFilter: stdDevHighPassFilter,
        modulatingFactor: 0.25
      }
    );
    return ee.Image()
      .expression(
        'image20mResampled + (HPF * W)', {
          image20mResampled: image20mResampled.select(bandName),
          HPF: highPassFilter,
          W: W
      }
    )
    .uint16();
  }
  
  var output = ee.ImageCollection(
      bands20m.map(calculateOutput)
    )
    .toBands()
    .regexpRename('.*_(.*)', '$1');
    
  var statsOutput = output
    .reduceRegion({
      reducer: ee.Reducer.stdDev().combine(
        ee.Reducer.mean(), null, true
      ),
      geometry: geometry,
      scale: 10,
      crs: crs, 
      bestEffort: bestEffort, 
      maxPixels: maxPixels, 
      tileScale: tileScale
    })
    .toImage();
    
  var meanOutput = statsOutput
    .select('.*_mean')
    .regexpRename('(.*)_mean', '$1');
    
  var stdDevOutput = statsOutput
    .select('.*_stdDev')
    .regexpRename('(.*)_stdDev', '$1');
  
  var sharpened = ee.Image()
    .expression(
      '(output - meanOutput) / stdDevOutput * stdDev20m + mean20m', {
        output: output,
        meanOutput: meanOutput,
        stdDevOutput: stdDevOutput,
        stdDev20m: stdDev20m,
        mean20m: mean20m
      }
    )
    .uint16() ;
  
  return image
    .addBands(sharpened, null, true)
    .select(image.bandNames());
}  
function addIndex(img, exp, dict, bname) {
  /*  Add an index band to an existing image. */
  var inx = img.expression(exp, dict).rename(bname).select(bname);
  return img.addBands(inx);
}
function addAWEIsh(img){
  /* 
    Define AWEI (remember to pan-sharpen and divide by 10000)
    https://www.sciencedirect.com/science/article/abs/pii/S0034425713002873
  */
  // SPEED BOOKMARK: Ignore pan-sharpening for Otsu thresholding to improve memory capacity
  // var sharp = img.divide(10000) // 40sec
  var sharp = panSharpen({image: img}).divide(10000); // 2m
  var AWEIsh = sharp.expression(
    'B2 + (2.5 * B3) - (1.5 * (B8 + B11)) - (0.25 * B12)', {
      'B2': sharp.select('B2'),
      'B3': sharp.select('B3'),
      'B8': sharp.select('B8'),
      'B11': sharp.select('B11'),
      'B12': sharp.select('B12')
    }).rename('AWEIsh').select('AWEIsh');
  return img.addBands(AWEIsh);
}

//..........................................................................
//   MASKING AND FILTERING FUNCTIONS
//..........................................................................

function s2Prep (dates, bounds) {
  /* Prepare your sentinel-2 data by masking out clouds, roads, etc., adding AWEI, 
      and taking median */
  var d1 = dates.get('d1'); 
  var d2 = dates.get('d2');
  var s2_masked = s2Cloudless(d1, d2, bounds) 
    .map(addAWEIsh)       // add critical AWEI band
    .map(maskBurns)       // dynamic
    // .map(maskRoads)       // static
    .map(maskSettlements) // static
    .map(maskBuildings)   // static
  var s2_median = ee.Image(s2_masked.reduce(ee.Reducer.median()));
  return(s2_median);
}
function s2Cloudless (d1, d2, bounds) {
  /* 
    Filters clouds out of Sentinel-2 imagery; use the parameters to tune to your liking.
    https://code.earthengine.google.co.in/76d9e10ea376b73fb16776a6d9d14588
    More about s2Cloudless: https://developers.google.com/earth-engine/tutorials/community/sentinel-2-s2cloudless
    
    SPEED BOOKMARK: I created maskImageSimple to improve performance; maskImage takes a long time.
  */
  var CLD_FILTER_THRESH = 10;
  var CLOUD_PROB_THRESH = 20;
  function maskImageSimple(image) {
    var s2c = image.select('probability');
    var cirrus = image.select('B10').multiply(0.0001);
    var isCloud = s2c.gt(CLOUD_PROB_THRESH).or(cirrus.gt(0.01));
    return image.updateMask(isCloud.not());
  }
  // Filter dates/bounds and join Sentinel-2 data
  function indexJoin(colA, colB, pname) {
    /* Join two collections on their 'system:index' property. */
    var joined = ee.ImageCollection(ee.Join.saveFirst(pname).apply({
      primary: colA, secondary: colB,
      condition: ee.Filter.equals({ leftField: 'system:index', rightField: 'system:index'})}));
    return joined.map(function(image) { return image.addBands(ee.Image(image.get(pname))) });
  }
  var s2_filtered = s2.filterDate(d1, d2).filterBounds(bounds);
  var s2c_filtered = s2c.filterDate(d1, d2).filterBounds(bounds);
  var s2sr_filtered = s2sr.filterDate(d1, d2).filterBounds(bounds);
  var withCloudProbability = indexJoin(s2sr_filtered, s2c_filtered, 'cloud_probability');
  var withS2L1C = indexJoin(withCloudProbability, s2_filtered, 'l1c');

  // Filter out cloud threshold
  var s2_cloudthresh = withS2L1C.filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CLD_FILTER_THRESH))
  var s2_cloudsmasked = ee.ImageCollection(s2_cloudthresh.map(maskImageSimple));
  return(s2_cloudsmasked);
}
function maskBurns(img) {
  /* 
    Use MODIS burn scarring data to remove burn scars within the last two months.
    I don't know why, but filterDate() doesn't quite work with this dataset so 
    I had to create a bit of my own to help out.
  */
  // this is the most updated version fixing the MODIS issue
  var endDate = img.date().advance(5, 'day');
  var startDate = ee.Date(endDate).advance(-2, 'month');
  // get the year and start of earliest year
  var endY = endDate.get('year');
  var startY = startDate.get('year');
  var Jan01 = ee.Date(ee.String(startY).cat('-01-01'))
  // function to wrap day of year
  var wrapDOY = function(b_img) {
    var year = b_img.date().get('year')
    var bd = b_img.select('BurnDate')
    if (year.gt(Jan01.get('year')) == 1) { bd = bd.add(365) }
    return(bd)
  }
  var burns_filtered = burns
    .filterDate(startDate.advance(-1, 'month'), endDate)
    .filterBounds(img.geometry())
    .select('BurnDate')
    .map(wrapDOY)
    .mosaic()
    
  // Find the start and end DOY for which we should mask
  var startDOY = startDate.difference(Jan01, 'day')
  var endDOY = endDate.difference(Jan01, 'day')
  var burnMask = burns_filtered
    .where(burns_filtered.lte(endDOY).and(burns_filtered.gte(startDOY)), 1)
    .where(burns_filtered.gt(endDOY).or(burns_filtered.lt(startDOY)), 0)
    .unmask()
  var burnMaskFlip = burnMask.where(burnMask.eq(0), 1).where(burnMask.eq(1), 0)
  return(img.updateMask(burnMaskFlip));
}
function maskRoads (image) { 
  /*  Mask out roads */
  var roadMask  = ee.Image(1).int().paint(roads,0,5);
  var masked = image.updateMask(roadMask)
  return(masked);
}
function maskSettlements (image) {
  /* Mask out human settlements */
  var wsfMask = wsf.where(wsf.eq(0), 1).where(wsf.neq(0), 0);
  return(image.updateMask(wsfMask));
}
function maskBuildings (image) {
  /* 
    Mask out anything with a 50% chance or higher of being a building.
  */
  var bt = 0.5;
  var buildMask = buildings.where(buildings.lte(bt), 0).where(buildings.gt(bt), 1);
  var buildMaskBuffer = ee.Image(1)
    .cumulativeCost({
      source: buildMask, 
      maxDistance: 10,
    }).lt(10); // keep this for mapping if you want.
  var buildMaskReverseBuffer = ee.Image(1)
    .where(buildMaskBuffer.eq(1), 0).where(buildMaskBuffer.eq(0), 1)
  return(image.updateMask(buildMaskReverseBuffer));
}

//..........................................................................
//   SORTING YEARS BY REGIONAL MEAN ANNUAL RAINFALL
//..........................................................................
function sortYearsByRainfall (bounds) {
  /* Uses CHIRPS data to give mmMAR for the region */
  var geom = ee.Geometry(bounds.geometry());
  var computeTotalMM = function(y) {
    var d2 = ee.Date.fromYMD(y,5,1);
    var d1 = d2.advance(-6, 'month');
    var chirps = ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY")
      .filterDate(d1, d2)
      .sum()
      .clip(geom);
    var data = chirps
      .reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: geom,
        maxPixels: 70000000,
        scale: 10,  // meters
        bestEffort : true
    }).get('precipitation');
    return(data);
  };
  var mm_dict = ee.Dictionary({
    2019:computeTotalMM(2019),
    2020:computeTotalMM(2020),
    2021:computeTotalMM(2021),
    2022:computeTotalMM(2022),
    2023:computeTotalMM(2023),
    2024:computeTotalMM(2024)
  });
  // sort years by mm of rainfall
  var sorted_keys = mm_dict.keys().sort(mm_dict.values()).reverse(); //descending order
  return(sorted_keys);
}

//..........................................................................
//   OTSU ALGORITHM AND THRESHOLD CALCULATION
//..........................................................................

// DONCHYTS 2016 paper: http://www.mdpi.com/2072-4292/8/5/386
// CODE: https://code.earthengine.google.com/e9c699d7b5ef811d0a67c02756473e9d
function otsu (histogram) {
  histogram = ee.Dictionary(histogram);
  var counts = ee.Array(histogram.get('histogram'));
  var means = ee.Array(histogram.get('bucketMeans'));
  var size = means.length().get([0]);
  var total = counts.reduce(ee.Reducer.sum(), [0]).get([0]);
  var sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]);
  var mean = sum.divide(total);
  var indices = ee.List.sequence(1, size);

  // Compute between sum of squares, where each mean partitions the data.
  var bss = indices.map(function(i) {
      var aCounts = counts.slice(0, 0, i);
      var aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0]);
      var aMeans = means.slice(0, 0, i);
      var aMean = aMeans.multiply(aCounts)
          .reduce(ee.Reducer.sum(), [0]).get([0])
          .divide(aCount);
      var bCount = total.subtract(aCount);
      var bMean = sum.subtract(aCount.multiply(aMean)).divide(bCount);
      return aCount.multiply(aMean.subtract(mean).pow(2)).add(
          bCount.multiply(bMean.subtract(mean).pow(2)));
  });

  // Return the mean value corresponding to the maximum BSS.
  return means.sort(bss).get([-1]);
};
function computeThresholdUsingOtsu( image, bounds, debug ) { 
  /*** Compute a threshold using Otsu method (bimodal) */
  var DEBUG = debug || false;
  var CANNY_THRESH = 0.7;
  var CANNY_SIGMA = 1;
  var SCALE = 100;
  var BUCKETS = 100;
  var MIN_VALUE = -1;
  var MIN_EDGE_LEN = 4;0
  
  // clip image edges
  var mask = image.mask().gt(0).clip(bounds).focal_min(ee.Number(SCALE).multiply(3), 'circle', 'meters');
  
  // detect sharp changes
  var edge = ee.Algorithms.CannyEdgeDetector(image, CANNY_THRESH, CANNY_SIGMA);
  edge = edge.multiply(mask);
  if(MIN_EDGE_LEN) {
      var connected = edge.mask(edge).lt(CANNY_THRESH).connectedPixelCount(200, true);
      var edgeLong = connected.gte(MIN_EDGE_LEN);
      edge = edgeLong;
  }
  // buffer around NDWI edges
  var edgeBuffer = edge.focal_max(ee.Number(SCALE), 'square', 'meters');
  edge = edge.updateMask(edgeBuffer)
  edgeBuffer = edge.focal_max(ee.Number(SCALE).multiply(1), 'square', 'meters');
  var imageEdge = image.mask(edgeBuffer);
  
  // compute threshold using Otsu thresholding
  var hist = ee.Dictionary(ee.Dictionary(imageEdge.reduceRegion(ee.Reducer.histogram(BUCKETS), bounds, SCALE)).values().get(0));
  if (hist.contains('bucketMeans')) { var threshold = ee.Number(otsu(hist));
  } else { print("WARNING: BUCKETMEANS ARE EMPTY!"); threshold=null; }
  if(DEBUG) {
    print("-----------------------------------------")
    print("EXAMPLE OTSU THRESHOLD")
    print("-----------------------------------------")
      // Map.addLayer(edge.mask(edge), {palette:['ff0000']}, 'edges', true);
      print('Otsu Threshold: ', threshold);
      print('Image values:', ui.Chart.image.histogram(image, bounds, SCALE, BUCKETS));
      var watermap = image
        .where(image.lte(threshold), 0) // non-water
        .where(image.gt(threshold), 1)  // water
        .rename('isWater')
        .selfMask();
      Map.addLayer(watermap, {'palette':'yellow'}, "Water detected by AWEIsh", false)
      print("-----------------------------------------")
  }
  return threshold
}

//..........................................................................
//   WATERHOLE CREATION WORKHORSE
//..........................................................................
function createWaterholesbyYear (bounds, y) {
  /* Create a waterhole dataset for the wet season of the given year. */
  
  // Get median s2 image over these dates and the bounds
  var dates_dict = defineWetSeason(y)
  var s2_median = s2Prep(dates_dict, bounds).select("AWEIsh_median");
  print(ee.String(y).cat(" wet season dates: "), dates_dict)
  
  // Calculate Otsu threshold. This can take a while.
  var ot_water = computeThresholdUsingOtsu(s2_median, bounds);
  print(ee.String(y).cat(" Otsu threshold: "))
  print(ot_water)
  
  // Mask out any AWEI value that is below the calculated Otsu threshold, 
  //  then self-mask and return.
  var watermap = s2_median
    .where(s2_median.lte(ot_water), 0) // non-water
    .where(s2_median.gt(ot_water), 1)  // water
    .rename('isWater')
    .selfMask();
  return(watermap)
}
function createMaxPotentialFill (bounds) {
  /* 
    This function combines thresholded AWEI rasters AND the Pekel 2016 GSW product
      into one "maximum water fill" dataset that can be used to mask other imagery.
  */
  var feat = ee.Feature(bounds)

  // Sort years by most rainfall -> least rainfall for the region
  var years_lst = sortYearsByRainfall(feat);
  
  // For two wettest years, create waterhole map and save to collection
  var y1 = ee.Number.parse(years_lst.get(0));
  var wh1 = createWaterholesbyYear(bounds, y1);
  var y2 = ee.Number.parse(years_lst.get(1));
  var wh2 = createWaterholesbyYear(bounds, y2);
  var wh_coll = ee.ImageCollection.fromImages([ wh1, wh2 ]);
  
  // create waterhole collection and mosaic into one
  var waterholes = ee.ImageCollection( wh_coll.sum() );
  var water_msc = waterholes.mosaic();
  
  // Add GSW recurrence layer to waterholes layer and sum, then set any "isWater"
  //  value that is greater than 0 to be 1.
  var rec = JRC.select('recurrence').rename('isWater')
  var water_sum = ee.ImageCollection.fromImages([ water_msc, rec ]).sum();
  var all_water = water_sum.where(water_sum.gt(0), 1).select('isWater');
  return(all_water);
}

//..........................................................................
//   BATCH EXPORT FUNCTIONS
//..........................................................................
function exportFunc(id, slug) {
  var aoi_id = aoi.filter(ee.Filter.eq(id_name, id));
  var waterholes = createMaxPotentialFill(aoi_id);
  var xname = slug + id;
  Export.image.toCloudStorage({
    image: waterholes,
    description: xname,
    bucket: bname,
    region: aoi_id,
    scale: 10,
    maxPixels: 10000000000000
  }); 
}
var batchExport = function(inx) {
  var inx_lo = ids_list.get(inx-1)
  var inx_hi = ids_list.get(inx)
  var bth_lo = all_ids.get(inx_lo)
  var bth_hi = all_ids.get(inx_hi) 
  var batch_ids = aoi
    .filter(ee.Filter.gt(id_name, bth_lo))
    .filter(ee.Filter.lte(id_name, bth_hi))
    .aggregate_array(id_name)
  batch_ids.evaluate(function (ids) {
    ids.map(function (id) exportFunc(id, xname_slug) )
  });
}



// EOF */
