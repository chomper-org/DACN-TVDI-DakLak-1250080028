// =============================================================================
// ĐỒ ÁN CHUYÊN NGÀNH — CHẶNG 1 (PHIÊN BẢN NÂNG CẤP V2)
// Sinh viên: Trương Thành Danh (1250080028) — Lớp 12_THMT
// Đề tài: Đánh giá hạn hán từ ảnh vệ tinh sử dụng NDVI và LST (TVDI)
// Khu vực: Tỉnh Đắk Lắk mới (Hợp nhất Đắk Lắk + Phú Yên)
//
// CÁC NÂNG CẤP ĐÁP ỨNG PHIẾU NHẬN XÉT PHẢN BIỆN:
// 1. Xuất toàn bộ raster ở hệ tọa độ chiếu phẳng UTM Zone 49N (EPSG:32649), scale 60m.
// 2. In thống kê số lượng cảnh ảnh từng năm (tách rõ Landsat 8 vs. Landsat 9 năm 2026).
// 3. Bổ sung lọc chất lượng nhiệt ST_QA (sai số <= 4.0 K) & giữ nguyên NoData (không nội suy).
// 4. Trích xuất bảng CSV khí tượng mùa khô (CHIRPS mưa + ERA5-Land nhiệt độ 2m & độ ẩm đất)
//    tại 10 điểm/trạm đại diện để kiểm chứng định lượng với TVDI.
// 5. Xuất lớp độ cao SRTM DEM, lớp phủ ESA WorldCover và CSV 10 địa danh cho QGIS.
// =============================================================================

// --- 1. RANH GIỚI TỈNH ĐẮK LẮK MỚI (ĐẮK LẮK + PHÚ YÊN) ---
var admin1 = ee.FeatureCollection('FAO/GAUL/2015/level1');
var daklak_fc = admin1.filter(ee.Filter.inList('ADM1_NAME', ['Dak Lak', 'Phu Yen']));
var roi = daklak_fc.union().geometry();

Map.centerObject(roi, 8);
Map.addLayer(roi, {color: 'red'}, 'Ranh gioi Dak Lak moi');

var TARGET_CRS = 'EPSG:32649'; // UTM Zone 49N (WGS 84) — đơn vị mét
var TARGET_SCALE = 60;         // 60m x 60m = 3.600 m2/pixel chuẩn xác
var DRIVE_FOLDER = 'DACN_DakLak_TVDI_V2';

// --- 2. HÀM TIỀN XỬ LÝ LANDSAT 8/9 C2 L2 (QA_PIXEL + ST_QA + HỆ SỐ HIỆU CHUẨN) ---
function preprocessLandsat(image) {
  var qa = image.select('QA_PIXEL');
  var cloudMask = qa.bitwiseAnd(1 << 3).eq(0);       // Bit 3: Cloud
  var cloudShadowMask = qa.bitwiseAnd(1 << 4).eq(0); // Bit 4: Cloud Shadow

  // Kiểm soát chất lượng nhiệt độ bề mặt ST_QA (đơn vị 0.01 K -> lọc sai số <= 4.0 K)
  var stQa = image.select('ST_QA').multiply(0.01);
  var stQaMask = stQa.lte(4.0);

  // Hiệu chuẩn bức xạ theo chuẩn USGS Collection 2 Level 2
  var red = image.select('SR_B4').multiply(0.0000275).add(-0.2);
  var nir = image.select('SR_B5').multiply(0.0000275).add(-0.2);
  var ndvi = nir.subtract(red).divide(nir.add(red)).rename('NDVI');

  var lst = image.select('ST_B10').multiply(0.00341802).add(149.0).subtract(273.15).rename('LST');

  // Lọc ngưỡng vật lý hợp lệ (giữ nguyên NoData ngoài ngưỡng, không nội suy lấp mây)
  var physMask = ndvi.gte(-1.0).and(ndvi.lte(1.0))
                     .and(lst.gt(0.0)).and(lst.lt(60.0));

  var validMask = cloudMask.and(cloudShadowMask).and(stQaMask).and(physMask);

  return image.addBands([ndvi, lst]).select(['NDVI', 'LST']).updateMask(validMask);
}

// --- 3. DANH SÁCH 10 ĐIỂM / TRẠM ĐẠI DIỆN ĐẮK LẮK MỚI (ĐỐI CHIẾU KHÍ TƯỢNG & ĐỊA DANH) ---
var stations = ee.FeatureCollection([
  ee.Feature(ee.Geometry.Point([108.0500, 12.6667]), {Station: 'Buon Ma Thuot', Region: 'Trung tam cao nguyen'}),
  ee.Feature(ee.Geometry.Point([107.8833, 13.0833]), {Station: 'Ea Sup',        Region: 'Tay Bac (rung khop)'}),
  ee.Feature(ee.Geometry.Point([107.8000, 12.8833]), {Station: 'Buon Don',      Region: 'Tay Bac (bien gioi)'}),
  ee.Feature(ee.Geometry.Point([108.2000, 13.2167]), {Station: 'Ea Hleo',       Region: 'Bac cao nguyen'}),
  ee.Feature(ee.Geometry.Point([108.7667, 12.7500]), {Station: 'MDrak',         Region: 'Dong cao nguyen'}),
  ee.Feature(ee.Geometry.Point([108.1667, 12.4167]), {Station: 'Lak',           Region: 'Nam cao nguyen'}),
  ee.Feature(ee.Geometry.Point([108.3167, 12.6833]), {Station: 'Krong Pac',     Region: 'Trung tam nong nghiep'}),
  ee.Feature(ee.Geometry.Point([109.2950, 13.0955]), {Station: 'Tuy Hoa',       Region: 'Duyen hai phia Dong'}),
  ee.Feature(ee.Geometry.Point([108.9667, 13.0500]), {Station: 'Son Hoa',       Region: 'Thung lung song Ba'}),
  ee.Feature(ee.Geometry.Point([109.0333, 13.3667]), {Station: 'Dong Xuan',     Region: 'Mien nui Dong Bac'})
]);

Map.addLayer(stations, {color: 'blue'}, '10 Tram / Dia danh dai dien');

// --- 4. VÒNG LẶP XỬ LÝ 3 MỐC NĂM (2016, 2021, 2026) ---
var years = [2016, 2021, 2026];
var meteoFeaturesAll = ee.FeatureCollection([]);

years.forEach(function(year) {
  var startDate = year + '-02-01';
  var endDate   = year + '-05-01'; // Tới hết ngày 30/04

  // Lọc Landsat 8 và Landsat 9 (Landsat 9 phóng tháng 9/2021 nên có ảnh từ 2022 trở đi)
  var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
             .filterBounds(roi)
             .filterDate(startDate, endDate);
  var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
             .filterBounds(roi)
             .filterDate(startDate, endDate);

  var merged = l8.merge(l9);

  // In thông tin số lượng cảnh ảnh ra Console (Bổ sung vào Bảng tiền xử lý Chương 2)
  print('=== NĂM ' + year + ' (01/02 - 30/04) ===');
  print('  Số cảnh Landsat 8:', l8.size());
  print('  Số cảnh Landsat 9:', l9.size());
  print('  Tổng số cảnh đầu vào:', merged.size());

  // Tiền xử lý và tổng hợp Median mùa khô
  var processed = merged.map(preprocessLandsat);
  var composite = processed.median().clip(roi);

  var ndviImg = composite.select('NDVI').toFloat();
  var lstImg  = composite.select('LST').toFloat();

  // Tính nhanh tỷ lệ che phủ hợp lệ (%) trên vùng nghiên cứu
  var validCount = ndviImg.mask().reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: roi,
    scale: 120,
    maxPixels: 1e9
  }).get('NDVI');
  print('  Tỷ lệ điểm ảnh hợp lệ (sau lọc mây & ST_QA):', ee.Number(validCount).multiply(100), '%');

  // Xuất NDVI và LST ở hệ tọa độ phẳng UTM Zone 49N (EPSG:32649), scale 60m
  Export.image.toDrive({
    image: ndviImg,
    description: 'NDVI_' + year + '_DakLak_UTM49N',
    folder: DRIVE_FOLDER,
    fileNamePrefix: 'NDVI_' + year + '_DakLak_UTM49N',
    region: roi,
    scale: TARGET_SCALE,
    crs: TARGET_CRS,
    maxPixels: 1e10
  });

  Export.image.toDrive({
    image: lstImg,
    description: 'LST_' + year + '_DakLak_UTM49N',
    folder: DRIVE_FOLDER,
    fileNamePrefix: 'LST_' + year + '_DakLak_UTM49N',
    region: roi,
    scale: TARGET_SCALE,
    crs: TARGET_CRS,
    maxPixels: 1e10
  });

  // --- 5. TRÍCH XUẤT KHÍ TƯỢNG ĐỐI CHIẾU (CHIRPS MƯA + ERA5-LAND NHIỆT ĐỘ & ĐỘ ẨM ĐẤT) ---
  var chirpsSum = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY')
                    .filterBounds(roi)
                    .filterDate(startDate, endDate)
                    .select('precipitation')
                    .sum()
                    .rename('Rain_CHIRPS_mm');

  var era5 = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR')
               .filterBounds(roi)
               .filterDate(startDate, endDate);

  var temp2m = era5.select('temperature_2m').mean().subtract(273.15).rename('Temp2m_ERA5_C');
  var soilM  = era5.select('volumetric_soil_water_layer_1').mean().rename('SoilMoisture_ERA5');

  var stackForValidation = ndviImg.addBands([lstImg, chirpsSum, temp2m, soilM]);

  // Lấy giá trị trung bình trong vùng đệm bán kính 3km quanh mỗi trạm/điểm đại diện
  var stationStats = stations.map(function(feat) {
    var buf = feat.geometry().buffer(3000);
    var stats = stackForValidation.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: buf,
      scale: 1000,
      maxPixels: 1e8
    });
    return feat.set(stats).set({
      'Year': year,
      'Lon': feat.geometry().coordinates().get(0),
      'Lat': feat.geometry().coordinates().get(1)
    });
  });

  meteoFeaturesAll = meteoFeaturesAll.merge(stationStats);
});

// Xuất 1 bảng CSV duy nhất chứa toàn bộ dữ liệu khí tượng đối chiếu 3 năm (30 dòng)
Export.table.toDrive({
  collection: meteoFeaturesAll,
  description: 'KhiTuong_DoiChieu_CHIRPS_ERA5_DakLak',
  folder: DRIVE_FOLDER,
  fileNamePrefix: 'KhiTuong_DoiChieu_CHIRPS_ERA5_DakLak',
  fileFormat: 'CSV',
  selectors: ['Year', 'Station', 'Region', 'Lon', 'Lat', 'Rain_CHIRPS_mm', 'Temp2m_ERA5_C', 'SoilMoisture_ERA5', 'NDVI', 'LST']
});

// --- 6. XUẤT LỚP ĐỊA HÌNH (SRTM DEM) & LỚP PHỦ ĐẤT (ESA WORLDCOVER) CÙNG UTM 49N ---
var dem = ee.Image('USGS/SRTMGL1_003').select('elevation').clip(roi).toInt16();
Export.image.toDrive({
  image: dem,
  description: 'DEM_SRTM_DakLak_UTM49N',
  folder: DRIVE_FOLDER,
  fileNamePrefix: 'DEM_SRTM_DakLak_UTM49N',
  region: roi,
  scale: TARGET_SCALE,
  crs: TARGET_CRS,
  maxPixels: 1e10
});

var worldcover = ee.ImageCollection('ESA/WorldCover/v200').first().select('Map').clip(roi).toUint8();
Export.image.toDrive({
  image: worldcover,
  description: 'LandCover_ESA_DakLak_UTM49N',
  folder: DRIVE_FOLDER,
  fileNamePrefix: 'LandCover_ESA_DakLak_UTM49N',
  region: roi,
  scale: TARGET_SCALE,
  crs: TARGET_CRS,
  maxPixels: 1e10
});

// --- 7. XUẤT RANH GIỚI TỈNH & TỌA ĐỘ 10 ĐỊA DANH CHO QGIS ---
Export.table.toDrive({
  collection: ee.FeatureCollection([ee.Feature(roi, {Name: 'Dak Lak Moi'})]),
  description: 'RanhGioi_DakLakMoi_SHP',
  folder: DRIVE_FOLDER,
  fileNamePrefix: 'RanhGioi_DakLakMoi_SHP',
  fileFormat: 'SHP'
});

Export.table.toDrive({
  collection: stations.map(function(f) {
    return f.set({
      'Lon': f.geometry().coordinates().get(0),
      'Lat': f.geometry().coordinates().get(1)
    });
  }),
  description: 'DiaDanh_DakLakMoi_QGIS',
  folder: DRIVE_FOLDER,
  fileNamePrefix: 'DiaDanh_DakLakMoi_QGIS',
  fileFormat: 'CSV',
  selectors: ['Station', 'Region', 'Lon', 'Lat']
});

print('>>> HOÀN TẤT THIẾT LẬP! Vui lòng chuyển sang tab Tasks bên phải và bấm RUN tất cả 11 tác vụ.');