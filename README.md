# Đánh giá hạn hán từ ảnh vệ tinh sử dụng NDVI và Nhiệt độ bề mặt (TVDI)

Đồ án Chuyên ngành Tin học Môi trường — Khoa Công nghệ Thông tin, Trường Đại học Tài nguyên và Môi trường Thành phố Hồ Chí Minh (HCMUNRE).

- **Sinh viên thực hiện:** Trương Thành Danh (MSSV: `1250080028` · Lớp: `12_ĐH_THMT` · Khóa: `2023 – 2027`)
- **Giảng viên hướng dẫn:** ThS. Từ Thanh Trí
- **Khu vực nghiên cứu:** Tỉnh Đắk Lắk mới (hợp nhất Đắk Lắk và Phú Yên từ 01/07/2025, diện tích tự nhiên ~18.096,40 km²)
- **Mốc thời gian:** Cao điểm mùa khô (`01/02 – 30/04`) các năm `2016`, `2021` và `2026`

---

## 1. Cấu trúc thư mục dự án

```text
DACN_TruongThanhDanh_1250080028/
├── 01_Code_Scripts/
│   ├── 01_gee_codeeditor_daklak_v2.js        # Tiền xử lý Landsat 8/9, lọc mây, xuất dữ liệu trên GEE
│   ├── 02_PhanTich_TVDI_Va_BieuDo.py         # Hồi quy OLS (Dry/Wet Edge), TVDI, vùng hợp lệ chung, kiểm chứng
│   ├── 03_Ve_Scatter_LST_NDVI.py             # Vẽ 3 đồ thị không gian LST–NDVI chuẩn in trắng đen
│   └── 04_Ve_Boxplot_TVDI.py                 # Vẽ biểu đồ hộp (Boxplot) phân bố TVDI qua 3 năm
├── 02_DuLieu_GeoTIFF/
│   └── V2_UTM49N/                            # Ảnh GeoTIFF (EPSG:32649, 60 m), DEM, ESA WorldCover, CSV khí tượng
│       └── KetQua_V2/                        # Bảng thống kê Bang_1..10 (.csv), ảnh Delta_TVDI, biểu đồ (.png)
├── 03_BanDo_QGIS/                            # Dự án BanDo_TVDI_DakLak_V2.qgz, 3 ảnh TVDI (.tif) và 6 bản đồ (.png)
├── README.md                                 # Hướng dẫn tổng quan, truy xuất dữ liệu và vận hành
├── requirements.txt                          # Danh sách thư viện Python cần thiết
└── .gitignore                                # Loại trừ ảnh GeoTIFF dung lượng lớn khỏi Git
```

---

## 2. Môi trường thực thi & Cài đặt thư viện

- **Hệ điều hành:** Windows 10/11 (64-bit)
- **Google Earth Engine (GEE):** Tài khoản [Google Earth Engine Code Editor](https://code.earthengine.google.com/)
- **Python:** Python `3.10+` (đã kiểm thử trên Python `3.14`)
- **QGIS:** QGIS `3.34 LTR` trở lên

Cài đặt toàn bộ thư viện Python từ file `requirements.txt`:

```bash
pip install -r requirements.txt
```

---

## 3. Dữ liệu (Dataset) & Hướng dẫn truy xuất

Do các tệp raster GeoTIFF có dung lượng lớn (~150 MB), bộ dữ liệu đầu vào (`02_DuLieu_GeoTIFF/V2_UTM49N/`) và dự án QGIS (`03_BanDo_QGIS/`) được lưu trữ kèm theo trên Google Drive:

- **Link tải trọn bộ Dataset & Dự án QGIS:** [Google Drive — Hồ sơ & Dữ liệu ĐACN (1250080028)](https://drive.google.com/drive/folders/1l9maQTkbSHpKw8UxGsa7vVWiKZlA5YOF?usp=sharing)

Hoặc tái lập dữ liệu gốc trực tiếp từ Google Earth Engine:
1. Mở [Google Earth Engine Code Editor](https://code.earthengine.google.com/).
2. Dán nội dung file `01_Code_Scripts/01_gee_codeeditor_daklak_v2.js` và bấm **Run**.
3. Chuyển sang tab **Tasks**, chạy các tác vụ xuất để tải về Google Drive:
   - 6 ảnh `NDVI_*_DakLak_UTM49N.tif` và `LST_*_DakLak_UTM49N.tif` (Landsat 8/9 Collection 2 Level 2, lọc `QA_PIXEL` Bit 3, 4 = 0 và `ST_QA ≤ 4,0 K`, chiếu phẳng `EPSG:32649`, độ phân giải `60 m`).
   - Ảnh địa hình `DEM_SRTM_DakLak_UTM49N.tif` (SRTM 30 m) và lớp phủ `LandCover_ESA_DakLak_UTM49N.tif` (ESA WorldCover 2021 v200).
   - Bảng khí tượng `KhiTuong_DoiChieu_CHIRPS_ERA5_DakLak.csv` (lượng mưa CHIRPS & nhiệt độ 2m ERA5-Land tại 10 khu vực đại diện).
   - Bộ ranh giới hành chính tỉnh/huyện Đắk Lắk mới (`GeoJSON` / `Shapefile`).

---

## 4. Hướng dẫn chạy hệ thống

### Bước 1: Chạy phân tích TVDI, thống kê vùng hợp lệ chung và kiểm chứng khí tượng
1. Mở file `01_Code_Scripts/02_PhanTich_TVDI_Va_BieuDo.py` và kiểm tra biến `DATA_DIR` trỏ đúng tới thư mục `02_DuLieu_GeoTIFF\V2_UTM49N` trên máy.
2. Chạy lệnh:
   ```bash
   python 01_Code_Scripts/02_PhanTich_TVDI_Va_BieuDo.py
   ```
   Toàn bộ 9 bảng kết quả (`Bang_1` đến `Bang_9`), 5 ảnh GeoTIFF (`TVDI` 3 năm + `Delta_TVDI` 2 giai đoạn) và các biểu đồ sẽ tự động xuất ra thư mục `KetQua_V2/` trong khoảng ~43 giây.

### Bước 2: Xuất đồ thị phân tán LST–NDVI và biểu đồ hộp Boxplot
```bash
python 01_Code_Scripts/03_Ve_Scatter_LST_NDVI.py
python 01_Code_Scripts/04_Ve_Boxplot_TVDI.py
```

### Bước 3: Mở bản đồ chuyên đề trên QGIS
- Mở file `03_BanDo_QGIS/BanDo_TVDI_DakLak_V2.qgz` bằng phần mềm QGIS để xem các lớp TVDI, ΔTVDI và bố cục in ấn trong **Print Layout** (300 DPI).

---

## 5. Tóm tắt kết quả chính

| Chỉ tiêu | Năm 2016 (El Niño) | Năm 2021 | Năm 2026 |
| :--- | :---: | :---: | :---: |
| Số cảnh Landsat 8/9 hợp lệ | 33 (95,23%) | 32 (98,51%) | 60 (99,16%) |
| Đường biên khô Dry Edge (`R²` / `RMSE`) | `0,658` / `3,30°C` | `0,790` / `2,14°C` | `0,695` / `2,44°C` |
| Tỷ lệ **Hạn nặng** (trên 13.999,74 km²) | `29,55%` (4.136,30 km²) | `22,88%` (3.203,22 km²) | `21,29%` (2.980,74 km²) |
| Tỷ lệ **Ẩm ướt** (trên 13.999,74 km²) | `11,06%` (1.547,84 km²) | `15,65%` (2.190,70 km²) | `19,62%` (2.746,27 km²) |

- **Kiểm chứng khí tượng (n = 30, TVDI biên gộp):** tương quan nghịch chặt với lượng mưa mùa khô CHIRPS (`r = −0,689; R² = 0,475; p < 0,0001`) và tương quan thuận chặt với nhiệt độ 2m ERA5-Land (`r = +0,661; R² = 0,437; p < 0,0001`).
