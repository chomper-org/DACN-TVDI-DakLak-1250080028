# =============================================================================
# ĐỒ ÁN CHUYÊN NGÀNH — CHẶNG 2: PHÂN TÍCH TOÀN DIỆN & KIỂM CHỨNG TVDI (V2)
# Sinh viên: Trương Thành Danh (1250080028) — Lớp 12_THMT
#
# HƯỚNG DẪN CHẠY:
# 1. Tải toàn bộ 11 tệp từ thư mục Google Drive 'DACN_DakLak_TVDI_V2' về máy,
#    đặt chung vào 1 thư mục (ví dụ: D:\DACN_TruongThanhDanh_1250080028\Data_V2).
# 2. Chỉnh đường dẫn DATA_DIR bên dưới đúng với thư mục trên máy bạn.
# 3. Chạy lệnh: python 04_PhanTich_ToanDien_TVDI_V2.py
# =============================================================================

import os
import warnings
import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import reproject, Resampling, transform as crs_transform
from scipy import stats
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# >>> CHỈNH ĐƯỜNG DẪN THƯ MỤC CHỨA DỮ LIỆU TẢI TỪ GEE TẠI ĐÂY <<<
DATA_DIR = r"D:\DACN_TruongThanhDanh_1250080028\02_DuLieu_GeoTIFF\V2_UTM49N"
OUT_DIR  = os.path.join(DATA_DIR, "KetQua_V2")
os.makedirs(OUT_DIR, exist_ok=True)

YEARS = [2016, 2021, 2026]
NODATA_OUT = -9999.0

# -----------------------------------------------------------------------------
# HÀM HỖ TRỢ ĐỌC & CĂN CHỈNH LƯỚI RASTER ĐỒNG NHẤT (EPSG:32649, 60x60m)
# -----------------------------------------------------------------------------
def load_and_align(path, ref_profile=None, resampling=Resampling.bilinear):
    with rasterio.open(path) as src:
        if ref_profile is None:
            arr = src.read(1).astype(np.float32)
            nodata = src.nodata
            if nodata is not None:
                arr[arr == nodata] = np.nan
            return arr, src.profile.copy()

        # Căn chỉnh về đúng kích thước & transform của ảnh chuẩn
        dst = np.full((ref_profile["height"], ref_profile["width"]), np.nan, dtype=np.float32)
        src_arr = src.read(1).astype(np.float32)
        if src.nodata is not None:
            src_arr[src_arr == src.nodata] = np.nan
        reproject(
            source=src_arr,
            destination=dst,
            src_transform=src.transform,
            src_crs=src.crs,
            src_nodata=np.nan,
            dst_transform=ref_profile["transform"],
            dst_crs=ref_profile["crs"],
            dst_nodata=np.nan,
            resampling=resampling,
        )
        return dst, ref_profile

def save_geotiff(arr, profile, out_path):
    p = profile.copy()
    p.update(dtype="float32", count=1, nodata=NODATA_OUT, compress="lzw")
    out_arr = np.where(np.isnan(arr), NODATA_OUT, arr).astype(np.float32)
    with rasterio.open(out_path, "w", **p) as dst:
        dst.write(out_arr, 1)

def fit_dry_wet_edges(ndvi_vals, lst_vals, p_high=95, p_low=5, bin_step=0.01):
    bins = np.arange(0.10, 0.85 + bin_step/2, bin_step)
    dry_x, dry_y, wet_x, wet_y = [], [], [], []
    for i in range(len(bins) - 1):
        m = (ndvi_vals >= bins[i]) & (ndvi_vals < bins[i + 1])
        if np.sum(m) >= 50:
            mid = (bins[i] + bins[i + 1]) / 2.0
            dry_x.append(mid); dry_y.append(np.percentile(lst_vals[m], p_high))
            wet_x.append(mid); wet_y.append(np.percentile(lst_vals[m], p_low))
    dry_x, dry_y = np.array(dry_x), np.array(dry_y)
    wet_x, wet_y = np.array(wet_x), np.array(wet_y)

    rd = stats.linregress(dry_x, dry_y)
    rw = stats.linregress(wet_x, wet_y)

    dry_rmse = float(np.sqrt(np.mean((dry_y - (rd.slope * dry_x + rd.intercept)) ** 2)))
    wet_rmse = float(np.sqrt(np.mean((wet_y - (rw.slope * wet_x + rw.intercept)) ** 2)))

    return {
        "a": rd.slope, "b": rd.intercept, "dry_r2": rd.rvalue ** 2, "dry_rmse": dry_rmse, "dry_p": rd.pvalue,
        "c": rw.slope, "d": rw.intercept, "wet_r2": rw.rvalue ** 2, "wet_rmse": wet_rmse, "wet_p": rw.pvalue,
        "dry_x": dry_x, "dry_y": dry_y, "wet_x": wet_x, "wet_y": wet_y, "n_bins": len(dry_x)
    }

def compute_tvdi_matrix(ndvi, lst, a, b, c, d):
    valid = (~np.isnan(ndvi)) & (~np.isnan(lst)) & (ndvi >= 0.10) & (ndvi <= 0.85) & (lst > 0) & (lst < 60)
    tvdi = np.full(ndvi.shape, np.nan, dtype=np.float32)
    lst_max = a * ndvi[valid] + b
    lst_min = c * ndvi[valid] + d
    denom = lst_max - lst_min
    good = denom > 0.5
    vals = np.full(np.sum(valid), np.nan, dtype=np.float32)
    vals[good] = np.clip((lst[valid][good] - lst_min[good]) / denom[good], 0.0, 1.0)
    tvdi[valid] = vals
    return tvdi

def classify_areas(tvdi_arr, pixel_area_km2, mask=None):
    if mask is not None:
        v = tvdi_arr[mask & (~np.isnan(tvdi_arr))]
    else:
        v = tvdi_arr[~np.isnan(tvdi_arr)]
    total_px = len(v)
    total_km2 = total_px * pixel_area_km2
    classes = [
        ("Ẩm ướt (0,0 - 0,2)",     (v >= 0.0) & (v < 0.2)),
        ("Bình thường (0,2 - 0,4)", (v >= 0.2) & (v < 0.4)),
        ("Hạn nhẹ (0,4 - 0,6)",     (v >= 0.4) & (v < 0.6)),
        ("Hạn vừa (0,6 - 0,8)",     (v >= 0.6) & (v < 0.8)),
        ("Hạn nặng (0,8 - 1,0)",    (v >= 0.8) & (v <= 1.0)),
    ]
    rows = []
    for label, cond in classes:
        cnt = int(np.sum(cond))
        km2 = cnt * pixel_area_km2
        pct = (cnt / total_px * 100.0) if total_px > 0 else 0.0
        rows.append((label, cnt, round(km2, 2), round(pct, 2)))
    return rows, total_px, round(total_km2, 2)

# =============================================================================
# BƯỚC 1: NẠP DỮ LIỆU RASTER (EPSG:32649) & KIỂM TRA KÍCH THƯỚC PIXEL CHUẨN
# =============================================================================
print("1. Đang nạp và đồng bộ lưới ảnh UTM Zone 49N (EPSG:32649)...")
ndvi_dict, lst_dict = {}, {}
ref_profile = None

for y in YEARS:
    p_ndvi = os.path.join(DATA_DIR, f"NDVI_{y}_DakLak_UTM49N.tif")
    p_lst  = os.path.join(DATA_DIR, f"LST_{y}_DakLak_UTM49N.tif")
    ndvi_dict[y], ref_profile = load_and_align(p_ndvi, ref_profile, Resampling.bilinear)
    lst_dict[y], _            = load_and_align(p_lst,  ref_profile, Resampling.bilinear)

pixel_w = abs(ref_profile["transform"].a)
pixel_h = abs(ref_profile["transform"].e)
PIXEL_AREA_KM2 = (pixel_w * pixel_h) / 1e6
print(f"   CRS: {ref_profile['crs']} | Kích thước pixel: {pixel_w:.2f}m x {pixel_h:.2f}m = {PIXEL_AREA_KM2:.6f} km2/pixel")

# =============================================================================
# BƯỚC 2: HỒI QUY DRY/WET EDGE + R2/RMSE/P-VALUE + ĐỘ NHẠY PHÂN VỊ (LỖI #07)
# =============================================================================
print("\n2. Đang chạy hồi quy OLS, tính R², RMSE, p-value và phân tích độ nhạy...")
reg_rows = []
sens_rows = []
tvdi_dict = {}
pooled_ndvi_samples, pooled_lst_samples = [], []
rng = np.random.default_rng(42)

for y in YEARS:
    ndvi, lst = ndvi_dict[y], lst_dict[y]
    m_valid = (~np.isnan(ndvi)) & (~np.isnan(lst)) & (ndvi >= 0.10) & (ndvi <= 0.85) & (lst > 0) & (lst < 60)
    nv, lv = ndvi[m_valid], lst[m_valid]

    # Lấy mẫu 200.000 điểm mỗi năm cho bài toán kiểm chứng biên gộp chung (#02)
    idx_sample = rng.choice(len(nv), size=min(200000, len(nv)), replace=False)
    pooled_ndvi_samples.append(nv[idx_sample])
    pooled_lst_samples.append(lv[idx_sample])

    # Hồi quy chuẩn 95% / 5%
    res = fit_dry_wet_edges(nv, lv, p_high=95, p_low=5)
    tvdi = compute_tvdi_matrix(ndvi, lst, res["a"], res["b"], res["c"], res["d"])
    tvdi_dict[y] = tvdi

    # Lưu GeoTIFF TVDI chuẩn UTM 49N cho QGIS
    save_geotiff(tvdi, ref_profile, os.path.join(OUT_DIR, f"TVDI_{y}_UTM49N.tif"))

    reg_rows.append({
        "Nam": y, "So_Bin": res["n_bins"],
        "Dry_Slope_a": round(res["a"], 4), "Dry_Intercept_b": round(res["b"], 2),
        "Dry_R2": round(res["dry_r2"], 4), "Dry_RMSE_C": round(res["dry_rmse"], 2), "Dry_p_value": f"{res['dry_p']:.2e}",
        "Wet_Slope_c": round(res["c"], 4), "Wet_Intercept_d": round(res["d"], 2),
        "Wet_R2": round(res["wet_r2"], 4), "Wet_RMSE_C": round(res["wet_rmse"], 2), "Wet_p_value": f"{res['wet_p']:.2e}",
    })

    # Vẽ Scatter Plot có kèm R² & RMSE
    fig, ax = plt.subplots(figsize=(7.5, 5.2), dpi=200)
    ax.scatter(nv[idx_sample[:35000]], lv[idx_sample[:35000]], s=1, c="#94A3B8", alpha=0.25, rasterized=True)
    ax.scatter(res["dry_x"], res["dry_y"], s=14, c="#DC2626", label="Phân vị 95% (Biên khô)")
    ax.scatter(res["wet_x"], res["wet_y"], s=14, c="#2563EB", label="Phân vị 5% (Biên ẩm)")
    x_line = np.linspace(0.10, 0.85, 100)
    ax.plot(x_line, res["a"]*x_line + res["b"], "r-", lw=2,
            label=f"Dry: y={res['a']:.2f}x+{res['b']:.2f} (R²={res['dry_r2']:.3f}, RMSE={res['dry_rmse']:.2f}°C)")
    ax.plot(x_line, res["c"]*x_line + res["d"], "b-", lw=2,
            label=f"Wet: y={res['c']:.2f}x+{res['d']:.2f} (R²={res['wet_r2']:.3f}, RMSE={res['wet_rmse']:.2f}°C)")
    ax.set_title(f"Không gian đặc trưng LST–NDVI mùa khô năm {y} (Đắk Lắk mới)", fontsize=12, fontweight="bold")
    ax.set_xlabel("Chỉ số thực vật NDVI"); ax.set_ylabel("Nhiệt độ bề mặt LST (°C)")
    ax.set_xlim(0.05, 0.90); ax.set_ylim(15, 62)
    ax.grid(True, linestyle="--", alpha=0.35)
    ax.legend(fontsize=8.5, loc="upper right")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, f"Scatter_LST_NDVI_{y}_V2.png"))
    plt.close()

    # Kiểm tra độ nhạy phân vị (98/2 và 90/10 so với chuẩn 95/5)
    base_severe = np.nanmean(tvdi >= 0.8) * 100.0
    for ph, pl in [(98, 2), (90, 10)]:
        r_alt = fit_dry_wet_edges(nv, lv, p_high=ph, p_low=pl)
        t_alt = compute_tvdi_matrix(ndvi, lst, r_alt["a"], r_alt["b"], r_alt["c"], r_alt["d"])
        both = (~np.isnan(tvdi)) & (~np.isnan(t_alt))
        mae = float(np.mean(np.abs(tvdi[both] - t_alt[both])))
        corr_r = float(np.corrcoef(tvdi[both], t_alt[both])[0, 1])
        alt_severe = np.nanmean(t_alt >= 0.8) * 100.0
        sens_rows.append({
            "Nam": y, "Kich_Ban": f"P{ph}/P{pl} vs P95/P5",
            "Dry_b_C": round(r_alt["b"], 2), "Wet_d_C": round(r_alt["d"], 2),
            "MAE_TVDI": round(mae, 4), "Pearson_r": round(corr_r, 4),
            "HanNang_Base_%": round(base_severe, 2), "HanNang_Alt_%": round(alt_severe, 2),
            "ChenhLech_HanNang_%": round(alt_severe - base_severe, 2)
        })

df_reg = pd.DataFrame(reg_rows)
df_reg.to_csv(os.path.join(OUT_DIR, "Bang_1_HoiQuy_DryWet_ChatLuong.csv"), index=False, encoding="utf-8-sig")
df_sens = pd.DataFrame(sens_rows)
df_sens.to_csv(os.path.join(OUT_DIR, "Bang_2_PhanTich_DoNhay_PhanVi.csv"), index=False, encoding="utf-8-sig")
print(df_reg[["Nam", "Dry_Slope_a", "Dry_Intercept_b", "Dry_R2", "Dry_RMSE_C", "Wet_R2", "Wet_RMSE_C"]].to_string(index=False))

# =============================================================================
# BƯỚC 3: VÙNG HỢP LỆ CHUNG (LỖI #03) & SO SÁNH HỆ QUY CHIẾU BIÊN GỘP (LỖI #02)
# =============================================================================
print("\n3. Đang thống kê trên vùng hợp lệ chung và kiểm chứng biên gộp chung 3 năm...")
common_mask = (~np.isnan(tvdi_dict[2016])) & (~np.isnan(tvdi_dict[2021])) & (~np.isnan(tvdi_dict[2026]))
common_km2 = round(int(np.sum(common_mask)) * PIXEL_AREA_KM2, 2)
print(f"   >>> Diện tích VÙNG HỢP LỆ CHUNG đồng nhất cả 3 năm: {common_km2:,.2f} km²")

# Xuất raster biến động ΔTVDI trên vùng hợp lệ chung cho QGIS (Lỗi #08)
diff_21_16 = np.where(common_mask, tvdi_dict[2021] - tvdi_dict[2016], np.nan)
diff_26_16 = np.where(common_mask, tvdi_dict[2026] - tvdi_dict[2016], np.nan)
save_geotiff(diff_21_16, ref_profile, os.path.join(OUT_DIR, "Delta_TVDI_2021_minus_2016_UTM49N.tif"))
save_geotiff(diff_26_16, ref_profile, os.path.join(OUT_DIR, "Delta_TVDI_2026_minus_2016_UTM49N.tif"))

# Tính đường biên gộp chung 3 năm (Pooled Dry/Wet Edge)
pool_n = np.concatenate(pooled_ndvi_samples)
pool_l = np.concatenate(pooled_lst_samples)
res_pool = fit_dry_wet_edges(pool_n, pool_l, p_high=95, p_low=5)
print(f"   >>> Biên gộp chung 3 năm: Dry = {res_pool['a']:.2f}*NDVI + {res_pool['b']:.2f} (R²={res_pool['dry_r2']:.3f}) | "
      f"Wet = {res_pool['c']:.2f}*NDVI + {res_pool['d']:.2f} (R²={res_pool['wet_r2']:.3f})")

rows_indiv, rows_common, rows_ref = [], [], []
tvdi_pooled_dict = {}

for y in YEARS:
    # A. Toàn vùng hợp lệ riêng của năm (UTM 49N)
    cls_i, px_i, km2_i = classify_areas(tvdi_dict[y], PIXEL_AREA_KM2, mask=None)
    for label, cnt, km2, pct in cls_i:
        rows_indiv.append({"Nam": y, "Cap_Han": label, "So_Pixel": cnt, "DienTich_km2": km2, "TyLe_%": pct, "TongHopLe_km2": km2_i})

    # B. Vùng hợp lệ chung cả 3 năm (UTM 49N)
    cls_c, px_c, km2_c = classify_areas(tvdi_dict[y], PIXEL_AREA_KM2, mask=common_mask)
    for label, cnt, km2, pct in cls_c:
        rows_common.append({"Nam": y, "Cap_Han": label, "So_Pixel": cnt, "DienTich_km2": km2, "TyLe_%": pct, "TongHopLeChung_km2": km2_c})

    # C. So sánh hệ quy chiếu TVDI: Biên riêng vs. Biên gộp chung trên vùng hợp lệ chung
    t_pool = compute_tvdi_matrix(ndvi_dict[y], lst_dict[y], res_pool["a"], res_pool["b"], res_pool["c"], res_pool["d"])
    tvdi_pooled_dict[y] = t_pool
    m_c = common_mask & (~np.isnan(t_pool))
    rows_ref.append({
        "Nam": y,
        "TVDI_BienRieng_TB": round(float(np.mean(tvdi_dict[y][m_c])), 4),
        "HanNang_BienRieng_%": round(float(np.mean(tvdi_dict[y][m_c] >= 0.8) * 100), 2),
        "TVDI_BienGop_TB": round(float(np.mean(t_pool[m_c])), 4),
        "HanNang_BienGop_%": round(float(np.mean(t_pool[m_c] >= 0.8) * 100), 2),
        "MAE_Giua_2_HeQuyChieu": round(float(np.mean(np.abs(tvdi_dict[y][m_c] - t_pool[m_c]))), 4),
        "Pearson_r_2_HeQuyChieu": round(float(np.corrcoef(tvdi_dict[y][m_c], t_pool[m_c])[0, 1]), 4)
    })

df_indiv  = pd.DataFrame(rows_indiv)
df_common = pd.DataFrame(rows_common)
df_ref    = pd.DataFrame(rows_ref)
df_indiv.to_csv(os.path.join(OUT_DIR, "Bang_4_DienTich_UTM49N_TungNam.csv"), index=False, encoding="utf-8-sig")
df_common.to_csv(os.path.join(OUT_DIR, "Bang_5_DienTich_VungHopLeChung.csv"), index=False, encoding="utf-8-sig")
df_ref.to_csv(os.path.join(OUT_DIR, "Bang_3_SoSanh_HeQuyChieu_TVDI.csv"), index=False, encoding="utf-8-sig")
print("\n   So sánh hệ quy chiếu TVDI (Biên riêng vs. Biên gộp chung trên vùng hợp lệ chung):")
print(df_ref.to_string(index=False))

# Vẽ biểu đồ cột biến động trên vùng hợp lệ chung (chuẩn in trắng đen)
plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"]})
piv = df_common.pivot(index="Cap_Han", columns="Nam", values="TyLe_%")
order = ["Ẩm ướt (0,0 - 0,2)", "Bình thường (0,2 - 0,4)", "Hạn nhẹ (0,4 - 0,6)", "Hạn vừa (0,6 - 0,8)", "Hạn nặng (0,8 - 1,0)"]
short_labels = ["Ẩm ướt\n(0,0–0,2)", "Bình thường\n(0,2–0,4)", "Hạn nhẹ\n(0,4–0,6)", "Hạn vừa\n(0,6–0,8)", "Hạn nặng\n(0,8–1,0)"]
piv = piv.reindex(order)
fig, ax = plt.subplots(figsize=(9.5, 5.2), dpi=300)
x = np.arange(len(order))
w = 0.25
bar_styles = [
    dict(color="#2C3643", edgecolor="#111111", hatch=""),
    dict(color="#9BA4B5", edgecolor="#111111", hatch="///"),
    dict(color="#FFFFFF", edgecolor="#111111", hatch="..."),
]
for idx, y in enumerate(YEARS):
    bars = ax.bar(x + (idx - 1)*w, piv[y].values, width=w, label=f"Năm {y}", linewidth=0.9, **bar_styles[idx])
    for b in bars:
        ax.text(b.get_x() + b.get_width()/2, b.get_height() + 0.4, f"{b.get_height():.1f}%",
                ha="center", va="bottom", fontsize=9)
ax.set_xticks(x); ax.set_xticklabels(short_labels, fontsize=10)
ax.set_ylabel("Tỷ lệ diện tích (%) trên vùng hợp lệ chung", fontsize=11)
ax.set_ylim(0, max(piv.max()) + 5.5)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ax.legend(loc="upper left", ncol=3, frameon=False, fontsize=10.5)
ax.grid(axis="y", linestyle=":", alpha=0.45)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "BieuDo_BienDong_VungHopLeChung.png"), dpi=300)
plt.close()

# =============================================================================
# BƯỚC 4: KIỂM CHỨNG ĐỊNH LƯỢNG KHÍ TƯỢNG CHIRPS & ERA5-LAND (LỖI #01)
# =============================================================================
meteo_path = os.path.join(DATA_DIR, "KhiTuong_DoiChieu_CHIRPS_ERA5_DakLak.csv")
if os.path.exists(meteo_path):
    print("\n4. Đang chạy kiểm chứng định lượng TVDI với dữ liệu khí tượng CHIRPS & ERA5-Land...")
    df_m = pd.read_csv(meteo_path)
    # Chuyển tọa độ (Lon, Lat) WGS84 sang UTM 49N để trích TVDI quanh bán kính 3km (~50 pixel)
    xs, ys = crs_transform("EPSG:4326", ref_profile["crs"], df_m["Lon"].tolist(), df_m["Lat"].tolist())
    inv_tf = ~ref_profile["transform"]
    rad_px = int(round(3000.0 / pixel_w))

    tvdi_st_indiv, tvdi_st_pool = [], []
    for i, row in df_m.iterrows():
        y = int(row["Year"])
        col, r = inv_tf * (xs[i], ys[i])
        col, r = int(round(col)), int(round(r))
        r0, r1 = max(0, r - rad_px), min(ref_profile["height"], r + rad_px + 1)
        c0, c1 = max(0, col - rad_px), min(ref_profile["width"], col + rad_px + 1)
        patch_i = tvdi_dict[y][r0:r1, c0:c1]
        patch_p = tvdi_pooled_dict[y][r0:r1, c0:c1]
        tvdi_st_indiv.append(round(float(np.nanmean(patch_i)), 4) if np.any(~np.isnan(patch_i)) else np.nan)
        tvdi_st_pool.append(round(float(np.nanmean(patch_p)), 4) if np.any(~np.isnan(patch_p)) else np.nan)

    df_m["TVDI_BienRieng"] = tvdi_st_indiv
    df_m["TVDI_BienGop"]   = tvdi_st_pool
    df_m.to_csv(os.path.join(OUT_DIR, "Bang_6_KiemChung_KhiTuong_30Diem.csv"), index=False, encoding="utf-8-sig")

    corr_rows = []
    vars_check = [
        ("Rain_CHIRPS_mm", "Lượng mưa mùa khô CHIRPS (mm)"),
        ("Temp2m_ERA5_C", "Nhiệt độ không khí 2m ERA5 (°C)"),
        ("SoilMoisture_ERA5", "Độ ẩm đất lớp mặt ERA5 (m³/m³)")
    ]
    for tvdi_col in ["TVDI_BienRieng", "TVDI_BienGop"]:
        for vcol, vlabel in vars_check:
            sub = df_m[[tvdi_col, vcol]].dropna()
            r_val, p_val = stats.pearsonr(sub[tvdi_col], sub[vcol])
            corr_rows.append({
                "Chi_So_TVDI": tvdi_col, "Bien_Khi_Tuong": vlabel, "So_Mau_n": len(sub),
                "Pearson_r": round(r_val, 4), "R2": round(r_val**2, 4), "p_value": f"{p_val:.4e}"
            })
    df_corr = pd.DataFrame(corr_rows)
    df_corr.to_csv(os.path.join(OUT_DIR, "Bang_7_KetQua_TuongQuan_KhiTuong.csv"), index=False, encoding="utf-8-sig")
    print(df_corr.to_string(index=False))

    # Vẽ 2 biểu đồ kiểm chứng khí tượng xếp dọc (chuẩn in trắng đen)
    def _ve_kiem_chung_doc(x_col, bar_label, bar_unit, bar_fmt, bar_ylim, x_scatter_label, fig_title, trend_note, box_loc, out_name):
        yr_agg = df_m.groupby("Year")[[x_col, "TVDI_BienGop"]].mean().reindex(YEARS)
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.6, 7.6), dpi=300, gridspec_kw={"height_ratios": [1.0, 1.15], "hspace": 0.38})
        x_pos = np.arange(len(YEARS))
        bar_vals = yr_agg[x_col].values
        tvdi_vals = yr_agg["TVDI_BienGop"].values
        bars = ax1.bar(x_pos, bar_vals, width=0.44, color="#D9D9D9", edgecolor="#222222", linewidth=1.1, hatch="///", zorder=2, label=bar_label)
        ax1.set_xticks(x_pos); ax1.set_xticklabels(["2016 (El Niño)", "2021", "2026"], fontsize=10.5, fontweight="bold")
        ax1.set_ylabel(f"{bar_label} ({bar_unit})", fontsize=10.5, fontweight="bold")
        ax1.set_ylim(*bar_ylim); ax1.grid(axis="y", linestyle=":", linewidth=0.7, color="#BDBDBD", zorder=0); ax1.spines["top"].set_visible(False)
        for rect, val in zip(bars, bar_vals):
            ax1.text(rect.get_x() + rect.get_width()/2, bar_ylim[0] + 0.08*(bar_ylim[1]-bar_ylim[0]),
                     bar_fmt.format(val).replace(".", ","), ha="center", va="bottom", fontsize=10, fontweight="bold",
                     bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#666666", linewidth=0.7), zorder=5)
        ax1_r = ax1.twinx()
        ax1_r.plot(x_pos, tvdi_vals, color="#111111", lw=2.3, marker="o", ms=8.5, mfc="#111111", mec="white", mew=1.3, zorder=6, label="TVDI trung bình (biên gộp)")
        ax1_r.set_ylabel("Chỉ số TVDI trung bình", fontsize=10.5, fontweight="bold"); ax1_r.set_ylim(0.0, 0.96); ax1_r.spines["top"].set_visible(False)
        for xp, tv in zip(x_pos, tvdi_vals):
            ax1_r.annotate(f"TVDI = {tv:.2f}".replace(".", ","), xy=(xp, tv), xytext=(0, 11), textcoords="offset points",
                           ha="center", va="bottom", fontsize=9.8, fontweight="bold",
                           bbox=dict(boxstyle="round,pad=0.2", facecolor="#F7F7F7", edgecolor="#222222", linewidth=0.75), zorder=7)
        ax1.set_title("(a) Đối chiếu trung bình 10 trạm qua 3 năm", fontsize=11.2, fontweight="bold", pad=10)
        h1, l1 = ax1.get_legend_handles_labels(); h2, l2 = ax1_r.get_legend_handles_labels()
        ax1.legend(h1+h2, l1+l2, loc="upper right", fontsize=9, frameon=True, facecolor="white", edgecolor="#999999")

        x_all, y_all = df_m[x_col].values, df_m["TVDI_BienGop"].values
        sl, ic, rv, pv, _ = stats.linregress(x_all, y_all)
        xl = np.linspace(x_all.min()-0.04*(x_all.max()-x_all.min()), x_all.max()+0.04*(x_all.max()-x_all.min()), 120)
        yl = sl*xl + ic
        n = len(x_all); t_c = stats.t.ppf(0.975, df=n-2)
        s_e = np.sqrt(np.sum((y_all - (sl*x_all + ic))**2)/(n-2))
        ci = t_c * s_e * np.sqrt(1.0/n + (xl - np.mean(x_all))**2 / np.sum((x_all - np.mean(x_all))**2))
        ax2.fill_between(xl, np.clip(yl-ci, 0, 1.05), np.clip(yl+ci, 0, 1.05), color="#E0E0E0", alpha=0.75, zorder=1, label="Khoảng tin cậy 95%")
        ax2.plot(xl, yl, color="#111111", ls="--", lw=2.0, zorder=2, label="Đường xu thế tuyến tính")
        ax2.axhline(0.6, color="#888888", ls=":", lw=1.0, zorder=1)
        ax2.text(xl[-1], 0.615, "Ngưỡng hạn vừa–nặng (0,6)", ha="right", va="bottom", fontsize=8.8, style="italic", color="#555555")
        mk = {2016: dict(marker="o", s=68, fc="#1A1A1A", ec="#1A1A1A", label="Năm 2016 (El Niño)"),
              2021: dict(marker="s", s=64, fc="white", ec="#1A1A1A", lw=1.35, label="Năm 2021"),
              2026: dict(marker="^", s=74, fc="#8C8C8C", ec="#1A1A1A", label="Năm 2026")}
        for yr in YEARS:
            sub = df_m[df_m["Year"] == yr]
            ax2.scatter(sub[x_col], sub["TVDI_BienGop"], zorder=4, **mk[yr])
        ax2.set_xlabel(x_scatter_label, fontsize=10.5, fontweight="bold"); ax2.set_ylabel("Chỉ số khô hạn TVDI (Biên gộp chung)", fontsize=10.5, fontweight="bold")
        ax2.set_ylim(0.0, 1.08); ax2.grid(True, ls=":", lw=0.6, color="#CCCCCC", zorder=0); ax2.spines["top"].set_visible(False); ax2.spines["right"].set_visible(False)
        ax2.set_title("(b) Tương quan tuyến tính toàn bộ 30 mẫu (n = 30)", fontsize=11.2, fontweight="bold", pad=10)
        p_str = "p < 0,0001" if pv < 0.0001 else f"p = {pv:.4f}".replace(".", ",")
        st_txt = f"r = {rv:+.3f} | R² = {rv**2:.3f}\n{p_str} (n = 30)\n{trend_note}".replace(".", ",")
        bx, by, ha, leg_loc = (0.98, 0.96, "right", "lower left") if box_loc == "upper_right" else (0.02, 0.96, "left", "lower right")
        ax2.text(bx, by, st_txt, transform=ax2.transAxes, ha=ha, va="top", fontsize=9.3, fontweight="bold",
                 bbox=dict(boxstyle="round,pad=0.38", facecolor="white", edgecolor="#333333", lw=0.9), zorder=6)
        ax2.legend(loc=leg_loc, fontsize=8.8, frameon=True, facecolor="white", edgecolor="#999999")
        fig.suptitle(fig_title, fontsize=12.2, fontweight="bold", y=0.98)
        fig.subplots_adjust(top=0.91, bottom=0.08, left=0.11, right=0.89, hspace=0.38)
        plt.savefig(os.path.join(OUT_DIR, out_name), dpi=300, bbox_inches="tight", facecolor="white")
        plt.close(fig)

    _ve_kiem_chung_doc("Rain_CHIRPS_mm", "Lượng mưa CHIRPS", "mm", "{:.1f} mm", (0, 130),
                       "Tổng lượng mưa mùa khô 01/02 – 30/04 (mm)",
                       "Kiểm chứng chỉ số TVDI với lượng mưa mùa khô CHIRPS (10 trạm × 3 năm)",
                       "→ Mưa ít hơn, TVDI cao hơn", "upper_right", "BieuDo_KiemChung_1_LuongMua_CHIRPS.png")
    _ve_kiem_chung_doc("Temp2m_ERA5_C", "Nhiệt độ 2m ERA5", "°C", "{:.1f} °C", (22.5, 27.2),
                       "Nhiệt độ không khí trung bình mùa khô ở độ cao 2m (°C)",
                       "Kiểm chứng chỉ số TVDI với nhiệt độ không khí 2m ERA5-Land (10 trạm × 3 năm)",
                       "→ Nóng hơn, TVDI cao hơn", "upper_left", "BieuDo_KiemChung_2_NhietDo_ERA5.png")
else:
    print("   [Cảnh báo] Không tìm thấy file KhiTuong_DoiChieu_CHIRPS_ERA5_DakLak.csv — bỏ qua Bước 4.")

# =============================================================================
# BƯỚC 5: PHÂN TÍCH ẢNH HƯỞNG ĐỊA HÌNH (DEM) & LỚP PHỦ ĐẤT (ESA) (LỖI #05)
# =============================================================================
dem_path = os.path.join(DATA_DIR, "DEM_SRTM_DakLak_UTM49N.tif")
lc_path  = os.path.join(DATA_DIR, "LandCover_ESA_DakLak_UTM49N.tif")

if os.path.exists(dem_path) and os.path.exists(lc_path):
    print("\n5. Đang phân tích ảnh hưởng địa hình (SRTM DEM) và lớp phủ (ESA WorldCover)...")
    dem_arr, _ = load_and_align(dem_path, ref_profile, Resampling.bilinear)
    lc_arr, _  = load_and_align(lc_path,  ref_profile, Resampling.nearest)

    # A. Theo đai độ cao
    dem_bands = [
        ("< 200 m (Đồng bằng/Thung lũng)", (dem_arr < 200)),
        ("200 - 500 m (Đồi thấp/Cao nguyên thấp)", (dem_arr >= 200) & (dem_arr < 500)),
        ("500 - 800 m (Cao nguyên trung bình)", (dem_arr >= 500) & (dem_arr < 800)),
        ("> 800 m (Núi cao)", (dem_arr >= 800))
    ]
    dem_rows = []
    for y in YEARS:
        t = tvdi_dict[y]
        for label, cond in dem_bands:
            m = common_mask & cond & (~np.isnan(t))
            if np.sum(m) > 0:
                dem_rows.append({
                    "Nam": y, "Dai_Do_Cao": label,
                    "DienTich_km2": round(int(np.sum(m)) * PIXEL_AREA_KM2, 2),
                    "TVDI_TB": round(float(np.mean(t[m])), 4),
                    "TVDI_SD": round(float(np.std(t[m])), 4),
                    "TyLe_HanNang_%": round(float(np.mean(t[m] >= 0.8) * 100), 2)
                })
    df_dem = pd.DataFrame(dem_rows)
    df_dem.to_csv(os.path.join(OUT_DIR, "Bang_8_TVDI_Theo_DiaHinh_DEM.csv"), index=False, encoding="utf-8-sig")

    # B. Theo lớp phủ ESA WorldCover
    lc_classes = {
        10: "Rừng / Cây thân gỗ (Tree cover)",
        20: "Cây bụi (Shrubland)",
        30: "Đồng cỏ / Rừng thưa (Grassland)",
        40: "Đất nông nghiệp (Cropland)",
        50: "Đô thị / Xây dựng (Built-up)",
        60: "Đất trống (Bare / sparse)"
    }
    lc_rows = []
    for y in YEARS:
        t = tvdi_dict[y]
        for code, label in lc_classes.items():
            m = common_mask & (np.round(lc_arr) == code) & (~np.isnan(t))
            if np.sum(m) > 100:
                lc_rows.append({
                    "Nam": y, "Ma_ESA": code, "Lop_Phu": label,
                    "DienTich_km2": round(int(np.sum(m)) * PIXEL_AREA_KM2, 2),
                    "TVDI_TB": round(float(np.mean(t[m])), 4),
                    "TVDI_SD": round(float(np.std(t[m])), 4),
                    "TyLe_HanNang_%": round(float(np.mean(t[m] >= 0.8) * 100), 2)
                })
    df_lc = pd.DataFrame(lc_rows)
    df_lc.to_csv(os.path.join(OUT_DIR, "Bang_9_TVDI_Theo_LopPhu_ESA.csv"), index=False, encoding="utf-8-sig")
    print("   Đã xuất Bang_8_TVDI_Theo_DiaHinh_DEM.csv và Bang_9_TVDI_Theo_LopPhu_ESA.csv.")

    # Vẽ biểu đồ DEM & Lớp phủ ESA (chuẩn in trắng đen)
    dem_order = [b[0] for b in dem_bands]
    dem_short = ["< 200 m\nĐồng bằng / Thung lũng", "200 – 500 m\nĐồi / Cao nguyên thấp", "500 – 800 m\nCao nguyên trung bình", "> 800 m\nVùng núi cao"]
    piv_d = df_dem.pivot(index="Dai_Do_Cao", columns="Nam", values="TVDI_TB").reindex(dem_order)
    fig, ax = plt.subplots(figsize=(9.2, 4.9), dpi=300)
    xd = np.arange(len(dem_order))
    for idx, y in enumerate(YEARS):
        bs = ax.bar(xd + (idx - 1)*w, piv_d[y].values, width=w, label=f"Năm {y}", lw=0.9, **bar_styles[idx])
        for b in bs:
            ax.text(b.get_x()+b.get_width()/2, b.get_height()+0.015, f"{b.get_height():.2f}", ha="center", va="bottom", fontsize=9)
    ax.set_xticks(xd); ax.set_xticklabels(dem_short, fontsize=9.8); ax.set_ylabel("Giá trị TVDI trung bình", fontsize=11)
    ax.set_ylim(0, 0.85); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.legend(loc="upper right", ncol=3, frameon=False, fontsize=10); ax.grid(axis="y", ls=":", alpha=0.45)
    plt.tight_layout(); plt.savefig(os.path.join(OUT_DIR, "BieuDo_DiaHinh_DEM_TVDI.png"), dpi=300); plt.close()

    lc_keep = [10, 60, 50, 40, 30]
    lc_names = {30: "Đồng cỏ / Rừng thưa", 40: "Đất nông nghiệp", 50: "Đô thị / Xây dựng", 60: "Đất trống", 10: "Rừng / Cây thân gỗ"}
    piv_l = df_lc[df_lc["Ma_ESA"].isin(lc_keep)].pivot(index="Ma_ESA", columns="Nam", values="TVDI_TB").reindex(lc_keep)
    areas_l = df_lc[df_lc["Nam"] == 2016].set_index("Ma_ESA")["DienTich_km2"].to_dict()
    yl_labels = [f"{lc_names[c]}\n({areas_l.get(c, 0):,.0f} km²)".replace(",", ".") for c in lc_keep]
    fig, ax = plt.subplots(figsize=(9.2, 5.2), dpi=300)
    yl = np.arange(len(lc_keep))
    for idx, y in enumerate(YEARS):
        bs = ax.barh(yl + (1 - idx)*w, piv_l[y].values, height=w, label=f"Năm {y}", lw=0.9, **bar_styles[idx])
        for b in bs:
            ax.text(b.get_width()+0.012, b.get_y()+b.get_height()/2, f"{b.get_width():.2f}", va="center", ha="left", fontsize=9)
    ax.set_yticks(yl); ax.set_yticklabels(yl_labels, fontsize=9.8); ax.set_xlabel("Giá trị TVDI trung bình trên vùng hợp lệ chung", fontsize=11)
    ax.set_xlim(0, 0.78); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.legend(loc="lower right", frameon=False, fontsize=10); ax.grid(axis="x", ls=":", alpha=0.45)
    plt.tight_layout(); plt.savefig(os.path.join(OUT_DIR, "BieuDo_LopPhu_ESA_TVDI.png"), dpi=300); plt.close()
else:
    print("   [Cảnh báo] Chưa có DEM hoặc LandCover — bỏ qua Bước 5.")

print(f"\n>>> HOÀN TẤT CHẶNG 2! Toàn bộ bảng CSV, đồ thị PNG và GeoTIFF đã lưu tại:\n    {OUT_DIR}")