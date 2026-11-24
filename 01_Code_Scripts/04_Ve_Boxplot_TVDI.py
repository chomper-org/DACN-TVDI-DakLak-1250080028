# =============================================================================
# VẼ BIỂU ĐỒ HỘP (BOXPLOT) TVDI 3 NĂM — PHONG CÁCH TỐI GIẢN CHUẨN BÀI BÁO
# (Mô phỏng đúng bố cục tối giản của Hình 6 trong bài báo NSVI.1909)
# =============================================================================

import os
import warnings
import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import reproject, Resampling
from scipy import stats
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
plt.rcParams["font.sans-serif"] = ["Times New Roman", "DejaVu Serif", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

DATA_DIR = r"D:\DACN_TruongThanhDanh_1250080028\02_DuLieu_GeoTIFF\V2_UTM49N"
OUT_DIR  = os.path.join(DATA_DIR, "KetQua_V2")
os.makedirs(OUT_DIR, exist_ok=True)
YEARS = [2016, 2021, 2026]

def load_aligned(path, ref_profile=None):
    with rasterio.open(path) as src:
        arr = src.read(1).astype(np.float32)
        if src.nodata is not None:
            arr[arr == src.nodata] = np.nan
        if ref_profile is None:
            return arr, src.profile.copy()
        dst = np.full((ref_profile["height"], ref_profile["width"]), np.nan, dtype=np.float32)
        reproject(
            source=arr, destination=dst,
            src_transform=src.transform, src_crs=src.crs, src_nodata=np.nan,
            dst_transform=ref_profile["transform"], dst_crs=ref_profile["crs"], dst_nodata=np.nan,
            resampling=Resampling.bilinear,
        )
        return dst, ref_profile

def get_tvdi_year(year, ref_profile=None):
    cand_paths = [
        os.path.join(OUT_DIR, f"TVDI_{year}_UTM49N.tif"),
        os.path.join(r"D:\DACN_TruongThanhDanh_1250080028\03_BanDo_QGIS", f"TVDI_{year}_UTM49N.tif"),
    ]
    for p in cand_paths:
        if os.path.exists(p):
            return load_aligned(p, ref_profile)

    ndvi, prof = load_aligned(os.path.join(DATA_DIR, f"NDVI_{year}_DakLak_UTM49N.tif"), ref_profile)
    lst, _     = load_aligned(os.path.join(DATA_DIR, f"LST_{year}_DakLak_UTM49N.tif"), prof)
    valid = (~np.isnan(ndvi)) & (~np.isnan(lst)) & (ndvi >= 0.10) & (ndvi <= 0.85) & (lst > 0) & (lst < 60)
    nv, lv = ndvi[valid], lst[valid]
    bins = np.arange(0.10, 0.855, 0.01)
    dx, dy, wx, wy = [], [], [], []
    for i in range(len(bins) - 1):
        m = (nv >= bins[i]) & (nv < bins[i + 1])
        if np.sum(m) >= 50:
            mid = (bins[i] + bins[i + 1]) / 2.0
            dx.append(mid); dy.append(np.percentile(lv[m], 95))
            wx.append(mid); wy.append(np.percentile(lv[m], 5))
    rd = stats.linregress(dx, dy)
    rw = stats.linregress(wx, wy)
    tvdi = np.full(ndvi.shape, np.nan, dtype=np.float32)
    denom = (rd.slope * nv + rd.intercept) - (rw.slope * nv + rw.intercept)
    tvdi[valid] = np.clip((lv - (rw.slope * nv + rw.intercept)) / np.maximum(denom, 0.5), 0.0, 1.0)
    return tvdi, prof

print("Đang nạp dữ liệu TVDI 3 năm trên vùng hợp lệ chung...")
tvdi_map, ref_prof = {}, None
for y in YEARS:
    tvdi_map[y], ref_prof = get_tvdi_year(y, ref_prof)

common_mask = (~np.isnan(tvdi_map[2016])) & (~np.isnan(tvdi_map[2021])) & (~np.isnan(tvdi_map[2026]))

data_list, stat_rows = [], []
for y in YEARS:
    vals = tvdi_map[y][common_mask]
    data_list.append(vals)
    q1, med, q3 = np.percentile(vals, [25, 50, 75])
    stat_rows.append({
        "Nam": y, "Min": round(float(np.min(vals)), 3),
        "Q1_25%": round(float(q1), 3), "TrungVi_Median": round(float(med), 3),
        "TrungBinh_Mean": round(float(np.mean(vals)), 3),
        "Q3_75%": round(float(q3), 3), "IQR": round(float(q3 - q1), 3),
        "Max": round(float(np.max(vals)), 3)
    })

df_box = pd.DataFrame(stat_rows)
df_box.to_csv(os.path.join(OUT_DIR, "Bang_10_ThongKe_Boxplot_TVDI.csv"), index=False, encoding="utf-8-sig")
print("\nBảng thống kê phân bố TVDI (dùng viết nhận xét dưới hình):")
print(df_box.to_string(index=False))

# Vẽ Boxplot tối giản đúng chuẩn Hình 6 (NSVI.1909): hộp xám đồng nhất, sạch chữ thừa
fig, ax = plt.subplots(figsize=(6.4, 4.4), dpi=300)
bp = ax.boxplot(
    data_list,
    labels=[f"Năm {y}" for y in YEARS],
    patch_artist=True,
    widths=0.34,
    showfliers=False,
    medianprops=dict(color="#111827", linewidth=1.6),
    whiskerprops=dict(color="#4B5563", linewidth=1.1),
    capprops=dict(color="#4B5563", linewidth=1.1),
    boxprops=dict(facecolor="#A6A6A6", edgecolor="#374151", linewidth=1.1)
)

ax.set_xlabel("Thời gian", fontsize=11, color="#111827", labelpad=8)
ax.set_ylabel("Chỉ số khô hạn TVDI", fontsize=11, color="#111827", labelpad=8)
ax.set_ylim(0.0, 1.05)
ax.set_yticks(np.arange(0.0, 1.01, 0.1))
ax.tick_params(labelsize=10, colors="#111827")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.spines["bottom"].set_color("#9CA3AF")
ax.tick_params(axis="y", length=0)
ax.grid(axis="y", linestyle="-", linewidth=0.6, color="#E5E7EB")
ax.set_axisbelow(True)

out_png = os.path.join(OUT_DIR, "BieuDo_Boxplot_TVDI_3Nam.png")
plt.tight_layout()
plt.savefig(out_png, dpi=300, facecolor="white")
plt.close()
print(f"\n>>> Đã lưu biểu đồ Boxplot tối giản tại: {out_png}")