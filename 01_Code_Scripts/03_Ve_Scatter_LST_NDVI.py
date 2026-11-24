# =============================================================================
# VẼ LẠI 3 ĐỒ THỊ TÁN XẠ LST–NDVI — TỐI ƯU CẢ BẢN PDF MÀU LẪN IN TRẮNG ĐEN
# - Đám mây điểm xám dịu (không bị bệt khi photo trắng đen)
# - Biên khô (Dry Edge): đường LIỀN NÉT đậm (đỏ sẫm)
# - Biên ẩm (Wet Edge): đường NÉT ĐỨT (xanh đậm) -> photo trắng đen phân biệt ngay
# - Đồng bộ phông Times New Roman với Boxplot và bài báo mẫu NSVI.1909
# =============================================================================

import os
import warnings
import numpy as np
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

rng = np.random.default_rng(42)

for y in YEARS:
    print(f"Đang vẽ đồ thị LST–NDVI (chuẩn in trắng đen) cho năm {y}...")
    ndvi, prof = load_aligned(os.path.join(DATA_DIR, f"NDVI_{y}_DakLak_UTM49N.tif"))
    lst, _     = load_aligned(os.path.join(DATA_DIR, f"LST_{y}_DakLak_UTM49N.tif"), prof)

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
    dx, dy, wx, wy = np.array(dx), np.array(dy), np.array(wx), np.array(wy)

    rd = stats.linregress(dx, dy)
    rw = stats.linregress(wx, wy)
    a, b, dry_r2 = rd.slope, rd.intercept, rd.rvalue ** 2
    c, d, wet_r2 = rw.slope, rw.intercept, rw.rvalue ** 2
    dry_rmse = float(np.sqrt(np.mean((dy - (a * dx + b)) ** 2)))
    wet_rmse = float(np.sqrt(np.mean((wy - (c * wx + d)) ** 2)))

    idx = rng.choice(len(nv), size=min(35000, len(nv)), replace=False)
    ns, ls = nv[idx], lv[idx]

    fig, ax = plt.subplots(figsize=(6.8, 4.6), dpi=300)

    # Đám mây điểm xám nhạt sạch sẽ (in trắng đen cực nét, không bị bệt)
    ax.scatter(ns, ls, c="#9CA3AF", s=1.4, alpha=0.28, rasterized=True, zorder=2)

    # Biên khô = nét liền đậm | Biên ẩm = nét đứt (phân biệt ngay cả khi photo trắng đen)
    x_line = np.linspace(0.10, 0.85, 100)
    ax.plot(x_line, a * x_line + b, color="#991B1B", linestyle="-",  lw=2.1, zorder=4, label="Biên khô — Dry Edge (nét liền)")
    ax.plot(x_line, c * x_line + d, color="#1E3A8A", linestyle="--", lw=2.1, zorder=4, label="Biên ẩm — Wet Edge (nét đứt)")

    # Phương trình + R² + RMSE góc trên phải
    eq_dry = f"Dry Edge: y = {a:.2f}x + {b:.2f}  (R² = {dry_r2:.3f}, RMSE = {dry_rmse:.2f}°C)"
    eq_wet = f"Wet Edge: y = {c:.2f}x + {d:.2f}  (R² = {wet_r2:.3f}, RMSE = {wet_rmse:.2f}°C)"
    ax.text(0.98, 0.96, eq_dry, transform=ax.transAxes, ha="right", va="top", fontsize=9.3, color="#111827")
    ax.text(0.98, 0.89, eq_wet, transform=ax.transAxes, ha="right", va="top", fontsize=9.3, color="#111827")

    ax.set_xlabel("Chỉ số thực vật NDVI", fontsize=11, color="#111827", labelpad=6)
    ax.set_ylabel("Nhiệt độ bề mặt đất LST (°C)", fontsize=11, color="#111827", labelpad=6)
    ax.set_xlim(0.08, 0.88)
    ax.set_ylim(16, 63)
    ax.tick_params(labelsize=10, colors="#111827")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#6B7280")
    ax.spines["bottom"].set_color("#6B7280")
    ax.grid(True, linestyle="-", linewidth=0.5, color="#E5E7EB")
    ax.set_axisbelow(True)
    ax.legend(loc="lower left", frameon=False, fontsize=9.3)

    out_png = os.path.join(OUT_DIR, f"Scatter_LST_NDVI_{y}_V2.png")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, facecolor="white")
    plt.close()
    print(f"   Đã lưu: {out_png}")

print("\n>>> HOÀN TẤT! 3 ảnh Scatter chuẩn cả bản màu PDF lẫn in trắng đen.")