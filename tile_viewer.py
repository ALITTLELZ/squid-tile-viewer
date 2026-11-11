# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = [
#     "napari==0.10.0",
#     "PyQt6>=6.7",
#     "numpy",
#     "scipy",
#     "tifffile",
#     "imagecodecs",
#     "pandas",
#     "pyyaml",
#     "pillow",
# ]
# ///
"""
Tile Viewer — stitch & browse tiled fluorescence acquisitions (Squid / Cephla format)
显微镜 tile 图像浏览器

Features
  * Reads coordinates.csv + ome_tiff/R0_*.ome.tiff and stitches all tiles into one
    mosaic using image registration, with flat-field (vignetting) correction
  * napari viewer: per-channel color / contrast / gamma / on-off, z slider or
    maximum-intensity projection (MIP)
  * Export: screenshot of current view; full-resolution composite PNG of any
    rectangle (free or fixed size, optional scale bar); 16-bit multichannel TIFF
    that opens in Fiji as a colored composite

Usage
  uv run tile_viewer.py [ACQ_FOLDER] [--lang en|zh]     # open viewer (folder picker if omitted)
  uv run tile_viewer.py ACQ_FOLDER --overview            # stitch + export overview PNG only (no GUI)
"""
from __future__ import annotations

import csv
import datetime as _dt
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import tifffile

CACHE_VERSION = 2
DEFAULT_COLORS = ["#FF00FF", "#00FF00", "#00BFFF", "#FFFFFF", "#FFD000", "#FF4040"]
PREFS_PATH = Path.home() / ".tile_viewer.json"


# ----------------------------------------------------------------------------
# Language / 语言
# ----------------------------------------------------------------------------
STRINGS = {
    # errors
    "no_ome_dir": ("找不到 ome_tiff 文件夹: {p}", "Cannot find the ome_tiff folder: {p}"),
    "no_coords": ("找不到 coordinates.csv（这次采集可能没有正常结束）",
                  "Cannot find coordinates.csv (the acquisition may have been aborted)"),
    "no_tiles": ("ome_tiff 里没有可用的 tile", "No usable tiles in ome_tiff"),
    "bad_folder_title": ("无法打开这个文件夹", "Cannot open this folder"),
    "bad_folder_body": ("{err}\n\n请选择一次采集的文件夹（里面有 ome_tiff 和 coordinates.csv）。",
                        "{err}\n\nPlease choose an acquisition folder (one that contains ome_tiff and coordinates.csv)."),
    "pick_folder": ("选择一次采集的文件夹（里面有 ome_tiff）", "Choose an acquisition folder (contains ome_tiff)"),
    # stitching log
    "log_first": ("首次打开：正在拼接 {n} 个 tile（{c} 通道 × {z} 层），结果会缓存，下次秒开",
                  "First time opening: stitching {n} tiles ({c} channels x {z} z-planes). "
                  "The result is cached, so next time it opens instantly."),
    "log_pairs": ("  配准 {n} 对相邻 tile ...", "  Registering {n} neighbouring tile pairs ..."),
    "log_good": ("  可靠配对 {g}/{n}", "  Reliable pairs: {g}/{n}"),
    "log_size": ("  拼接后尺寸 {w} × {h} 像素 ({wmm:.2f} × {hmm:.2f} mm)",
                 "  Mosaic size {w} x {h} px ({wmm:.2f} x {hmm:.2f} mm)"),
    "log_flat": ("  平场校正 通道{c}：边缘/中心亮度比 {a:.2f}–{b:.2f}",
                 "  Flat-field channel {c}: relative illumination {a:.2f}-{b:.2f}"),
    "log_z": ("  第 {i}/{n} 层完成", "  z-plane {i}/{n} done"),
    "log_done": ("拼接完成，用时 {t:.0f} 秒；缓存在 {p}", "Stitching finished in {t:.0f} s; cache: {p}"),
    "overview_saved": ("总览图:", "Overview image:"),
    # panel
    "info": ("{n} tiles · {c} 通道 · {z} 层 · {px:g} µm/px · Δz {dz:g} µm<br>拼接尺寸 {w} × {h} px",
             "{n} tiles · {c} channels · {z} z-planes · {px:g} µm/px · Δz {dz:g} µm<br>Mosaic {w} x {h} px"),
    "grp_display": ("显示", "Display"),
    "mip": ("最大投影 (MIP，所有 Z 层合并)", "Max intensity projection (MIP, all z-planes)"),
    "auto": ("自动对比度（按当前画面）", "Auto contrast (current view)"),
    "reset": ("重置对比度（按整张图）", "Reset contrast (whole image)"),
    "grp_export": ("导出", "Export"),
    "roi_layer": ("ROI 选区", "ROI"),
    "hint_free": ("自由选区：点下面按钮后在图上拖出矩形（可画多个）",
                  "Free-size ROI: click the button below, then drag a rectangle on the image (multiple allowed)"),
    "bt_free": ("开始画矩形选区（自由大小）", "Draw rectangle ROI (free size)"),
    "grp_fixed": ("固定大小选框", "Fixed-size ROI"),
    "w": ("宽", "W"),
    "h": ("高", "H"),
    "unit_px": ("像素", "pixels"),
    "bt_fixed": ("添加固定大小选框（放在画面中央）", "Add fixed-size box (at view center)"),
    "fixed_tip": ("添加后直接拖动框内部即可移动；固定框的尺寸是锁定的，误拖边角也会自动恢复。选中框按 Delete 删除。",
                  "Drag inside the box to move it. Its size is locked: accidental resizing snaps back. "
                  "Select a box and press Delete to remove it."),
    "cb_bar": ("叠加图加比例尺", "Add scale bar to composite"),
    "bar_len": ("比例尺长度 µm（0=自动）", "Scale bar length µm (0 = auto)"),
    "cb_tif": ("同时保存 16-bit TIFF（Fiji 可直接打开，带颜色）", "Also save 16-bit TIFF (opens in Fiji as colored composite)"),
    "cb_stack": ("TIFF 包含全部 Z 层（否则只存当前层/投影）", "TIFF includes all z-planes (otherwise current plane / MIP only)"),
    "bt_roi": ("导出选区（全分辨率）", "Export ROIs (full resolution)"),
    "bt_all": ("导出整张大图（全分辨率）", "Export whole mosaic (full resolution)"),
    "bt_shot": ("导出当前画面截图", "Save screenshot of current view"),
    "bt_open": ("打开导出文件夹", "Open export folder"),
    "tip": ("调颜色/亮度：在左侧图层列表选中一个通道，用左上方 layer controls 调整 contrast limits（亮度对比度）、"
            "gamma、colormap（颜色）、opacity；点图层前面的眼睛开关通道。",
            "Color / brightness: select a channel in the layer list (left), then use the layer controls above it: "
            "contrast limits, gamma, colormap, opacity. Click the eye icon to show/hide a channel."),
    "lang_label": ("界面语言（重启生效）", "Language (applies on restart)"),
    "lang_saved": ("语言已保存，下次启动生效", "Language saved; it applies next time you start the viewer"),
    # notifications
    "n_auto": ("已按当前画面自动调整对比度", "Contrast adjusted to the current view"),
    "n_empty": ("选区为空或在图像外", "ROI is empty or outside the image"),
    "n_no_roi": ("还没有选区：先添加固定大小选框，或点「开始画矩形选区」在图上拖出一个框",
                 "No ROI yet: add a fixed-size box, or click 'Draw rectangle ROI' and drag on the image"),
    "n_exported": ("已导出 {n} 个文件到 {p}", "Exported {n} files to {p}"),
    "n_exported_all": ("已导出整张图到 {p}", "Exported the whole mosaic to {p}"),
    "n_shot": ("截图已保存：{f}", "Screenshot saved: {f}"),
}
LANG = "zh"


def tr(key, **kw):
    s = STRINGS[key][0 if LANG == "zh" else 1]
    return s.format(**kw) if kw else s


def load_prefs():
    try:
        return json.loads(PREFS_PATH.read_text())
    except Exception:
        return {}


def save_prefs(prefs):
    try:
        PREFS_PATH.write_text(json.dumps(prefs))
    except Exception:
        pass


def detect_lang(argv_lang=None):
    """--lang flag > TILE_VIEWER_LANG env > saved preference > system language"""
    for v in (argv_lang, os.environ.get("TILE_VIEWER_LANG"), load_prefs().get("lang")):
        if v:
            return "zh" if v.lower().startswith("zh") else "en"
    sysl = os.environ.get("LC_ALL") or os.environ.get("LANG") or ""
    if sys.platform == "darwin":
        try:
            out = subprocess.run(["defaults", "read", "-g", "AppleLanguages"],
                                 capture_output=True, text=True, timeout=3).stdout
            m = re.search(r'"?([a-zA-Z]{2})[-_"]', out)
            if m:
                sysl = m.group(1)
        except Exception:
            pass
    elif sys.platform.startswith("win"):
        try:
            import locale
            sysl = locale.getlocale()[0] or ""
        except Exception:
            pass
    return "zh" if sysl.lower().startswith(("zh", "chinese")) else "en"


# ----------------------------------------------------------------------------
# Reading an acquisition
# ----------------------------------------------------------------------------
def _load_yaml(p: Path):
    try:
        import yaml
        return yaml.safe_load(p.read_text()) if p.exists() else None
    except Exception:
        return None


class Acquisition:
    def __init__(self, folder: str | Path):
        self.folder = Path(folder).expanduser().resolve()
        tdir = self.folder / "ome_tiff"
        if not tdir.is_dir():
            raise FileNotFoundError(tr("no_ome_dir", p=tdir))
        coord = self.folder / "coordinates.csv"
        if not coord.exists():
            coord = self.folder / "0" / "coordinates.csv"
        if not coord.exists():
            raise FileNotFoundError(tr("no_coords"))

        with open(coord) as f:
            rows = list(csv.DictReader(f))
        if not rows:
            raise FileNotFoundError(tr("no_coords"))
        region = rows[0].get("region", "R0")  # first region only
        rows = [r for r in rows if r.get("region", region) == region]
        self.region = region

        def num(r, key_start):
            for k, v in r.items():
                if k.strip().lower().startswith(key_start):
                    return float(v)
            raise KeyError(key_start)

        self.stage_xy_mm = np.array([[num(r, "x"), num(r, "y")] for r in rows])

        files = {}
        for p in tdir.glob(f"{region}_*.ome.tif*"):
            m = re.match(rf"{re.escape(region)}_(\d+)\.ome\.tiff?$", p.name)
            if m:
                files[int(m.group(1))] = p
        self.tiles = {i: files[i] for i in sorted(files) if i < len(rows)}
        if not self.tiles:
            raise FileNotFoundError(tr("no_tiles"))

        first = next(iter(self.tiles.values()))
        with tifffile.TiffFile(first) as t:
            s = t.series[0]
            self.axes = s.axes
            self.tile_shape = s.shape
            self.dtype = s.dtype
            ome = t.ome_metadata or ""
        ax = self.axes
        self.nz = self.tile_shape[ax.index("Z")] if "Z" in ax else 1
        self.nc = self.tile_shape[ax.index("C")] if "C" in ax else 1
        self.th = self.tile_shape[ax.index("Y")]
        self.tw = self.tile_shape[ax.index("X")]

        self.channel_names = re.findall(r'<Channel [^>]*Name="([^"]*)"', ome)[: self.nc]
        while len(self.channel_names) < self.nc:
            self.channel_names.append(f"Channel {len(self.channel_names)}")

        acq = _load_yaml(self.folder / "acquisition.yaml") or {}
        try:
            px = float(acq["objective"]["pixel_size_um"])
        except Exception:
            m = re.search(r'PhysicalSizeX="([\d.eE+-]+)"', ome)
            px = float(m.group(1)) if m else 1.0
        self.px_um = px
        try:
            dz = float(acq["z_stack"]["delta_z_mm"]) * 1000
        except Exception:
            m = re.search(r'PhysicalSizeZ="([\d.eE+-]+)"', ome)
            dz = float(m.group(1)) if m else 1.0
        self.dz_um = dz

        colors = {}
        chy = _load_yaml(self.folder / "acquisition_channels.yaml") or {}
        for ch in (chy.get("channels") or acq.get("channels") or []):
            if isinstance(ch, dict) and ch.get("name") and ch.get("display_color"):
                colors[ch["name"]] = ch["display_color"]
        self.channel_colors = []
        for i, n in enumerate(self.channel_names):
            c = colors.get(n)
            if c is None:
                for k, v in colors.items():
                    if k in n or n in k:
                        c = v
            self.channel_colors.append(c or DEFAULT_COLORS[i % len(DEFAULT_COLORS)])

        self.cache_dir = self.folder / "_viewer_cache"
        self._mm = {}

    def tile(self, i) -> np.ndarray:
        """Return tile i as a (Z, C, Y, X) array (memory-mapped when possible)."""
        if i in self._mm:
            return self._mm[i]
        p = self.tiles[i]
        try:
            a = tifffile.memmap(p)
        except Exception:
            a = tifffile.imread(p)
        if a.shape != self.tile_shape:
            a = a.reshape(self.tile_shape)
        order = [self.axes.index(k) for k in "ZCYX" if k in self.axes]
        a = np.transpose(a, order) if order != sorted(order) else a
        if "Z" not in self.axes:
            a = a[None]
        if "C" not in self.axes:
            a = a[:, None]
        self._mm[i] = a
        return a

    def key(self):
        h = hashlib.sha1()
        for i, p in self.tiles.items():
            st = p.stat()
            h.update(f"{i}:{p.name}:{st.st_size}:{int(st.st_mtime)}".encode())
        h.update(f"v{CACHE_VERSION}".encode())
        return h.hexdigest()[:16]


# ----------------------------------------------------------------------------
# Registration + stitching
# ----------------------------------------------------------------------------
def _reg_image(acq: Acquisition, i, ds=2):
    """Sum of normalized per-channel MIPs, downsampled; used for registration."""
    a = acq.tile(i)
    out = np.zeros((acq.th // ds, acq.tw // ds), np.float32)
    for c in range(acq.nc):
        m = a[:, c, : out.shape[0] * ds : ds, : out.shape[1] * ds : ds].max(axis=0).astype(np.float32)
        lo, hi = np.percentile(m, (1, 99.5))
        out += np.clip((m - lo) / max(hi - lo, 1e-6), 0, 1)
    return out


def _phase_corr(a, b):
    """Return d such that a(x) ≈ b(x - d). Zero-padded to avoid wrap-around."""
    a = a - a.mean()
    b = b - b.mean()
    sh = (a.shape[0] * 2, a.shape[1] * 2)
    Fa = np.fft.rfft2(a, sh)
    Fb = np.fft.rfft2(b, sh)
    R = Fa * np.conj(Fb)
    R /= np.abs(R) ** 0.7 + 1e-9  # partial whitening: more robust than pure phase correlation
    r = np.fft.irfft2(R, sh)
    py, px = np.unravel_index(np.argmax(r), r.shape)
    if py > sh[0] // 2:
        py -= sh[0]
    if px > sh[1] // 2:
        px -= sh[1]
    return py, px


def _ncc(a, b):
    a = a - a.mean()
    b = b - b.mean()
    d = np.sqrt((a * a).sum() * (b * b).sum())
    return float((a * b).sum() / d) if d > 0 else 0.0


def compute_positions(acq: Acquisition, log=print):
    ds = 2
    idx = list(acq.tiles)
    xy = acq.stage_xy_mm[idx]

    def step(v):
        d = np.diff(np.unique(np.round(v, 4)))
        d = d[d > 1e-3]
        return float(np.median(d)) if len(d) else None

    sx, sy = step(xy[:, 0]), step(xy[:, 1])
    stage_px = acq.stage_xy_mm * 1000.0 / acq.px_um

    pairs = []
    for a_i in idx:
        for b_i in idx:
            if a_i == b_i:
                continue
            dx, dy = acq.stage_xy_mm[b_i] - acq.stage_xy_mm[a_i]
            if sx and abs(dx - sx) < 0.2 * sx and abs(dy) < 0.2 * sx:
                pairs.append((a_i, b_i, "x"))
            if sy and abs(dy - sy) < 0.2 * sy and abs(dx) < 0.2 * sy:
                pairs.append((a_i, b_i, "y"))

    log(tr("log_pairs", n=len(pairs)))
    imgs = {i: _reg_image(acq, i, ds) for i in idx}
    H, W = imgs[idx[0]].shape
    S = min(W, int(1300 / ds))  # strip width used for registration (downsampled px)
    results = []
    for a_i, b_i, d in pairs:
        A, B = imgs[a_i], imgs[b_i]
        if d == "x":
            sa, sb = A[:, W - S:], B[:, :S]
        else:
            sa, sb = A[H - S:, :], B[:S, :]
        py, px = _phase_corr(sa, sb)
        if d == "x":
            off = np.array([py, (W - S) + px], float)
        else:
            off = np.array([(H - S) + py, px], float)
        oy, ox = int(off[0]), int(off[1])
        ya0, ya1 = max(0, oy), min(H, H + oy)
        xa0, xa1 = max(0, ox), min(W, W + ox)
        ok = (ya1 - ya0) > 20 and (xa1 - xa0) > 20
        score = _ncc(A[ya0:ya1, xa0:xa1], B[ya0 - oy:ya1 - oy, xa0 - ox:xa1 - ox]) if ok else 0
        nominal = (stage_px[b_i] - stage_px[a_i])[::-1] / ds
        plausible = ok and np.all(np.abs(off - nominal) < 0.25 * max(H, W))
        results.append((a_i, b_i, d, off * ds, score if plausible else 0.0))

    good = [r for r in results if r[4] > 0.3]
    log(tr("log_good", g=len(good), n=len(results)))
    med = {}
    for d in "xy":
        g = [r[3] for r in good if r[2] == d]
        if g:
            med[d] = np.median(np.array(g), axis=0)
    eqs = []
    for a_i, b_i, d, off, sc in results:
        if sc > 0.3:
            eqs.append((a_i, b_i, off, sc))
        else:  # unreliable pair (e.g. empty glass): fall back to the median offset
            fallback = med.get(d, (stage_px[b_i] - stage_px[a_i])[::-1])
            eqs.append((a_i, b_i, np.asarray(fallback), 0.05))
    # global least-squares placement
    n = len(idx)
    pos_of = {t: k for k, t in enumerate(idx)}
    rows, rhs, w = [], [], []
    for a_i, b_i, off, sc in eqs:
        r = np.zeros(n)
        r[pos_of[b_i]] = 1
        r[pos_of[a_i]] = -1
        rows.append(r)
        rhs.append(off)
        w.append(sc)
    scale_est = 1.0
    if "x" in med and sx:
        scale_est = med["x"][1] / (sx * 1000 / acq.px_um)
    for t in idx:  # weak anchors so isolated tiles still get a position
        r = np.zeros(n)
        r[pos_of[t]] = 1
        rows.append(r)
        rhs.append(stage_px[t][::-1] * scale_est)
        w.append(1e-4)
    M = np.array(rows) * np.array(w)[:, None]
    Y = np.array(rhs) * np.array(w)[:, None]
    P, *_ = np.linalg.lstsq(M, Y, rcond=None)
    P -= P.min(axis=0)
    pos = {t: [int(round(P[pos_of[t], 0])), int(round(P[pos_of[t], 1]))] for t in idx}
    info = {
        "pairs_total": len(results),
        "pairs_good": len(good),
        "median_offset_x": med.get("x", np.zeros(2)).tolist(),
        "median_offset_y": med.get("y", np.zeros(2)).tolist(),
    }
    return pos, info


def estimate_flatfield(acq: Acquisition, log=print, ds=8):
    """Estimate per-channel illumination profile from all tiles.
    Returns (flats[C,H,W] normalized to mean 1, darks[C])."""
    from scipy.ndimage import gaussian_filter, zoom
    flats, darks = [], []
    for c in range(acq.nc):
        imgs = []
        for i in acq.tiles:
            a = acq.tile(i)
            for z in range(acq.nz):
                imgs.append(a[z, c, ::ds, ::ds].astype(np.float32))
        imgs = np.stack(imgs)
        dark = float(np.percentile(imgs, 0.5))
        means = imgs.mean(axis=(1, 2))
        keep = imgs[means >= np.median(means)]  # use tissue-rich images only
        keep = (keep - dark) / np.maximum(keep.mean(axis=(1, 2), keepdims=True) - dark, 1e-3)
        f = np.median(keep, axis=0)
        f = gaussian_filter(f, sigma=max(f.shape) / 12, mode="nearest")
        f = np.clip(f / f.mean(), 0.2, 5)
        F = zoom(f, (acq.th / f.shape[0], acq.tw / f.shape[1]), order=1)[: acq.th, : acq.tw]
        flats.append(F.astype(np.float32))
        darks.append(dark)
        log(tr("log_flat", c=c, a=f.min(), b=f.max()))
    return np.stack(flats), darks


def _feather(h, w):
    y = np.minimum(np.arange(h), np.arange(h)[::-1]).astype(np.float32) + 1
    x = np.minimum(np.arange(w), np.arange(w)[::-1]).astype(np.float32) + 1
    return np.minimum(y[:, None], x[None, :])


def _downsample2(a):
    """2x2 mean over the last two axes."""
    h, w = a.shape[-2] // 2 * 2, a.shape[-1] // 2 * 2
    a = a[..., :h, :w].astype(np.float32)
    return (a[..., 0::2, 0::2] + a[..., 1::2, 0::2] + a[..., 0::2, 1::2] + a[..., 1::2, 1::2]) / 4


def build_mosaic(acq: Acquisition, log=print, force=False, flatfield=True):
    """Stitch and write the cache. Returns dict(levels=[ZCYX memmaps], mip=[CYX memmaps], meta)."""
    cd = acq.cache_dir
    meta_p = cd / "meta.json"
    key = acq.key() + ("-ff" if flatfield else "-raw")
    if not force and meta_p.exists():
        meta = json.loads(meta_p.read_text())
        if meta.get("key") == key:
            return _open_cache(cd, meta)

    cd.mkdir(exist_ok=True)
    for p in cd.glob("*.npy"):
        p.unlink()
    t0 = time.time()
    log(tr("log_first", n=len(acq.tiles), c=acq.nc, z=acq.nz))
    pos, info = compute_positions(acq, log)
    H = max(p[0] for p in pos.values()) + acq.th
    W = max(p[1] for p in pos.values()) + acq.tw
    log(tr("log_size", w=W, h=H, wmm=W * acq.px_um / 1000, hmm=H * acq.px_um / 1000))

    if flatfield:
        flats, darks = estimate_flatfield(acq, log)
    L0 = np.lib.format.open_memmap(cd / "L0.npy", "w+", np.uint16, (acq.nz, acq.nc, H, W))
    fw = _feather(acq.th, acq.tw)
    wsum = np.zeros((H, W), np.float32)
    for i, (y, x) in pos.items():
        wsum[y:y + acq.th, x:x + acq.tw] += fw
    wsum[wsum == 0] = 1
    canvas = np.zeros((H, W), np.float32)
    for z in range(acq.nz):
        for c in range(acq.nc):
            canvas[:] = 0
            for i, (y, x) in pos.items():
                t = acq.tile(i)[z, c].astype(np.float32)
                if flatfield:
                    t = (t - darks[c]) / flats[c] + darks[c]
                canvas[y:y + acq.th, x:x + acq.tw] += t * fw
            L0[z, c] = np.clip(canvas / wsum + 0.5, 0, 65535).astype(np.uint16)
        log(tr("log_z", i=z + 1, n=acq.nz))
    L0.flush()

    M0 = np.lib.format.open_memmap(cd / "M0.npy", "w+", np.uint16, (acq.nc, H, W))
    for c in range(acq.nc):
        M0[c] = L0[:, c].max(axis=0)
    M0.flush()

    # resolution pyramid for smooth zooming
    nlev = 1
    prevL, prevM = L0, M0
    while max(prevL.shape[-2:]) > 1024:
        h2, w2 = prevL.shape[-2] // 2, prevL.shape[-1] // 2
        Lk = np.lib.format.open_memmap(cd / f"L{nlev}.npy", "w+", np.uint16, (acq.nz, acq.nc, h2, w2))
        for z in range(acq.nz):
            Lk[z] = _downsample2(prevL[z]).astype(np.uint16)
        Mk = np.lib.format.open_memmap(cd / f"M{nlev}.npy", "w+", np.uint16, (acq.nc, h2, w2))
        Mk[:] = _downsample2(prevM).astype(np.uint16)
        Lk.flush()
        Mk.flush()
        prevL, prevM = Lk, Mk
        nlev += 1

    meta = dict(key=key, nlev=nlev, flatfield=flatfield, positions={str(k): v for k, v in pos.items()},
                info=info, shape=[acq.nz, acq.nc, H, W], built=_dt.datetime.now().isoformat(timespec="seconds"))
    meta_p.write_text(json.dumps(meta, indent=1))
    log(tr("log_done", t=time.time() - t0, p=cd))
    return _open_cache(cd, meta)


def _open_cache(cd, meta):
    levels = [np.load(cd / f"L{k}.npy", mmap_mode="r") for k in range(meta["nlev"])]
    mips = [np.load(cd / f"M{k}.npy", mmap_mode="r") for k in range(meta["nlev"])]
    return dict(levels=levels, mip=mips, meta=meta)


# ----------------------------------------------------------------------------
# Rendering & export (no napari needed)
# ----------------------------------------------------------------------------
def hex_to_rgb(h):
    h = h.lstrip("#")
    return np.array([int(h[i:i + 2], 16) for i in (0, 2, 4)], np.float32) / 255


def auto_limits(a, lo=0.5, hi=99.8):
    a = np.asarray(a)
    a = a[a > 0] if (a > 0).any() else a.ravel()
    if a.size > 2_000_000:
        a = a[:: a.size // 2_000_000]
    l, h = np.percentile(a, (lo, hi))
    return float(l), float(max(h, l + 1))


def composite_rgb(chans, settings):
    """chans: list of 2D arrays; settings: list of dict(visible, lo, hi, gamma, opacity, cmap).
    cmap: callable mapping (N,) in [0,1] -> (N,3), or an RGB color."""
    h, w = chans[0].shape
    out = np.zeros((h, w, 3), np.float32)
    for a, s in zip(chans, settings):
        if not s["visible"]:
            continue
        step = max(1, 4_000_000 // max(w, 1))
        for r0 in range(0, h, step):
            blk = a[r0:r0 + step].astype(np.float32)
            n = np.clip((blk - s["lo"]) / max(s["hi"] - s["lo"], 1e-9), 0, 1)
            if s.get("gamma", 1) != 1:
                n = n ** s["gamma"]
            cm = s["cmap"]
            if callable(cm):
                rgb = cm(n.ravel()).reshape(n.shape + (3,))
            else:
                rgb = n[..., None] * np.asarray(cm, np.float32)
            out[r0:r0 + step] += rgb * s.get("opacity", 1)
    return (np.clip(out, 0, 1) * 255 + 0.5).astype(np.uint8)


def _nice_bar(width_um):
    target = width_um / 5
    for v in [1, 2, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000, 2000, 5000, 10000]:
        if v >= target * 0.6:
            return v
    return 10000


def _font(size):
    from PIL import ImageFont
    for fp in ("/System/Library/Fonts/Supplemental/Arial.ttf", "/System/Library/Fonts/Helvetica.ttc",
               "/Library/Fonts/Arial.ttf", "C:/Windows/Fonts/arial.ttf",
               "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"):
        try:
            return ImageFont.truetype(fp, size)
        except Exception:
            pass
    try:
        return ImageFont.load_default(size=size)
    except TypeError:
        return ImageFont.load_default()


def add_scale_bar(rgb, px_um, bar_um=None):
    from PIL import Image, ImageDraw
    h, w = rgb.shape[:2]
    bar_um = bar_um or _nice_bar(w * px_um)
    L = int(round(bar_um / px_um))
    if L < 4 or L > w * 0.9:
        return rgb
    im = Image.fromarray(rgb)
    d = ImageDraw.Draw(im)
    th = max(3, int(min(h, w) / 120))
    m = max(10, int(min(h, w) / 30))
    fs = max(12, int(min(h, w) / 30))
    font = _font(fs)
    x1, y1 = w - m, h - m
    d.rectangle([x1 - L, y1 - th, x1, y1], fill=(255, 255, 255))
    label = f"{bar_um:g} µm" if bar_um < 1000 else f"{bar_um / 1000:g} mm"
    tb = d.textbbox((0, 0), label, font=font)
    tw = tb[2] - tb[0]
    d.text((x1 - L / 2 - tw / 2, y1 - th - (tb[3] - tb[1]) - fs * 0.4), label, fill=(255, 255, 255), font=font)
    return np.asarray(im)


def save_imagej(path, zcyx, px_um, dz_um, colors_rgb, ranges, names):
    luts = []
    for c in colors_rgb:
        ramp = np.linspace(0, 1, 256)[:, None] * np.asarray(c)[None, :]
        luts.append((np.clip(ramp, 0, 1) * 255).astype(np.uint8).T.copy())
    tifffile.imwrite(
        path, np.ascontiguousarray(zcyx), imagej=True,
        resolution=(1 / px_um, 1 / px_um),
        metadata={"axes": "ZCYX", "mode": "composite", "unit": "um", "spacing": dz_um,
                  "LUTs": luts, "Ranges": tuple(v for r in ranges for v in r),
                  "Labels": [n for _ in range(zcyx.shape[0]) for n in names]},
    )


def export_overview(acq, mosaic, out_path=None, max_px=4000):
    """Non-interactive: export a MIP overview PNG (default colors + auto contrast)."""
    levels = mosaic["mip"]
    k = 0
    while k + 1 < len(levels) and max(levels[k].shape[-2:]) > max_px:
        k += 1
    m = levels[k]
    small = mosaic["mip"][-1]
    settings = []
    for c in range(acq.nc):
        lo, hi = auto_limits(small[c])
        settings.append(dict(visible=True, lo=lo, hi=hi, gamma=1, opacity=1, cmap=hex_to_rgb(acq.channel_colors[c])))
    rgb = composite_rgb([np.asarray(m[c]) for c in range(acq.nc)], settings)
    rgb = add_scale_bar(rgb, acq.px_um * (2 ** k))
    out = Path(out_path) if out_path else acq.folder / "exports" / "overview_MIP.png"
    out.parent.mkdir(exist_ok=True)
    from PIL import Image
    Image.fromarray(rgb).save(out)
    return out


# ----------------------------------------------------------------------------
# napari GUI
# ----------------------------------------------------------------------------
def run_viewer(folder: str | None):
    import napari
    import pandas as pd
    from napari.utils.colormaps import Colormap
    from napari.utils.notifications import show_info, show_warning
    from qtpy.QtCore import Qt, QUrl
    from qtpy.QtGui import QDesktopServices
    from qtpy.QtWidgets import (QApplication, QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QFrame,
                                QGroupBox, QHBoxLayout, QLabel, QMessageBox, QPushButton, QScrollArea,
                                QSizePolicy, QSpinBox, QVBoxLayout, QWidget)

    prefs = load_prefs()
    app = QApplication.instance() or QApplication(sys.argv)  # noqa: F841

    acq = None
    while acq is None:
        if not folder:
            start = prefs.get("last_parent", str(Path.home() / "Documents"))
            folder = QFileDialog.getExistingDirectory(None, tr("pick_folder"), start)
            if not folder:
                return None
        try:
            acq = Acquisition(folder)
        except Exception as e:
            QMessageBox.warning(None, tr("bad_folder_title"), tr("bad_folder_body", err=e))
            folder = None
    prefs["last_parent"] = str(acq.folder.parent)
    save_prefs(prefs)

    mosaic = build_mosaic(acq)
    levels, mips = mosaic["levels"], mosaic["mip"]

    viewer = napari.Viewer(title=f"Tile Viewer — {acq.folder.name}")
    scale = (acq.dz_um, acq.px_um, acq.px_um)
    small_mip = mips[-1]
    ch_layers = []
    for c in range(acq.nc):
        col = acq.channel_colors[c]
        cmap = Colormap([[0, 0, 0, 1], list(hex_to_rgb(col)) + [1]], name=f"{acq.channel_names[c]} ({col})")
        lo, hi = auto_limits(small_mip[c])
        dmax = float(np.asarray(small_mip[c]).max())
        lay = viewer.add_image(
            [lv[:, c] for lv in levels], multiscale=True, name=acq.channel_names[c],
            colormap=cmap, blending="additive", scale=scale, contrast_limits=(lo, hi),
        )
        lay.contrast_limits_range = (0, max(dmax * 1.5, hi * 1.2, 10))
        lay.contrast_limits = (lo, hi)
        ch_layers.append(lay)

    roi = viewer.add_shapes(
        ndim=2, name=tr("roi_layer"), scale=scale[1:], edge_color="white", face_color=[0, 0, 0, 0], edge_width=6,
        features=pd.DataFrame({"fw": pd.Series(dtype=float), "fh": pd.Series(dtype=float),
                               "label": pd.Series(dtype=str)}),
        text={"string": "{label}", "color": "white", "size": 11, "anchor": "upper_left"},
    )
    roi.feature_defaults = {"fw": 0.0, "fh": 0.0, "label": ""}
    roi.mode = "pan_zoom"
    viewer.scale_bar.visible = True
    for l in ch_layers + [roi]:
        try:
            l.units = ("um",) * l.ndim  # napari >= 0.6
        except Exception:
            pass
    try:
        viewer.scale_bar.unit = "um"  # older napari
    except Exception:
        pass
    viewer.dims.axis_labels = ("z", "y", "x")
    viewer.dims.set_current_step(0, acq.nz // 2)
    viewer.layers.selection.active = ch_layers[0]

    state = {"mip": False, "z": acq.nz // 2}
    out_dir = acq.folder / "exports"

    # ---------------- panel ----------------
    w = QWidget()
    lay_v = QVBoxLayout(w)

    info = QLabel(f"<b>{acq.folder.name}</b><br>" + tr(
        "info", n=len(acq.tiles), c=acq.nc, z=acq.nz, px=acq.px_um, dz=acq.dz_um,
        w=levels[0].shape[-1], h=levels[0].shape[-2]))
    info.setWordWrap(True)
    lay_v.addWidget(info)

    g1 = QGroupBox(tr("grp_display"))
    v1 = QVBoxLayout(g1)
    cb_mip = QCheckBox(tr("mip"))
    v1.addWidget(cb_mip)
    bt_auto = QPushButton(tr("auto"))
    v1.addWidget(bt_auto)
    bt_reset = QPushButton(tr("reset"))
    v1.addWidget(bt_reset)
    lay_v.addWidget(g1)

    g2 = QGroupBox(tr("grp_export"))
    v2 = QVBoxLayout(g2)

    gf = QGroupBox(tr("grp_fixed"))
    vf = QVBoxLayout(gf)
    hf = QHBoxLayout()
    sp_w = QDoubleSpinBox()
    sp_h = QDoubleSpinBox()
    for sp in (sp_w, sp_h):
        sp.setRange(1, 1e6)
        sp.setDecimals(1)
        sp.setValue(prefs.get("fixed_w", 500.0) if sp is sp_w else prefs.get("fixed_h", 500.0))
    cb_unit = QComboBox()
    cb_unit.addItems(["µm", tr("unit_px")])
    cb_unit.setCurrentIndex(prefs.get("fixed_unit", 0))
    hf.addWidget(QLabel(tr("w")))
    hf.addWidget(sp_w)
    hf.addWidget(QLabel(tr("h")))
    hf.addWidget(sp_h)
    hf.addWidget(cb_unit)
    vf.addLayout(hf)
    bt_fixed = QPushButton(tr("bt_fixed"))
    vf.addWidget(bt_fixed)
    lb_fixed = QLabel(tr("fixed_tip"))
    lb_fixed.setWordWrap(True)
    lb_fixed.setStyleSheet("color: gray")
    vf.addWidget(lb_fixed)
    v2.addWidget(gf)

    hint = QLabel(tr("hint_free"))
    hint.setWordWrap(True)
    v2.addWidget(hint)
    bt_rect = QPushButton(tr("bt_free"))
    v2.addWidget(bt_rect)

    cb_bar = QCheckBox(tr("cb_bar"))
    cb_bar.setChecked(True)
    hb = QHBoxLayout()
    hb.addWidget(QLabel(tr("bar_len")))
    sp_bar = QSpinBox()
    sp_bar.setRange(0, 100000)
    sp_bar.setValue(0)
    hb.addWidget(sp_bar)
    v2.addWidget(cb_bar)
    v2.addLayout(hb)
    cb_tif = QCheckBox(tr("cb_tif"))
    cb_tif.setChecked(True)
    v2.addWidget(cb_tif)
    cb_stack = QCheckBox(tr("cb_stack"))
    v2.addWidget(cb_stack)
    bt_roi = QPushButton(tr("bt_roi"))
    bt_all = QPushButton(tr("bt_all"))
    bt_shot = QPushButton(tr("bt_shot"))
    bt_open = QPushButton(tr("bt_open"))
    for b in (bt_roi, bt_all, bt_shot, bt_open):
        v2.addWidget(b)
    lay_v.addWidget(g2)

    tip = QLabel(tr("tip"))
    tip.setWordWrap(True)
    tip.setStyleSheet("color: gray")
    lay_v.addWidget(tip)

    hl = QHBoxLayout()
    hl.addWidget(QLabel(tr("lang_label")))
    cb_lang = QComboBox()
    cb_lang.addItems(["中文", "English"])
    cb_lang.setCurrentIndex(0 if LANG == "zh" else 1)
    hl.addWidget(cb_lang)
    lay_v.addLayout(hl)
    lay_v.addStretch(1)
    for g in (g1, g2, gf):
        g.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Minimum)
    scroll = QScrollArea()
    scroll.setWidgetResizable(True)
    scroll.setFrameShape(QFrame.Shape.NoFrame)
    scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
    scroll.setWidget(w)
    scroll.setMinimumWidth(w.sizeHint().width() + 20)
    viewer.window.add_dock_widget(scroll, name="Tile Viewer", area="right")

    # ---------------- logic ----------------
    def on_lang(i):
        p = load_prefs()
        p["lang"] = "zh" if i == 0 else "en"
        save_prefs(p)
        show_info(STRINGS["lang_saved"][i])

    cb_lang.currentIndexChanged.connect(on_lang)

    def mip_data(c):
        # Broadcast the MIP to the same (Z, Y, X) shape as the z-stack (no extra memory),
        # so swapping layer data never changes its shape. Changing shape made napari draw
        # one frame with stale tile coordinates -> empty texture -> OpenGL error.
        return [np.broadcast_to(m[c][None], (acq.nz,) + m.shape[-2:]) for m in mips]

    def set_mip(on):
        if on == state["mip"]:
            return
        for c, l in enumerate(ch_layers):
            cl, cr = l.contrast_limits, l.contrast_limits_range
            l.data = mip_data(c) if on else [lv[:, c] for lv in levels]
            l.contrast_limits_range = cr
            l.contrast_limits = cl
        state["mip"] = on

    cb_mip.toggled.connect(set_mip)

    def visible_settings():
        out = []
        for l in ch_layers:
            cm = l.colormap
            out.append(dict(visible=l.visible, lo=l.contrast_limits[0], hi=l.contrast_limits[1],
                            gamma=l.gamma, opacity=l.opacity,
                            cmap=(lambda n, cm=cm: cm.map(n)[:, :3])))
        return out

    def plane(c, y0, y1, x0, x1, level=0):
        if state["mip"]:
            return np.asarray(mips[level][c, y0:y1, x0:x1])
        return np.asarray(levels[level][viewer.dims.current_step[0], c, y0:y1, x0:x1])

    def on_auto():
        l0 = ch_layers[0]
        k = l0.data_level
        cp = l0.corner_pixels
        y0, x0 = int(cp[0, -2]), int(cp[0, -1])
        y1, x1 = int(cp[1, -2]) + 1, int(cp[1, -1]) + 1
        for c, l in enumerate(ch_layers):
            a = plane(c, y0, y1, x0, x1, level=k)
            if a.size:
                l.contrast_limits = auto_limits(a)
        show_info(tr("n_auto"))

    def on_reset():
        for c, l in enumerate(ch_layers):
            l.contrast_limits = auto_limits(mips[-1][c])

    bt_auto.clicked.connect(on_auto)
    bt_reset.clicked.connect(on_reset)

    def on_rect():
        roi.selected_data = set()
        roi.feature_defaults = {"fw": 0.0, "fh": 0.0, "label": ""}
        viewer.layers.selection.active = roi
        roi.mode = "add_rectangle"

    bt_rect.clicked.connect(on_rect)

    H0, W0 = levels[0].shape[-2:]

    def fixed_rect(cy, cx, fh, fw):
        """fh x fw px rectangle centred on (cy, cx), shifted to stay inside the image."""
        fh, fw = min(fh, H0), min(fw, W0)
        y0 = int(round(cy - fh / 2))
        x0 = int(round(cx - fw / 2))
        y0 = min(max(y0, 0), int(H0 - fh))
        x0 = min(max(x0, 0), int(W0 - fw))
        y1, x1 = y0 + fh, x0 + fw
        return np.array([[y0, x0], [y0, x1], [y1, x1], [y1, x0]], float)

    snapping = {"on": False, "adding_fixed": False, "n": 0}

    def add_fixed():
        k = acq.px_um if cb_unit.currentIndex() == 0 else 1.0
        fw = int(round(sp_w.value() / k))
        fh = int(round(sp_h.value() / k))
        if cb_unit.currentIndex() == 0:
            label = f"{sp_w.value():g}×{sp_h.value():g} µm"
        else:
            label = f"{fw}×{fh} px"
        p = load_prefs()
        p.update(fixed_w=sp_w.value(), fixed_h=sp_h.value(), fixed_unit=cb_unit.currentIndex())
        save_prefs(p)
        cy, cx = np.asarray(viewer.camera.center)[-2:] / acq.px_um
        roi.selected_data = set()
        roi.feature_defaults = {"fw": float(fw), "fh": float(fh), "label": label}
        snapping["adding_fixed"] = True
        try:
            roi.add_rectangles(fixed_rect(cy, cx, fh, fw))
        finally:
            snapping["adding_fixed"] = False
        roi.feature_defaults = {"fw": 0.0, "fh": 0.0, "label": ""}
        snapping["n"] = len(roi.data)
        roi.refresh_text()
        viewer.layers.selection.active = roi
        roi.mode = "select"
        roi.selected_data = {len(roi.data) - 1}

    bt_fixed.clicked.connect(add_fixed)

    def snap_fixed(event=None):
        """Restore the size of fixed boxes after accidental resizing (keeps centre);
        any shape not added via 'Add fixed-size box' is treated as free-size."""
        if snapping["on"]:
            return
        n_prev, n = snapping["n"], len(roi.data)
        action = str(getattr(event, "action", "changed")).lower()
        if n > n_prev and not snapping["adding_fixed"] and "add" in action:
            f = roi.features.copy()
            f.loc[f.index[n_prev:n], "fw"] = 0.0
            f.loc[f.index[n_prev:n], "fh"] = 0.0
            f.loc[f.index[n_prev:n], "label"] = ""
            snapping["on"] = True
            try:
                roi.features = f
                roi.feature_defaults = {"fw": 0.0, "fh": 0.0, "label": ""}
                roi.refresh_text()
            finally:
                snapping["on"] = False
        if "add" in action or "remov" in action or n < n_prev:
            snapping["n"] = n
        if n == 0 or not ("chang" in action or "add" in action) or action.endswith("changing"):
            return
        feats = roi.features
        fixes = {}
        for i, s in enumerate(roi.data):
            fw = float(feats["fw"].iloc[i]) if i < len(feats) else 0
            fh = float(feats["fh"].iloc[i]) if i < len(feats) else 0
            s = np.asarray(s)[:, -2:]
            if fw > 0 and fh > 0:
                c = s.mean(axis=0)
                r = fixed_rect(c[0], c[1], int(fh), int(fw))
                if s.shape != r.shape or not np.allclose(np.sort(s, 0), np.sort(r, 0), atol=0.5):
                    fixes[i] = r
        if fixes:
            snapping["on"] = True
            try:
                for i, r in fixes.items():
                    roi._data_view.edit(i, r, new_type="rectangle")
                roi.refresh()
                roi.refresh_text()
            finally:
                snapping["on"] = False

    roi.events.data.connect(snap_fixed)

    def tag():
        z = "MIP" if state["mip"] else f"z{viewer.dims.current_step[0]:02d}"
        return z, _dt.datetime.now().strftime("%Y%m%d-%H%M%S")

    def export_region(y0, y1, x0, x1, name):
        from PIL import Image
        out_dir.mkdir(exist_ok=True)
        H, W = levels[0].shape[-2:]
        y0, x0 = max(0, y0), max(0, x0)
        y1, x1 = min(H, y1), min(W, x1)
        if y1 - y0 < 2 or x1 - x0 < 2:
            show_warning(tr("n_empty"))
            return []
        chans = [plane(c, y0, y1, x0, x1) for c in range(acq.nc)]
        rgb = composite_rgb(chans, visible_settings())
        if cb_bar.isChecked():
            rgb = add_scale_bar(rgb, acq.px_um, sp_bar.value() or None)
        z, ts = tag()
        base = out_dir / f"{name}_{z}_{ts}"
        Image.fromarray(rgb).save(f"{base}.png")
        files = [f"{base}.png"]
        if cb_tif.isChecked():
            if cb_stack.isChecked() and not state["mip"]:
                arr = np.asarray(levels[0][:, :, y0:y1, x0:x1])
            else:
                arr = np.stack(chans)[None]
            cols = [l.colormap.map(np.array([1.0]))[0, :3] for l in ch_layers]
            rng = [tuple(map(float, l.contrast_limits)) for l in ch_layers]
            save_imagej(f"{base}.tif", arr, acq.px_um, acq.dz_um, cols, rng, acq.channel_names)
            files.append(f"{base}.tif")
        return files

    def on_roi():
        if len(roi.data) == 0:
            show_warning(tr("n_no_roi"))
            return
        feats = roi.features
        files = []
        for k, s in enumerate(roi.data):
            yx = np.asarray(s)[:, -2:]
            fw = float(feats["fw"].iloc[k]) if k < len(feats) else 0
            fh = float(feats["fh"].iloc[k]) if k < len(feats) else 0
            if fw > 0 and fh > 0:  # fixed box: exact size
                r = fixed_rect(*yx.mean(0), int(fh), int(fw))
                y0, x0 = r[0].astype(int)
                y1, x1 = r[2].astype(int)
                lab = str(feats["label"].iloc[k]).replace("×", "x").replace(" ", "").replace("µ", "u")
                name = f"ROI{k + 1}_{lab}"
            else:
                y0, x0 = np.floor(yx.min(0)).astype(int)
                y1, x1 = np.ceil(yx.max(0)).astype(int)
                name = f"ROI{k + 1}"
            files += export_region(y0, y1, x0, x1, name)
        if files:
            show_info(tr("n_exported", n=len(files), p=out_dir))

    def on_all():
        H, W = levels[0].shape[-2:]
        if export_region(0, H, 0, W, "full"):
            show_info(tr("n_exported_all", p=out_dir))

    def on_shot():
        out_dir.mkdir(exist_ok=True)
        z, ts = tag()
        p = out_dir / f"screenshot_{z}_{ts}.png"
        viewer.screenshot(str(p), canvas_only=True, flash=True, scale=2)
        show_info(tr("n_shot", f=p.name))

    def on_open():
        out_dir.mkdir(exist_ok=True)
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(out_dir)))

    bt_roi.clicked.connect(on_roi)
    bt_all.clicked.connect(on_all)
    bt_shot.clicked.connect(on_shot)
    bt_open.clicked.connect(on_open)

    viewer.reset_view()
    viewer._tv = dict(acq=acq, set_mip=set_mip, export_region=export_region, on_auto=on_auto,
                      ch_layers=ch_layers, roi=roi, on_roi=on_roi, on_shot=on_shot,
                      add_fixed=add_fixed, snap_fixed=snap_fixed, sp_w=sp_w, sp_h=sp_h, cb_unit=cb_unit)
    return viewer


def main(argv):
    global LANG
    args, flags, lang = [], set(), None
    it = iter(argv[1:])
    for a in it:
        if a == "--lang":
            lang = next(it, None)
        elif a.startswith("--lang="):
            lang = a.split("=", 1)[1]
        elif a.startswith("--"):
            flags.add(a)
        else:
            args.append(a)
    LANG = detect_lang(lang)
    folder = args[0] if args else None
    if "--overview" in flags:
        acq = Acquisition(folder)
        m = build_mosaic(acq, force="--rebuild" in flags, flatfield="--raw" not in flags)
        print(tr("overview_saved"), export_overview(acq, m))
        return
    import napari
    v = run_viewer(folder)
    if v is not None:
        napari.run()


if __name__ == "__main__":
    main(sys.argv)
