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
#     "matplotlib",
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
    "log_mip": ("  最大投影 通道 {c}/{n} 完成", "  max projection channel {c}/{n} done"),
    "log_qc": ("  计算拼接质控指标 ...", "  computing stitching QC ..."),
    "log_rebuild": ("检测到旧版本的拼接缓存，按改进后的照明校正重新拼接一次",
                    "Found a cache made with the older method; re-stitching with the improved illumination correction"),
    # QC
    "grp_qc": ("拼接质控", "Stitching QC"),
    "bt_qc_fig": ("查看质控图", "Open QC figure"),
    "cb_tiles": ("显示 tile 边界和编号", "Show tile borders and numbers"),
    "qc_reg": ("配准：{g}/{n} 对相邻 tile 按图像内容对齐，位置误差中位 {med:.1f} px（最大 {mx:.1f} px）",
               "Registration: {g}/{n} neighbour pairs matched by image content; position error median "
               "{med:.1f} px (max {mx:.1f} px)"),
    "qc_reg_none": ("配准：{n} 对相邻 tile 都无法按图像内容对齐（组织太少或信号太弱），全部按网格位置放置，无法验证对齐",
                    "Registration: none of the {n} neighbour pairs could be matched by image content (too little "
                    "tissue or signal); all tiles placed on the grid, alignment not verifiable"),
    "qc_fallback": ("  其中 {n}/{t} 个 tile 没有可对齐的邻居（多为空白区域），按网格位置放置",
                    "  {n}/{t} tiles had no matchable neighbour (mostly empty areas) and were placed on the grid"),
    "qc_agree": ("{name}：同一块组织在相邻两个 tile 里亮度比 {m:.3f}，离散 ±{s:.1f}%（未校正时 {rm:.3f}，±{rs:.1f}%）",
                 "{name}: same tissue in neighbouring tiles, brightness ratio {m:.3f}, spread ±{s:.1f}% "
                 "(uncorrected {rm:.3f}, ±{rs:.1f}%)"),
    "qc_agree_na": ("{name}：重叠区组织太少，无法检验照明校正", "{name}: too little tissue in the overlaps to check the illumination correction"),
    "qc_seam": ("{name}：接缝处亮度 / 两侧 = 单层 {z:.3f}，最大投影 {m:.3f}（{n} 条接缝）",
                "{name}: brightness at seams / beside = single plane {z:.3f}, max projection {m:.3f} ({n} seams)"),
    "qc_help": ("✅ 通过　⚠️ 需要留意　❌ 有问题　ℹ️ 仅供参考。照明校正看“相邻 tile 亮度比”：在 0.95–1.05 且离散不超过 ±5% 为通过。"
                "接缝亮度在接缝少于 20 条时受组织结构影响大，只作参考。",
                "✅ pass  ⚠️ check  ❌ problem  ℹ️ for information. Illumination correction: the neighbour brightness "
                "ratio should be within 0.95–1.05 with spread ≤ ±5%. Seam brightness is only indicative with fewer "
                "than 20 seams (tissue structure dominates)."),
    "qc_missing": ("（这个缓存没有质控数据）", "(no QC data in this cache)"),
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
    "fixed_tip": ("添加后直接拖动框内部即可移动；固定框的尺寸是锁定的，误拖边角也会自动恢复。"
                  "在图像上点选一个框后按 delete 键删除这个框。",
                  "Drag inside the box to move it. Its size is locked: accidental resizing snaps back. "
                  "Click a box on the image, then press delete to remove that box."),
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
    "grp_orient": ("方向（翻转 / 旋转）", "Orientation (flip / rotate)"),
    "flip_h": ("水平翻转 ⇆", "Flip horizontal ⇆"),
    "flip_v": ("垂直翻转 ⇅", "Flip vertical ⇅"),
    "rot_l": ("向左旋转 90° ⟲", "Rotate 90° left ⟲"),
    "rot_r": ("向右旋转 90° ⟳", "Rotate 90° right ⟳"),
    "orient_reset": ("恢复原始方向", "Reset orientation"),
    "orient_cur": ("当前：", "Current: "),
    "orient_none": ("原始方向", "original"),
    "d_flip": ("水平翻转", "flipped horizontally"),
    "d_rot1": ("逆时针旋转 90°", "rotated 90° counter-clockwise"),
    "d_rot2": ("旋转 180°", "rotated 180°"),
    "d_rot3": ("顺时针旋转 90°", "rotated 90° clockwise"),
    "autosave_note": ("对比度、颜色、方向、选框等显示设置会在关闭窗口（和每次导出）时自动保存，"
                      "下次打开同一文件夹自动恢复。",
                      "Display settings (contrast, colors, orientation, ROIs, ...) are saved automatically when "
                      "you close the window (and on every export) and restored next time you open this folder."),
    "grp_compare": ("多样本比较", "Compare samples"),
    "bt_open_other": ("打开另一个采集（沿用当前颜色和对比度）", "Open another acquisition (same colors & contrast)"),
    "bt_sync": ("把当前颜色和对比度同步到其他窗口", "Copy colors & contrast to other windows"),
    "cb_link": ("联动：调整时自动同步到其他窗口", "Link: changes here update other windows live"),
    "bt_tile": ("并排排列所有窗口", "Arrange windows side by side"),
    "compare_tip": ("通道按名称对应。只有成像条件（曝光、激光强度、增益）相同时，统一对比度后的亮度比较才有意义。",
                    "Channels are matched by name. Comparing brightness at the same contrast is only meaningful "
                    "if exposure, laser intensity and gain were the same."),
    "pick_other": ("选择要比较的另一个采集文件夹", "Choose another acquisition folder to compare"),
    "opening": ("正在打开 {name} …（第一次打开需要先拼接，请看终端进度）",
                "Opening {name} … (first opening stitches the tiles; see progress in the terminal)"),
    "n_inherited": ("已沿用「{src}」的颜色和对比度（{n} 个通道）", "Applied colors & contrast from '{src}' ({n} channels)"),
    "n_synced": ("已同步到 {n} 个窗口", "Copied to {n} window(s)"),
    "n_no_other": ("没有其他打开的窗口", "No other Tile Viewer windows are open"),
    "n_no_match": ("「{dst}」里没有与当前窗口同名的通道，未同步", "'{dst}' has no channels with matching names; nothing copied"),
    "w_imaging": ("注意：成像条件不同，亮度不能直接比较 —\n{diff}", "Warning: imaging settings differ, brightness is not directly comparable —\n{diff}"),
    "imaging_key": ({"exposure_ms": "曝光(ms)", "intensity": "激光强度", "gain": "增益"},
                    {"exposure_ms": "exposure (ms)", "intensity": "laser intensity", "gain": "gain"}),
    # notifications
    "n_restored": ("已恢复上次的显示设置", "Restored your last display settings for this folder"),
    "n_defaults": ("已沿用上次相同通道的颜色和对比度", "Applied your last color/contrast settings for these channels"),
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


# Every open Tile Viewer window registers itself here so settings can be copied between them.
OPEN_WINDOWS = []
_SYNC = {"busy": False}


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

        colors, imaging = {}, {}
        chy = _load_yaml(self.folder / "acquisition_channels.yaml") or {}
        for ch in (chy.get("channels") or acq.get("channels") or []):
            if not (isinstance(ch, dict) and ch.get("name")):
                continue
            if ch.get("display_color"):
                colors[ch["name"]] = ch["display_color"]
            cam = ch.get("camera_settings") or {}
            ill = ch.get("illumination_settings") or {}
            imaging[ch["name"]] = {"exposure_ms": cam.get("exposure_time_ms"), "gain": cam.get("gain_mode"),
                                   "intensity": ill.get("intensity")}

        def lookup(d, n):
            if n in d:
                return d[n]
            for k, v in d.items():
                if k in n or n in k:
                    return v
            return None

        self.channel_colors = [lookup(colors, n) or DEFAULT_COLORS[i % len(DEFAULT_COLORS)]
                               for i, n in enumerate(self.channel_names)]
        # imaging settings per channel (exposure, gain, laser intensity), used to warn when
        # comparing acquisitions that were not imaged the same way
        self.channel_imaging = {n: lookup(imaging, n) or {} for n in self.channel_names}

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
        # per pair: tile a, tile b, direction, match score (0 = unreliable), measured offset (y, x)
        "pairs": [[int(a_i), int(b_i), d, round(float(sc), 3), float(off[0]), float(off[1])]
                  for a_i, b_i, d, off, sc in results],
    }
    return pos, info


def _poly_surface(f, deg=6):
    """Smooth surface: least-squares fit of a 2D polynomial to log(f). Unlike Gaussian blurring,
    this does not pull the steep fall-off at the field edges towards the centre value."""
    h, w = f.shape
    yy, xx = np.mgrid[0:h, 0:w]
    yn, xn = (yy / (h - 1) * 2 - 1).ravel(), (xx / (w - 1) * 2 - 1).ravel()
    M = np.stack([yn ** i * xn ** j for i in range(deg + 1) for j in range(deg + 1 - i)], 1)
    coef, *_ = np.linalg.lstsq(M, np.log(np.clip(f, 1e-3, None)).ravel(), rcond=None)
    return np.exp(M @ coef).reshape(h, w)


def _overlap_pairs(pos, th, tw):
    """Neighbouring tile pairs (a, b, dy, dx) in full-resolution pixels."""
    out = []
    for a in sorted(pos):
        for b in sorted(pos):
            dy, dx = pos[b][0] - pos[a][0], pos[b][1] - pos[a][1]
            if (abs(dy) < th * 0.2 and tw * 0.5 < dx < tw) or (abs(dx) < tw * 0.2 and th * 0.5 < dy < th):
                out.append((a, b, dy, dx))
    return out


def estimate_flatfield(acq: Acquisition, log=print, ds=8, stats=None, pos=None):
    """Estimate each channel's illumination profile (flat field) from the tiles themselves.

    1. camera offset ("dark") = 0.5th percentile of all pixels.
    2. Shape: median of many images (each scaled to mean 1), smoothed with a 2D polynomial
       surface. Two versions: from images fully covered by tissue (an image that is half empty
       glass is dark at one edge, which would look like uneven illumination), and from all
       images with tissue.
    3. Tilt from the overlaps (needs tile positions): the same tissue is imaged at the left edge
       of one tile and the right edge of its neighbour (top/bottom likewise), so the brightness
       ratio there measures the left/right and top/bottom illumination ratio directly. Each shape
       is optionally corrected by a fitted tilt (x, y, and x*y terms).
    4. The candidate with the best overlap agreement (smallest spread of the brightness ratio of
       the same tissue in two tiles) is used. Its agreement is reported in the QC.
    Returns (flats[C,H,W] normalised to mean 1, darks[C]); fills `stats` (per channel) if given."""
    from scipy.ndimage import zoom
    ids = list(acq.tiles)
    pairs = _overlap_pairs(pos, acq.th, acq.tw) if pos else []
    flats, darks = [], []
    for c in range(acq.nc):
        stack = {i: acq.tile(i)[:, c, ::ds, ::ds].astype(np.float32) for i in ids}  # (Z, h, w) per tile
        imgs = np.concatenate([stack[i] for i in ids])
        n, h, w = imgs.shape
        dark = float(np.percentile(imgs, 0.5))
        # tissue coverage: fraction of 8x8 blocks clearly above background
        thr = dark + 0.25 * (np.percentile(imgs, 90) - dark)
        bh, bw = h // 8, w // 8
        blocks = imgs[:, :bh * 8, :bw * 8].reshape(n, 8, bh, 8, bw).transpose(0, 1, 3, 2, 4).reshape(n, 64, -1)
        cov = (np.median(blocks, -1) > thr).mean(1)
        n_full = int((cov >= 0.98).sum())
        means = imgs.mean(axis=(1, 2))

        def shape_from(sel):
            k = imgs[sel]
            k = (k - dark) / np.maximum(k.mean(axis=(1, 2), keepdims=True) - dark, 1e-3)
            f = _poly_surface(np.median(k, axis=0))
            return f / f.mean()

        shapes = {"covered": shape_from(cov >= 0.98 if n_full >= 10 else np.argsort(cov)[-max(10, n // 10):]),
                  "all": shape_from(means >= np.median(means))}

        # overlap data (z-mean, dark subtracted) for agreement and tilt fitting
        yy, xx = np.mgrid[0:h, 0:w]
        yn, xn = yy / (h - 1) * 2 - 1, xx / (w - 1) * 2 - 1
        ov = []
        for a, b, dy, dx in pairs:
            dy, dx = int(round(dy / ds)), int(round(dx / ds))
            ya0, ya1, xa0, xa1 = max(0, dy), min(h, h + dy), max(0, dx), min(w, w + dx)
            sa, sb = (slice(ya0, ya1), slice(xa0, xa1)), (slice(ya0 - dy, ya1 - dy), slice(xa0 - dx, xa1 - dx))
            A, B = stack[a][:, sa[0], sa[1]].mean(0) - dark, stack[b][:, sb[0], sb[1]].mean(0) - dark
            t = 0.25 * np.percentile(np.r_[A.ravel(), B.ravel()], 90)
            m = (A > t) & (B > t)
            if m.sum() >= 50:
                ov.append((A[m], B[m], sa, sb, m))

        def agreement(f):
            r = np.array([np.median(B / f[sb][m]) / np.median(A / f[sa][m]) for A, B, sa, sb, m in ov])
            return (float(np.median(r)), float(np.percentile(np.abs(np.log(r)), 75))) if len(r) else (None, None)

        def with_tilt(f0, terms):
            rows, rhs = [], []
            for A, B, sa, sb, m in ov:
                rhs.append(np.log(np.median(B / f0[sb][m]) / np.median(A / f0[sa][m])))
                xa, ya, xb, yb = xn[sa][m].mean(), yn[sa][m].mean(), xn[sb][m].mean(), yn[sb][m].mean()
                row = [xb - xa, yb - ya]
                if terms == 3:
                    row.append(xb * yb - xa * ya)
                rows.append(row)
            coef, *_ = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)
            g = coef[0] * xn + coef[1] * yn + (coef[2] * xn * yn if terms == 3 else 0)
            f = f0 * np.exp(g)
            return f / f.mean()

        cands = {}
        for sname, f0 in shapes.items():
            cands[sname] = f0
            if len(ov) >= 6:
                cands[sname + "+tilt"] = with_tilt(f0, 2)
                cands[sname + "+tilt+diag"] = with_tilt(f0, 3)
        scored = []
        for name, f in cands.items():
            if f.min() < 0.1 or f.max() > 5:
                continue  # implausible illumination profile
            med, spread = agreement(f)
            scored.append((spread if spread is not None else 0.0, name, f, med))
        if not scored:
            scored = [(0.0, "covered", shapes["covered"], None)]
        scored.sort(key=lambda s: s[0])
        spread, method, f, med = scored[0]
        raw_med, raw_spread = agreement(np.ones((h, w)))
        f = np.clip(f, 0.05, 20)

        # Max projection: the maximum over z is pushed up by noise. Dividing by the flat field
        # would amplify that push at the dim field edges, making tile corners bright in the MIP of
        # low-signal data. Part of the noise push (alpha * bias) is therefore treated like the camera
        # offset (not amplified); alpha is chosen so the MIP's centre-to-corner profile matches the
        # z-average's, which has no noise push.
        mx = np.stack([stack[i].max(0) for i in ids])
        mn = np.stack([stack[i].mean(0) for i in ids])
        bg = mn < np.percentile(mn, 5)
        bias = float(np.median((mx - mn)[bg]))
        tile_level = np.median(mn, axis=(1, 2))
        tissue = tile_level >= np.median(tile_level)
        kk = max(2, h // 20)

        def corner_ratio(img):
            img = np.median(img / np.maximum(np.median(img, axis=(1, 2), keepdims=True), 1e-6), 0)
            corners = np.mean([img[:kk, :kk].mean(), img[:kk, -kk:].mean(), img[-kk:, :kk].mean(), img[-kk:, -kk:].mean()])
            return corners / max(img[h // 2 - kk:h // 2 + kk, w // 2 - kk:w // 2 + kk].mean(), 1e-6)

        target = corner_ratio((mn[tissue] - dark) / f)
        err = lambda a: corner_ratio((mx[tissue] - dark - a * bias) / f + a * bias) - target
        if err(0.0) * err(1.0) > 0:
            alpha = 0.0 if abs(err(0.0)) <= abs(err(1.0)) else 1.0
        else:
            lo_, hi_ = 0.0, 1.0
            for _ in range(20):
                m_ = (lo_ + hi_) / 2
                lo_, hi_ = (lo_, m_) if err(lo_) * err(m_) <= 0 else (m_, hi_)
            alpha = (lo_ + hi_) / 2
        mip_bias = float(alpha * bias)
        del mx, mn
        F = zoom(f, (acq.th / f.shape[0], acq.tw / f.shape[1]), order=1)[: acq.th, : acq.tw]
        flats.append(F.astype(np.float32))
        darks.append(dark)
        if stats is not None:
            stats.append(dict(n_images=n, n_full=n_full, method=method, n_overlaps=len(ov),
                              agree_median=med, agree_spread=spread if ov else None,
                              raw_agree_median=raw_med, raw_agree_spread=raw_spread,
                              mip_bias=round(mip_bias, 2), mip_noise=round(bias, 2)))
        log(tr("log_flat", c=c, a=f.min(), b=f.max()) + (f"  [{method}]" if ov else ""))
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


def _correct(t, c, flats, darks):
    """Flat-field correct one tile plane (or stack) of channel c."""
    t = t.astype(np.float32)
    if flats is None:
        return t
    return (t - darks[c]) / flats[c] + darks[c]


def fuse_mip(acq, pos, flats, darks, out, log=print, mip_bias=None):
    """Maximum projection, computed per tile and then blended.
    Blending first and projecting afterwards made seams dark in low-signal data: averaging two
    noisy tiles lowers the noise, so the maximum over z is lower in the overlap than elsewhere.
    mip_bias[c]: part of the noise push of the maximum that is not amplified by the flat field
    (see estimate_flatfield)."""
    H, W = out.shape[-2:]
    fw = _feather(acq.th, acq.tw)
    wsum = np.zeros((H, W), np.float32)
    for i, (y, x) in pos.items():
        wsum[y:y + acq.th, x:x + acq.tw] += fw
    wsum[wsum == 0] = 1
    canvas = np.zeros((H, W), np.float32)
    for c in range(acq.nc):
        canvas[:] = 0
        extra = (mip_bias[c] if mip_bias else 0.0) if flats is not None else 0.0
        for i, (y, x) in pos.items():
            m = acq.tile(i)[:, c].max(axis=0).astype(np.float32)
            if flats is not None:
                base = darks[c] + extra
                m = (m - base) / flats[c] + base
            canvas[y:y + acq.th, x:x + acq.tw] += m * fw
        out[c] = np.clip(canvas / wsum + 0.5, 0, 65535).astype(np.uint16)
        log(tr("log_mip", c=c + 1, n=acq.nc))
    out.flush()


def _seam_profiles(acq, pos, img2d_fn, dark, max_pairs=60):
    """Brightness across tile seams of a fused image, compared with the bands right next to the
    overlap. img2d_fn(y0, y1, x0, x1) returns a region of the fused image.
    Returns (median per-seam ratio, number of seams, mean profile) - ratio 1.0 = no visible seam.
    Tissue structure varies across a seam too, so with few seams this is only indicative."""
    th, tw = acq.th, acq.tw
    ratios, profs = [], []
    for a, b, dy, dx in _overlap_pairs(pos, th, tw):
        horiz = abs(dy) < th * 0.2
        if horiz:
            ov = pos[a][1] + tw - pos[b][1]
            mid = (pos[b][1] + pos[a][1] + tw) // 2
            y0, y1 = max(pos[a][0], pos[b][0]), min(pos[a][0], pos[b][0]) + th
            x0, x1 = mid - int(1.1 * ov), mid + int(1.1 * ov)
            if x0 < 0:
                continue
            reg = img2d_fn(y0, y1, x0, x1)[::4]
            want = x1 - x0
        else:
            ov = pos[a][0] + th - pos[b][0]
            mid = (pos[b][0] + pos[a][0] + th) // 2
            x0, x1 = max(pos[a][1], pos[b][1]), min(pos[a][1], pos[b][1]) + tw
            y0, y1 = mid - int(1.1 * ov), mid + int(1.1 * ov)
            if y0 < 0:
                continue
            reg = img2d_fn(y0, y1, x0, x1)[:, ::4].T
            want = y1 - y0
        if reg.shape[1] != want:
            continue  # at the image border
        reg = reg.astype(np.float32) - dark
        rows = reg.mean(1)
        keep = rows > np.percentile(rows, 50)  # lines through tissue
        if keep.sum() < 10:
            continue
        p = np.median(reg[keep], axis=0)
        L = len(p)
        centre = p[int(L * 0.45):int(L * 0.55)].mean()
        side = 0.5 * (p[int(L * 0.05):int(L * 0.2)].mean() + p[int(L * 0.8):int(L * 0.95)].mean())
        if side <= 0:
            continue
        ratios.append(centre / side)
        profs.append(np.interp(np.linspace(0, L - 1, 111), np.arange(L), p / side))
        if len(ratios) >= max_pairs:
            break
    if not ratios:
        return None, 0, None
    return float(np.median(ratios)), len(ratios), np.mean(profs, 0).round(4).tolist()


def compute_qc(acq, pos, info, flats, darks, L0, M0, log=print):
    """Numbers that tell whether stitching and flat-field correction worked."""
    log(tr("log_qc"))
    pairs = info.get("pairs", [])
    good = [p for p in pairs if p[3] > 0.3]
    resid = [float(np.hypot(pos[p[1]][0] - pos[p[0]][0] - p[4], pos[p[1]][1] - pos[p[0]][1] - p[5])) for p in good]
    tiles_with_good = {p[0] for p in good} | {p[1] for p in good}
    qc = dict(
        pairs_total=len(pairs), pairs_good=len(good),
        resid_median=float(np.median(resid)) if resid else None,
        resid_max=float(np.max(resid)) if resid else None,
        tiles_total=len(pos), tiles_fallback=sorted(set(pos) - tiles_with_good),
        channels=[],
    )
    zmid = acq.nz // 2
    for c, name in enumerate(acq.channel_names):
        d = darks[c] if darks is not None else float(np.percentile(np.asarray(M0[c, ::64, ::64]), 0.5))
        r_mip, n_mip, p_mip = _seam_profiles(acq, pos, lambda y0, y1, x0, x1: np.asarray(M0[c, y0:y1, x0:x1]), d)
        r_z, n_z, p_z = _seam_profiles(acq, pos, lambda y0, y1, x0, x1: np.asarray(L0[zmid, c, y0:y1, x0:x1]), d)
        f = flats[c] if flats is not None else None
        qc["channels"].append(dict(
            name=name, dark=round(float(d), 1),
            illum_min=round(float(f.min()), 3) if f is not None else None,
            illum_max=round(float(f.max()), 3) if f is not None else None,
            seam_mip=r_mip, seam_z=r_z, n_seams=n_z, profile_mip=p_mip, profile_z=p_z))
    return qc


def qc_verdict(qc):
    """List of (status, text) lines for the panel; status True = pass, None = check,
    False = problem, "info" = for information only."""
    out = []
    if qc["pairs_good"] == 0:
        out.append((None, tr("qc_reg_none", n=qc["pairs_total"])))
    else:
        ok = qc["resid_median"] is not None and qc["resid_median"] <= 3 and qc["resid_max"] <= 15
        out.append((ok, tr("qc_reg", g=qc["pairs_good"], n=qc["pairs_total"],
                            med=qc["resid_median"], mx=qc["resid_max"])))
        if qc["tiles_fallback"]:
            out.append(("info", tr("qc_fallback", n=len(qc["tiles_fallback"]), t=qc["tiles_total"])))
    for ch in qc["channels"]:
        if ch.get("agree_spread") is not None:
            dev, sp = abs(ch["agree_median"] - 1), ch["agree_spread"]
            ok = True if (dev <= 0.05 and sp <= 0.05) else (None if (dev <= 0.10 and sp <= 0.10) else False)
            out.append((ok, tr("qc_agree", name=ch["name"], m=ch["agree_median"], s=sp * 100,
                               rm=ch["raw_agree_median"], rs=ch["raw_agree_spread"] * 100)))
        elif ch.get("illum_min") is not None:
            out.append(("info", tr("qc_agree_na", name=ch["name"])))
        if ch.get("seam_z") is not None:
            worst = max(abs(v - 1) for v in (ch["seam_z"], ch["seam_mip"]) if v is not None)
            if ch.get("n_seams", 0) >= 20:
                ok = True if worst <= 0.05 else (None if worst <= 0.10 else False)
            else:
                ok = "info"  # too few seams: tissue structure dominates the number
            out.append((ok, tr("qc_seam", name=ch["name"], z=ch["seam_z"],
                               m=ch["seam_mip"] if ch["seam_mip"] is not None else float("nan"), n=ch["n_seams"])))
    return out


def make_qc_figure(acq, meta, flats, out_path):
    """Save a one-page QC figure: tile layout + registration, illumination maps, seam profiles."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    qc, pos = meta["qc"], {int(k): v for k, v in meta["positions"].items()}
    pairs = meta["info"].get("pairs", [])
    nc = acq.nc
    fig = plt.figure(figsize=(6 + 3.2 * nc, 9), dpi=110)
    gs = fig.add_gridspec(2, 1 + nc, width_ratios=[2] + [1] * nc)
    # 1) layout
    ax = fig.add_subplot(gs[:, 0])
    th, tw = acq.th, acq.tw
    fallback = set(qc["tiles_fallback"])
    for i, (y, x) in pos.items():
        ax.add_patch(Rectangle((x, y), tw, th, fill=True, alpha=0.18,
                               color="tab:orange" if i in fallback else "tab:green", lw=0))
        ax.add_patch(Rectangle((x, y), tw, th, fill=False, ec="0.4", lw=0.5))
        ax.text(x + tw / 2, y + th / 2, str(i), ha="center", va="center", fontsize=7, color="0.2")
    for a, b, d, sc, *_ in pairs:
        if a not in pos or b not in pos:
            continue
        ya, xa = pos[a][0] + th / 2, pos[a][1] + tw / 2
        yb, xb = pos[b][0] + th / 2, pos[b][1] + tw / 2
        ax.plot([xa, xb], [ya, yb], color="tab:green" if sc > 0.3 else "0.75", lw=2 if sc > 0.3 else 0.8)
    H, W = meta["shape"][-2:]
    ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("Registration: green line = matched by image content,\n"
                 "grey = placed from the grid; orange tile = no matched neighbour", fontsize=9)
    # 2) illumination maps and 3) seam profiles
    for c, ch in enumerate(qc["channels"]):
        ax1 = fig.add_subplot(gs[0, 1 + c])
        if flats is not None:
            im = ax1.imshow(flats[c][::8, ::8], cmap="magma", vmin=0)
            fig.colorbar(im, ax=ax1, fraction=0.046, pad=0.04)
        ax1.set_title(f"{ch['name']}\nillumination (corrected for)", fontsize=9)
        ax1.set_xticks([]); ax1.set_yticks([])
        ax2 = fig.add_subplot(gs[1, 1 + c])
        x = np.linspace(-1.1, 1.1, 111)
        for key, lab, col in (("profile_z", "single z plane", "tab:blue"), ("profile_mip", "max projection", "tab:red")):
            if ch[key]:
                ax2.plot(x, ch[key], color=col, label=f"{lab}: {ch['seam_' + key.split('_')[1]]:.3f}")
        if ch.get("agree_median") is not None:
            ax1.set_xlabel(f"neighbour ratio {ch['agree_median']:.3f} ± {ch['agree_spread'] * 100:.1f}%  "
                           f"(uncorrected {ch['raw_agree_median']:.3f} ± {ch['raw_agree_spread'] * 100:.1f}%)\n"
                           f"method: {ch.get('method', '')}", fontsize=8)
        ax2.axhline(1, color="0.5", ls="--", lw=0.8)
        ax2.axvspan(-0.5, 0.5, color="0.9", zorder=0)
        ax2.set_ylim(0.6, 1.4)
        ax2.set_xlabel("distance from seam (overlap widths)", fontsize=8)
        ax2.set_title(f"brightness across seams, {ch.get('n_seams', 0)} seams (1 = no seam)", fontsize=9)
        ax2.legend(fontsize=7, loc="lower right")
    fig.suptitle(f"Stitching QC — {acq.folder.name}", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def _finish_qc(acq, cd, meta, pos, info, flats, darks, flat_stats, log):
    L0 = np.load(cd / "L0.npy", mmap_mode="r")
    M0 = np.load(cd / "M0.npy", mmap_mode="r")
    qc = compute_qc(acq, pos, info, flats, darks, L0, M0, log)
    for ch, st in zip(qc["channels"], flat_stats or []):
        ch.update(st)
    meta["qc"] = qc
    try:
        make_qc_figure(acq, meta, flats, cd / "qc_report.png")
    except Exception as e:  # the figure is a convenience; never fail the build for it
        log(f"QC figure failed: {e}")


def build_mosaic(acq: Acquisition, log=print, force=False, flatfield=True):
    """Stitch and write the cache. Returns dict(levels=[ZCYX memmaps], mip=[CYX memmaps], meta)."""
    cd = acq.cache_dir
    meta_p = cd / "meta.json"
    # "-ff2": flat field from tissue-covered images + polynomial surface, MIP per tile.
    # Caches made with the older method ("-ff") have a different key and are rebuilt.
    key = acq.key() + ("-ff2" if flatfield else "-raw2")
    if not force and meta_p.exists():
        meta = json.loads(meta_p.read_text())
        if meta.get("key") == key:
            return _open_cache(cd, meta)
        if str(meta.get("key", "")).startswith(acq.key()):
            log(tr("log_rebuild"))

    cd.mkdir(exist_ok=True)
    for p in cd.glob("*.npy"):
        p.unlink()
    t0 = time.time()
    log(tr("log_first", n=len(acq.tiles), c=acq.nc, z=acq.nz))
    pos, info = compute_positions(acq, log)
    H = max(p[0] for p in pos.values()) + acq.th
    W = max(p[1] for p in pos.values()) + acq.tw
    log(tr("log_size", w=W, h=H, wmm=W * acq.px_um / 1000, hmm=H * acq.px_um / 1000))

    flat_stats = []
    flats, darks = estimate_flatfield(acq, log, stats=flat_stats, pos=pos) if flatfield else (None, None)
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
                canvas[y:y + acq.th, x:x + acq.tw] += _correct(acq.tile(i)[z, c], c, flats, darks) * fw
            L0[z, c] = np.clip(canvas / wsum + 0.5, 0, 65535).astype(np.uint16)
        log(tr("log_z", i=z + 1, n=acq.nz))
    L0.flush()

    M0 = np.lib.format.open_memmap(cd / "M0.npy", "w+", np.uint16, (acq.nc, H, W))
    fuse_mip(acq, pos, flats, darks, M0, log, mip_bias=[s["mip_bias"] for s in flat_stats] if flat_stats else None)

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
                info=info, shape=[acq.nz, acq.nc, H, W], built=_dt.datetime.now().isoformat(timespec="seconds"),
                mip_method="per_tile")
    del L0, M0, prevL, prevM
    _finish_qc(acq, cd, meta, pos, info, flats, darks, flat_stats, log)
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
def run_viewer(folder: str | None, inherit: dict | None = None):
    """Open a viewer window. inherit = dict(src=name, channels={channel name: settings}, imaging=...)
    to apply another window's display settings to this one."""
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

    # tile borders (hidden by default; toggled from the QC panel)
    meta = mosaic["meta"]
    tile_pos = {int(k): v for k, v in meta["positions"].items()}
    qc = meta.get("qc")
    fallback_tiles = set((qc or {}).get("tiles_fallback", []))
    tile_ids = sorted(tile_pos)
    tiles_layer = viewer.add_shapes(
        [np.zeros((4, 2))] * len(tile_ids), shape_type="rectangle", ndim=2, name=tr("cb_tiles"),
        scale=scale[1:], face_color=[0, 0, 0, 0], edge_width=4,
        edge_color=[[1, 0.6, 0.1, 1] if i in fallback_tiles else [0.3, 1, 0.3, 1] for i in tile_ids],
        features=pd.DataFrame({"tile": [str(i) for i in tile_ids]}),
        text={"string": "{tile}", "color": "white", "size": 10, "anchor": "center"}, visible=False,
    )

    roi = viewer.add_shapes(
        ndim=2, name=tr("roi_layer"), scale=scale[1:], edge_color="white", face_color=[0, 0, 0, 0], edge_width=6,
        features=pd.DataFrame({"fw": pd.Series(dtype=float), "fh": pd.Series(dtype=float),
                               "label": pd.Series(dtype=str)}),
        text={"string": "{label}", "color": "white", "size": 11, "anchor": "upper_left"},
    )
    roi.feature_defaults = {"fw": 0.0, "fh": 0.0, "label": ""}
    roi.mode = "pan_zoom"
    viewer.scale_bar.visible = True
    for l in ch_layers + [roi, tiles_layer]:
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

    # flip/k: displayed image = rot90(fliplr?(data), k)  (k = counter-clockwise quarter turns)
    state = {"mip": False, "flip": False, "k": 0}
    out_dir = acq.folder / "exports"
    settings_file = acq.folder / "tile_viewer_settings.json"

    def orient(a):
        """Apply the current flip/rotation to the last two axes (returns a view, no copy)."""
        if state["flip"]:
            a = np.flip(a, -1)
        if state["k"]:
            a = np.rot90(a, state["k"], axes=(-2, -1))
        return a

    def cur_hw():
        h, w = levels[0].shape[-2:]
        return (w, h) if state["k"] % 2 else (h, w)

    def orient_pts(pts):
        """Map (y, x) points in original mosaic pixels to the displayed (flipped/rotated) frame."""
        H, W = levels[0].shape[-2:]
        y, x = pts[:, 0].astype(float), pts[:, 1].astype(float)
        if state["flip"]:
            x = W - x
        for _ in range(state["k"]):  # each counter-clockwise quarter turn
            y, x = W - x, y
            H, W = W, H
        return np.c_[y, x]

    def update_tiles_layer():
        rects = []
        for i in tile_ids:
            y, x = tile_pos[i]
            r = np.array([[y, x], [y, x + acq.tw], [y + acq.th, x + acq.tw], [y + acq.th, x]], float)
            o = orient_pts(r)
            (y0, x0), (y1, x1) = o.min(0), o.max(0)
            rects.append(np.array([[y0, x0], [y0, x1], [y1, x1], [y1, x0]]))
        for i, r in enumerate(rects):
            tiles_layer._data_view.edit(i, r)
        tiles_layer.refresh()
        tiles_layer.refresh_text()

    update_tiles_layer()

    # ---------------- panel ----------------
    w = QWidget()
    lay_v = QVBoxLayout(w)

    info = QLabel(f"<b>{acq.folder.name}</b><br>" + tr(
        "info", n=len(acq.tiles), c=acq.nc, z=acq.nz, px=acq.px_um, dz=acq.dz_um,
        w=levels[0].shape[-1], h=levels[0].shape[-2]))
    info.setWordWrap(True)
    lay_v.addWidget(info)

    gq = QGroupBox(tr("grp_qc"))
    vq = QVBoxLayout(gq)
    if qc:
        icon = {True: "✅", None: "⚠️", False: "❌", "info": "ℹ️"}
        txt = "<br>".join(f"{icon[ok]} {line}" for ok, line in qc_verdict(qc))
    else:
        txt = tr("qc_missing")
    lb_qc = QLabel(txt)
    lb_qc.setWordWrap(True)
    vq.addWidget(lb_qc)
    lb_qch = QLabel(tr("qc_help"))
    lb_qch.setWordWrap(True)
    lb_qch.setStyleSheet("color: gray")
    vq.addWidget(lb_qch)
    hq = QHBoxLayout()
    bt_qc = QPushButton(tr("bt_qc_fig"))
    bt_qc.setEnabled((acq.cache_dir / "qc_report.png").exists())
    cb_tiles = QCheckBox(tr("cb_tiles"))
    hq.addWidget(bt_qc)
    hq.addWidget(cb_tiles)
    vq.addLayout(hq)
    lay_v.addWidget(gq)

    g1 = QGroupBox(tr("grp_display"))
    v1 = QVBoxLayout(g1)
    cb_mip = QCheckBox(tr("mip"))
    v1.addWidget(cb_mip)
    bt_auto = QPushButton(tr("auto"))
    v1.addWidget(bt_auto)
    bt_reset = QPushButton(tr("reset"))
    v1.addWidget(bt_reset)
    lay_v.addWidget(g1)

    go = QGroupBox(tr("grp_orient"))
    vo = QVBoxLayout(go)
    ho1, ho2 = QHBoxLayout(), QHBoxLayout()
    bt_fh, bt_fv = QPushButton(tr("flip_h")), QPushButton(tr("flip_v"))
    bt_rl, bt_rr = QPushButton(tr("rot_l")), QPushButton(tr("rot_r"))
    ho1.addWidget(bt_fh)
    ho1.addWidget(bt_fv)
    ho2.addWidget(bt_rl)
    ho2.addWidget(bt_rr)
    vo.addLayout(ho1)
    vo.addLayout(ho2)
    bt_oreset = QPushButton(tr("orient_reset"))
    vo.addWidget(bt_oreset)
    lb_orient = QLabel()
    lb_orient.setWordWrap(True)
    lb_orient.setStyleSheet("color: gray")
    vo.addWidget(lb_orient)
    lay_v.addWidget(go)

    gc = QGroupBox(tr("grp_compare"))
    vc = QVBoxLayout(gc)
    bt_other = QPushButton(tr("bt_open_other"))
    bt_sync = QPushButton(tr("bt_sync"))
    cb_link = QCheckBox(tr("cb_link"))
    bt_tile = QPushButton(tr("bt_tile"))
    for wdg in (bt_other, bt_sync, cb_link, bt_tile):
        vc.addWidget(wdg)
    lb_cmp = QLabel(tr("compare_tip"))
    lb_cmp.setWordWrap(True)
    lb_cmp.setStyleSheet("color: gray")
    vc.addWidget(lb_cmp)
    lay_v.addWidget(gc)

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
    note = QLabel(tr("autosave_note"))
    note.setWordWrap(True)
    note.setStyleSheet("color: gray")
    lay_v.addWidget(note)

    hl = QHBoxLayout()
    hl.addWidget(QLabel(tr("lang_label")))
    cb_lang = QComboBox()
    cb_lang.addItems(["中文", "English"])
    cb_lang.setCurrentIndex(0 if LANG == "zh" else 1)
    hl.addWidget(cb_lang)
    lay_v.addLayout(hl)
    lay_v.addStretch(1)
    for g in (gq, g1, go, gc, g2, gf):
        g.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Minimum)
    scroll = QScrollArea()
    scroll.setWidgetResizable(True)
    scroll.setFrameShape(QFrame.Shape.NoFrame)
    scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
    scroll.setWidget(w)
    scroll.setMinimumWidth(w.sizeHint().width() + 20)
    tv_dock = viewer.window.add_dock_widget(scroll, name="Tile Viewer", area="right")

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

    def layer_data(c):
        src = mip_data(c) if state["mip"] else [lv[:, c] for lv in levels]
        return [orient(a) for a in src]

    def refresh_layers():
        """Swap layer data for the current MIP/orientation state, keeping display settings.
        If the shape changes (90° rotation), hide layers during the swap so napari never
        draws a frame with stale tile coordinates."""
        shape_changes = ch_layers[0].data[0].shape != layer_data(0)[0].shape
        vis = [l.visible for l in ch_layers]
        if shape_changes:
            for l in ch_layers:
                l.visible = False
        for c, l in enumerate(ch_layers):
            cl, cr = l.contrast_limits, l.contrast_limits_range
            l.data = layer_data(c)
            l.contrast_limits_range = cr
            l.contrast_limits = cl
        if shape_changes:
            for l, v in zip(ch_layers, vis):
                l.visible = v
        update_tiles_layer()

    def set_mip(on):
        if on == state["mip"]:
            return
        state["mip"] = on
        refresh_layers()

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
        """2D region (in displayed/oriented pixel coordinates) of channel c at the current z or MIP."""
        src = mips[level][c] if state["mip"] else levels[level][viewer.dims.current_step[0], c]
        return np.asarray(orient(src)[y0:y1, x0:x1])

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

    def fixed_rect(cy, cx, fh, fw):
        """fh x fw px rectangle centred on (cy, cx), shifted to stay inside the image."""
        H0, W0 = cur_hw()
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
        rot_tag = {1: "rot90ccw", 2: "rot180", 3: "rot90cw"}
        parts = (["flipH"] if state["flip"] else []) + ([rot_tag[state["k"]]] if state["k"] else [])
        if parts:
            z += "_" + "-".join(parts)
        return z, _dt.datetime.now().strftime("%Y%m%d-%H%M%S")

    def export_region(y0, y1, x0, x1, name):
        from PIL import Image
        out_dir.mkdir(exist_ok=True)
        H, W = cur_hw()
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
                arr = np.asarray(orient(levels[0])[:, :, y0:y1, x0:x1])
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
            save_settings()

    def on_all():
        H, W = cur_hw()
        if export_region(0, H, 0, W, "full"):
            show_info(tr("n_exported_all", p=out_dir))
            save_settings()

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

    # ---------------- orientation (flip / rotate) ----------------
    probe = np.arange(6).reshape(2, 3)
    ops = {"fliph": lambda a: np.flip(a, -1), "flipv": lambda a: np.flip(a, -2),
           "rotl": lambda a: np.rot90(a, 1), "rotr": lambda a: np.rot90(a, -1)}

    def canon(a):
        """Express any flip/rotation combination as (flip, k)."""
        for f in (False, True):
            for k in range(4):
                b = np.rot90(np.flip(probe, -1) if f else probe, k)
                if b.shape == a.shape and np.array_equal(b, a):
                    return f, k
        return False, 0

    def map_pts(op, pts, H, W):
        """Where points (y, x) of the current display land after applying op (H, W = display size)."""
        y, x = pts[:, 0], pts[:, 1]
        if op == "fliph":
            return np.c_[y, W - x]
        if op == "flipv":
            return np.c_[H - y, x]
        if op == "rotl":
            return np.c_[W - x, y]
        return np.c_[x, H - y]  # rotr

    def orient_desc():
        parts = ([tr("d_flip")] if state["flip"] else []) + ([tr(f"d_rot{state['k']}")] if state["k"] else [])
        return tr("orient_cur") + (" + ".join(parts) if parts else tr("orient_none"))

    def apply_orient_op(op):
        H, W = cur_hw()
        cur = np.rot90(np.flip(probe, -1) if state["flip"] else probe, state["k"])
        state["flip"], state["k"] = canon(ops[op](cur))
        center = list(viewer.camera.center)
        zoom = viewer.camera.zoom
        cyx = np.asarray(center[-2:]) / acq.px_um
        ny, nx = map_pts(op, cyx[None], H, W)[0]
        # move the ROIs with the image
        snapping["on"] = True
        try:
            for i, s in enumerate(list(roi.data)):
                roi._data_view.edit(i, map_pts(op, np.asarray(s)[:, -2:], H, W))
            if op in ("rotl", "rotr") and len(roi.data):
                f = roi.features.copy()
                f["fw"], f["fh"] = f["fh"].to_numpy().copy(), f["fw"].to_numpy().copy()
                f["label"] = [re.sub(r"^([\d.]+)×([\d.]+)", r"\2×\1", str(s)) for s in f["label"]]
                roi.features = f
            roi.refresh()
            roi.refresh_text()
        finally:
            snapping["on"] = False
        refresh_layers()
        center[-2:] = [ny * acq.px_um, nx * acq.px_um]
        viewer.camera.center = tuple(center)
        viewer.camera.zoom = zoom
        lb_orient.setText(orient_desc())

    def reset_orient():
        for _ in range(state["k"]):
            apply_orient_op("rotr")
        if state["flip"]:
            apply_orient_op("fliph")

    bt_fh.clicked.connect(lambda: apply_orient_op("fliph"))
    bt_fv.clicked.connect(lambda: apply_orient_op("flipv"))
    bt_rl.clicked.connect(lambda: apply_orient_op("rotl"))
    bt_rr.clicked.connect(lambda: apply_orient_op("rotr"))
    bt_oreset.clicked.connect(reset_orient)
    lb_orient.setText(orient_desc())

    # ---------------- save / restore display settings ----------------
    def channel_settings(l):
        cm = l.colormap
        return dict(contrast_limits=[float(v) for v in l.contrast_limits], gamma=float(l.gamma),
                    opacity=float(l.opacity), visible=bool(l.visible),
                    colormap=dict(name=cm.name, colors=np.asarray(cm.colors).tolist(),
                                  controls=np.asarray(cm.controls).tolist()))

    def collect_settings():
        f = roi.features
        rois = [dict(pts=np.asarray(s)[:, -2:].tolist(),
                     fw=float(f["fw"].iloc[i]) if i < len(f) else 0.0,
                     fh=float(f["fh"].iloc[i]) if i < len(f) else 0.0,
                     label=str(f["label"].iloc[i]) if i < len(f) else "")
                for i, s in enumerate(roi.data)]
        return dict(version=1, saved=_dt.datetime.now().isoformat(timespec="seconds"),
                    channels={l.name: channel_settings(l) for l in ch_layers},
                    mip=state["mip"], z=int(viewer.dims.current_step[0]),
                    flip=state["flip"], k=state["k"],
                    camera=dict(center=[float(v) for v in viewer.camera.center], zoom=float(viewer.camera.zoom)),
                    rois=rois,
                    export=dict(scale_bar=cb_bar.isChecked(), bar_um=sp_bar.value(),
                                tif=cb_tif.isChecked(), stack=cb_stack.isChecked()))

    def save_settings():
        try:
            s = collect_settings()
            settings_file.write_text(json.dumps(s, indent=1, ensure_ascii=False), encoding="utf-8")
            p = load_prefs()
            p.setdefault("channel_defaults", {}).update(s["channels"])
            p["export_opts"] = s["export"]
            save_prefs(p)
        except Exception as e:  # never block closing the window
            print("Could not save display settings:", e)

    def apply_channel(l, s):
        try:
            cm = s.get("colormap")
            if cm:
                l.colormap = Colormap(colors=cm["colors"], controls=cm["controls"], name=cm["name"])
            cl = s.get("contrast_limits")
            if cl:
                lo, hi = l.contrast_limits_range
                l.contrast_limits_range = (min(lo, cl[0]), max(hi, cl[1]))
                l.contrast_limits = tuple(cl)
            for k in ("gamma", "opacity", "visible"):
                if k in s:
                    setattr(l, k, s[k])
        except Exception as e:
            print(f"Could not restore settings for {l.name}:", e)

    def apply_export_opts(e):
        cb_bar.setChecked(e.get("scale_bar", True))
        sp_bar.setValue(int(e.get("bar_um", 0)))
        cb_tif.setChecked(e.get("tif", True))
        cb_stack.setChecked(e.get("stack", False))

    def restore_settings():
        s = None
        if settings_file.exists():
            try:
                s = json.loads(settings_file.read_text(encoding="utf-8"))
            except Exception as e:
                print("Could not read saved display settings:", e)
        if not s:
            cd = prefs.get("channel_defaults", {})
            used = [l for l in ch_layers if l.name in cd]
            for l in used:
                apply_channel(l, cd[l.name])
            apply_export_opts(prefs.get("export_opts", {}))
            if used:
                show_info(tr("n_defaults"))
            return
        state["flip"], state["k"] = bool(s.get("flip", False)), int(s.get("k", 0)) % 4
        state["mip"] = bool(s.get("mip", False))
        if state["flip"] or state["k"] or state["mip"]:
            refresh_layers()
        cb_mip.blockSignals(True)
        cb_mip.setChecked(state["mip"])
        cb_mip.blockSignals(False)
        lb_orient.setText(orient_desc())
        for l in ch_layers:
            if l.name in s.get("channels", {}):
                apply_channel(l, s["channels"][l.name])
        if "z" in s:
            viewer.dims.set_current_step(0, min(max(int(s["z"]), 0), acq.nz - 1))
        apply_export_opts(s.get("export", {}))
        for r in s.get("rois", []):
            roi.feature_defaults = {"fw": r.get("fw", 0.0), "fh": r.get("fh", 0.0), "label": r.get("label", "")}
            snapping["adding_fixed"] = True
            try:
                roi.add_rectangles(np.asarray(r["pts"], float))
            except Exception as e:
                print("Could not restore an ROI:", e)
            finally:
                snapping["adding_fixed"] = False
        roi.feature_defaults = {"fw": 0.0, "fh": 0.0, "label": ""}
        snapping["n"] = len(roi.data)
        roi.refresh_text()
        viewer.reset_view()
        cam = s.get("camera")
        if cam:
            try:
                viewer.camera.center = tuple(cam["center"])
                viewer.camera.zoom = cam["zoom"]
            except Exception:
                pass
        show_info(tr("n_restored"))

    # ---------------- compare samples: share settings between windows ----------------
    qwin = viewer.window._qt_window
    me = {"viewer": viewer, "qwin": qwin, "name": acq.folder.name, "imaging": acq.channel_imaging}

    def current_channels():
        return {l.name: channel_settings(l) for l in ch_layers}

    def imaging_diffs(src_imaging, names):
        keys = tr("imaging_key")
        out = []
        for n in names:
            a, b = (src_imaging or {}).get(n) or {}, acq.channel_imaging.get(n) or {}
            for k, label in keys.items():
                if a.get(k) is not None and b.get(k) is not None and a[k] != b[k]:
                    out.append(f"{n}: {label} {a[k]} ≠ {b[k]}")
        return out

    def apply_from(chans, src_imaging=None):
        """Apply {channel name: settings} to this window. Returns (channels applied, imaging differences)."""
        pairs = [(l, chans[l.name]) for l in ch_layers if l.name in chans]
        if not pairs and len(chans) == len(ch_layers):  # names differ: fall back to channel order
            pairs = list(zip(ch_layers, chans.values()))
        _SYNC["busy"] = True
        try:
            for l, s in pairs:
                apply_channel(l, s)
        finally:
            _SYNC["busy"] = False
        return len(pairs), imaging_diffs(src_imaging, [l.name for l, _ in pairs])

    me["apply"] = apply_from
    me["channels"] = current_channels
    OPEN_WINDOWS.append(me)

    def others():
        return [e for e in OPEN_WINDOWS if e is not me]

    def push_to_others(quiet=False):
        targets = others()
        if not targets:
            if not quiet:
                show_info(tr("n_no_other"))
            return
        chans, done, warn = current_channels(), 0, []
        for e in targets:
            try:
                n, diffs = e["apply"](chans, acq.channel_imaging)
            except Exception as ex:  # window may be closing
                print("Could not sync to", e.get("name"), ex)
                continue
            done += n > 0
            if n == 0 and not quiet:
                warn.append(tr("n_no_match", dst=e["name"]))
            if diffs and not quiet:
                warn.append(f"{e['name']}: " + tr("w_imaging", diff="\n".join(diffs)))
        if not quiet:
            show_info(tr("n_synced", n=done))
            for w_ in warn:
                show_warning(w_)

    def on_layer_change(event=None):
        if cb_link.isChecked() and not _SYNC["busy"]:
            push_to_others(quiet=True)

    for l in ch_layers:
        for ev in ("contrast_limits", "gamma", "colormap", "opacity", "visible"):
            getattr(l.events, ev).connect(on_layer_change)

    def on_link(checked):
        if checked:
            push_to_others()

    def compact():
        """Make the window narrower: move the Tile Viewer panel into a tab next to 'layer controls'."""
        from qtpy.QtWidgets import QDockWidget
        lc = [d for d in qwin.findChildren(QDockWidget) if d.windowTitle() == "layer controls"]
        if not lc or me.get("compact"):
            return
        qwin.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, tv_dock)
        qwin.tabifyDockWidget(lc[0], tv_dock)
        lc[0].raise_()
        me["compact"] = True

    me["compact_fn"] = compact

    def arrange_windows():
        """Tile all viewer windows side by side; on narrow screens, tab the panels so each fits."""
        g = QApplication.primaryScreen().availableGeometry()
        wdt = g.width() // max(len(OPEN_WINDOWS), 1)
        for i, e in enumerate(OPEN_WINDOWS):
            qw = e["qwin"]
            qw.showNormal()
            if qw.minimumSizeHint().width() > wdt:
                e["compact_fn"]()
                QApplication.processEvents()
            qw.setGeometry(g.x() + i * wdt, g.y(), wdt, g.height())
        QApplication.processEvents()
        for e in OPEN_WINDOWS:  # re-fit each image to its (now smaller) canvas
            e["viewer"].reset_view()

    def open_other():
        path = QFileDialog.getExistingDirectory(qwin, tr("pick_other"), str(acq.folder.parent))
        if not path:
            return
        save_settings()
        show_info(tr("opening", name=Path(path).name))
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        QApplication.processEvents()
        try:
            v2 = run_viewer(path, inherit=dict(src=acq.folder.name, channels=current_channels(),
                                               imaging=acq.channel_imaging))
        finally:
            QApplication.restoreOverrideCursor()
        if v2 is not None:
            arrange_windows()

    bt_other.clicked.connect(open_other)
    bt_sync.clicked.connect(lambda: push_to_others())
    cb_link.toggled.connect(on_link)
    bt_tile.clicked.connect(arrange_windows)

    from qtpy.QtCore import QEvent, QObject, QTimer

    class _SaveOnClose(QObject):
        def eventFilter(self, obj, ev):
            if ev.type() == QEvent.Type.Close:
                me["closing"] = True
                save_settings()
                if me in OPEN_WINDOWS:
                    OPEN_WINDOWS.remove(me)
            return False

    saver = _SaveOnClose(qwin)
    qwin.installEventFilter(saver)

    # ---------------- protect the channel and ROI layers from deletion ----------------
    # napari deletes the selected layer when the layer list has keyboard focus and you press
    # delete/backspace -- easy to trigger by accident (e.g. when trying to delete an ROI box).
    # 1) lock the layers against deletion (napari's own LayerLock: blocks the trash button,
    #    right-click "Delete" and Cmd+Backspace; does not affect contrast, colors or ROI editing);
    # 2) swallow delete/backspace in the layer list so nothing happens there at all.
    # (Re-inserting a deleted layer was tried and breaks its OpenGL drawing, so prevent instead.)
    bt_qc.clicked.connect(lambda: QDesktopServices.openUrl(QUrl.fromLocalFile(str(acq.cache_dir / "qc_report.png"))))
    cb_tiles.toggled.connect(lambda on: setattr(tiles_layer, "visible", on))

    for l in list(ch_layers) + [roi, tiles_layer]:
        try:
            l.locked = True
        except Exception:  # napari without layer locking
            pass

    class _NoDeleteKey(QObject):
        def eventFilter(self, obj, ev):
            if ev.type() == QEvent.Type.KeyPress and ev.key() in (Qt.Key.Key_Backspace, Qt.Key.Key_Delete):
                return True
            return False

    layer_list_view = getattr(viewer.window._qt_viewer, "layers", None)
    no_del = _NoDeleteKey(qwin)
    if layer_list_view is not None:
        layer_list_view.installEventFilter(no_del)

    viewer.reset_view()
    restore_settings()
    if inherit:
        n, diffs = apply_from(inherit.get("channels", {}), inherit.get("imaging"))
        if n:
            show_info(tr("n_inherited", src=inherit.get("src", "?"), n=n))
        else:
            show_warning(tr("n_no_match", dst=acq.folder.name))
        if diffs:
            show_warning(tr("w_imaging", diff="\n".join(diffs)))
    # Force a fresh slice. Without this, a folder opened for the first time (no saved settings,
    # so nothing else moves the z slider) stayed on napari's empty 1x1 placeholder slice and
    # showed a black canvas. layer.refresh() does not fix it; changing the z step does.
    z_now = viewer.dims.current_step[0]
    if acq.nz > 1:
        viewer.dims.set_current_step(0, (z_now + 1) % acq.nz)
        viewer.dims.set_current_step(0, z_now)
    else:
        for l in ch_layers:
            l.refresh()
    viewer._tv = dict(acq=acq, set_mip=set_mip, export_region=export_region, on_auto=on_auto,
                      ch_layers=ch_layers, roi=roi, on_roi=on_roi, on_shot=on_shot,
                      add_fixed=add_fixed, snap_fixed=snap_fixed, sp_w=sp_w, sp_h=sp_h, cb_unit=cb_unit,
                      apply_orient_op=apply_orient_op, reset_orient=reset_orient, state=state,
                      save_settings=save_settings, settings_file=settings_file, saver=saver,
                      plane=plane, cur_hw=cur_hw, me=me, push_to_others=push_to_others, cb_link=cb_link,
                      arrange_windows=arrange_windows, apply_from=apply_from, current_channels=current_channels,
                      layer_list_view=layer_list_view)
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
