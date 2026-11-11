# Tile Viewer

**English** | [简体中文](README.zh-CN.md)

Stitch and browse tiled fluorescence acquisitions from the Squid (Cephla) microscope,
adjust channel colors/contrast, and export composite images of fixed- or free-size regions.

## What it does
- Opens one acquisition folder (the one containing `ome_tiff/` and `coordinates.csv`).
- Stitches all tiles into a single mosaic using image registration (not just stage
  coordinates, which can be off by ~100 px), with flat-field correction to remove the
  dark tile borders (vignetting).
- Shows every channel as its own layer in [napari](https://napari.org): color, contrast,
  gamma and visibility per channel; z-slider or maximum-intensity projection (MIP).
- Flips and 90° rotations, so every section can be shown (and exported) in the same orientation.
- Remembers your display settings per folder (contrast, colors, orientation, ROIs, ...).
- Exports, into `<acquisition>/exports/`:
  - **ROIs** at full resolution: fixed-size boxes (enter width × height in µm or pixels,
    drag to position; the size stays locked) or free rectangles.
  - The **whole mosaic** at full resolution.
  - A **screenshot** of the current view.

  ROI and mosaic exports are saved as a colored composite PNG (optional scale bar) and,
  optionally, a 16-bit TIFF that opens in Fiji as a colored composite (for quantification).

## Install & start
Nothing to install by hand. The launcher uses [uv](https://docs.astral.sh/uv/) to download
Python and all packages into an isolated cache. It doesn't touch any existing
Python/conda setup.

**macOS:** double-click **`Tile Viewer.command`**.
- If macOS says it "cannot be opened because it is from an unidentified developer":
  right-click the file → **Open** → **Open** (only needed once).
- If it says it is not executable, run once in Terminal:
  `chmod +x "/path/to/Tile Viewer.command"`

**Windows:** double-click **`Tile Viewer (Windows).bat`**.

If uv isn't installed, the launcher asks before installing it from astral.sh.
The **first start downloads ~300 MB** and takes a few minutes; later starts are fast.
Requires an internet connection only for that first start.

Command line alternative (any OS with uv):
```
uv run --python 3.12 tile_viewer.py [ACQ_FOLDER] [--lang en|zh]
uv run --python 3.12 tile_viewer.py ACQ_FOLDER --overview   # just export a stitched overview PNG
```

## Using it
1. Pick an acquisition folder in the dialog. The first time you open a folder it's
   stitched (about 20 s for 12 tiles; progress shows in the terminal window). The result
   is cached in `<acquisition>/_viewer_cache/` (about the same size as the raw data;
   safe to delete, it's rebuilt on demand).
2. **Channels:** select a channel in the layer list (left) and use *layer controls*:
   *contrast limits* (brightness), *gamma*, *colormap*, *opacity*. The eye icon toggles a channel.
   Panel buttons: *Auto contrast (current view)* and *Reset contrast*.
3. **z / MIP:** use the slider under the image, or tick *Max intensity projection*.
4. **Orientation:** *Flip horizontal / vertical*, *Rotate 90° left / right*, *Reset orientation*.
   ROIs move with the image; exports (PNG, TIFF, scale bar) come out in the orientation you see, and the
   file name records it (e.g. `_flipH-rot90cw`). Nothing is copied: it's only a different view of the data.
5. **Fixed-size ROI:** enter W and H, choose µm or pixels, click *Add fixed-size box*.
   Drag inside the box to move it. Add as many as you like; select one and press Delete to remove it.
6. Click *Export ROIs* → files appear in `exports/` (button *Open export folder*).
   Exports use exactly the colors/contrast you see. Untick channels you don't want in the composite.

**Saved settings.** When you close the window (and after every export), the viewer saves the
per-channel contrast/gamma/colormap/opacity/visibility, MIP on/off, z-plane, orientation, view
position, ROIs and export options to `<acquisition>/tile_viewer_settings.json`, and restores them
the next time you open that folder. A folder you have never opened before gets the last colors and
contrast you used for channels with the same names, so similar samples look consistent.
Delete the JSON file to start fresh.

Language follows the system language; change it at the bottom of the panel (applies on restart).

## Notes
- Expected folder layout: `ome_tiff/R0_<i>.ome.tiff` (Z × C × Y × X per tile) and
  `coordinates.csv` with the stage position of tile *i* on row *i*. Folders from aborted
  acquisitions (missing tiles) still open, with gaps where tiles are missing.
- Pixel size and z-step come from `acquisition.yaml`, otherwise from the OME metadata.
- Brightness in exported PNGs reflects your display settings. For quantification, use the
  16-bit TIFF and keep identical settings across samples. Note that the stitched data is
  flat-field corrected and blended in tile overlaps, so for strict intensity measurements
  go back to the raw tiles in `ome_tiff/`.
- Tested on macOS (Apple Silicon) with napari 0.10. Windows launcher not yet tested.

## License
MIT, see [LICENSE](LICENSE).
