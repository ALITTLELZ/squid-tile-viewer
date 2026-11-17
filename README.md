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
- Compare samples at identical display settings: open another acquisition with the current
  colors/contrast, side by side, optionally linked live; warns when exposure/laser/gain differ.
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

**File names and what is recorded with each export.**
- The *Export file name* field shows the default `<folder>_<channels>_<YYYYMMDD-HHMMSS>`, where
  `<channels>` is `merge` when two or more channels are shown and the wavelength (e.g. `561nm`)
  when one is. Type any name instead; *Default* switches back. Several ROIs get `_ROI1`, `_ROI2`;
  screenshots get `_screenshot`; an existing file is never overwritten (`_2`, `_3` … is added).
- Every exported image carries a description: acquisition, time, single plane or max projection,
  orientation, ROI position and size (px and µm), pixel size, scale bar, and for each channel
  shown/hidden, colour, contrast limits, gamma, opacity and imaging settings (exposure, laser
  intensity, gain), plus the flat-field method. Where to see it: TIFF — Fiji *Image > Show Info*;
  PNG — stored as the PNG *Description*/*Comment* (e.g. Preview *Tools > Show Inspector*).
- The same information is appended to `exports/export_log.csv`, one row per file, so a folder's
  exports can be sorted and filtered in Excel when assembling figures or slides.

**Comparing samples at the same contrast.** In *Compare samples*:
- *Open another acquisition (same colors & contrast)* opens a second folder in a new window with
  this window's per-channel colormap, contrast limits, gamma, opacity and visibility (channels are
  matched by name), and arranges the windows side by side. On a laptop screen each window's panel
  moves into a tab next to *layer controls* so both fit.
- *Copy colors & contrast to other windows* pushes the current settings once; tick *Link* to have
  every change in this window update the others live.
- If the two acquisitions were imaged with different exposure time, laser intensity or gain
  (read from `acquisition_channels.yaml`), you get a warning: the same contrast does not make
  brightness comparable in that case.

**Saved settings.** When you close the window (and after every export), the viewer saves the
per-channel contrast/gamma/colormap/opacity/visibility, MIP on/off, z-plane, orientation, view
position, ROIs and export options to `<acquisition>/tile_viewer_settings.json`, and restores them
the next time you open that folder. A folder you have never opened before gets the last colors and
contrast you used for channels with the same names, so similar samples look consistent.
Delete the JSON file to start fresh.

Language follows the system language; change it at the bottom of the panel (applies on restart).

## How stitching works, and how to tell whether it worked
The first time a folder is opened, these steps run once and are cached:

1. **Registration (where each tile goes).** For every pair of neighbouring tiles, the overlapping
   strips are compared (phase correlation on a max projection of all channels) to find the exact
   offset. A pair counts as *matched* if the strips agree well (correlation > 0.3). Pairs that
   don't match — usually empty glass or very weak signal — get the median offset of the matched
   pairs (the grid is regular). All offsets are then solved together (least squares) for one
   position per tile.
2. **Illumination (flat-field) correction.** The microscope lights the centre of the field more
   than its edges and corners, and not always symmetrically. The illumination profile is
   estimated from the data itself:
   - *shape*: the median of many tile images (each scaled to the same mean), smoothed with a 2D
     polynomial surface; once from images fully covered by tissue, once from all images with tissue;
   - *tilt*: the same tissue is imaged at the right edge of one tile and the left edge of the next
     (and top/bottom), so the brightness ratio there measures left/right and top/bottom
     asymmetry directly; a tilt is fitted to all overlaps;
   - the variant with the best agreement in the overlaps is used (shown in the QC figure).
   Each tile is then divided by this profile (after subtracting the camera offset).
3. **Blending.** In overlaps, the two tiles are mixed with weights that fade towards each tile's
   edge, so there is no hard line. The max projection is computed per tile first and then
   blended; blending first made seams dark in low-signal data (averaging two noisy images
   lowers their maximum). In low-signal data the maximum over z is also pushed up by noise;
   dividing by the illumination profile would amplify that push at the dim field corners. So
   the part of it needed to make the max projection's centre-to-corner profile match the
   z-average's is handled like the camera offset (automatically ~0 for bright samples).

**Stitching QC** (panel at the top, *Open QC figure*, *Show tile borders and numbers*):
- *Registration*: how many neighbour pairs were matched by image content and how far the final
  positions are from those measurements (✅ if median ≤ 3 px). Tiles without a matched neighbour
  are orange in the QC figure and in the tile-border overlay; they are placed on the grid, which
  is normally fine for empty areas but can't be verified.
- *Neighbour brightness ratio* (main check of the illumination correction): the same tissue seen
  in two neighbouring tiles should be equally bright after correction. ✅ if the median ratio is
  within 0.95–1.05 and the spread ≤ ±5%. The uncorrected value is shown for comparison.
- *Seam brightness*: brightness at the seams relative to right beside them, in a single plane and
  in the max projection (1 = invisible). Tissue structure also changes across a seam, so with fewer
  than 20 seams this is only indicative (ℹ️).
- What correction can't fix: differences in *focus* between tiles (a slice that isn't flat puts
  different depths in the same z-plane), photobleaching, or real differences in staining. These
  show up as single tiles that look brighter/dimmer in one z-plane; the max projection is less
  affected.

## Notes
- Expected folder layout: `ome_tiff/R0_<i>.ome.tiff` (Z × C × Y × X per tile) and
  `coordinates.csv` with the stage position of tile *i* on row *i*. Folders from aborted
  acquisitions (missing tiles) still open, with gaps where tiles are missing.
- Pixel size and z-step come from `acquisition.yaml`, otherwise from the OME metadata.
- Brightness in exported PNGs reflects your display settings. For quantification, use the
  16-bit TIFF and keep identical settings across samples. Note that the stitched data is
  flat-field corrected and blended in tile overlaps, so for strict intensity measurements
  go back to the raw tiles in `ome_tiff/`.
- The channel and ROI layers are locked against deletion (napari would otherwise delete the
  selected layer when you press delete in the layer list). To delete an ROI box, click it on
  the image and press delete.
- Tested on macOS (Apple Silicon) with napari 0.10. Windows launcher not yet tested.

## License
MIT, see [LICENSE](LICENSE).
