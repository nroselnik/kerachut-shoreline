# Vegetation line from Google Earth historical imagery

Google Earth has the only free sub-metre imagery of Pantai Kerachut going back to 2002. There is no API for it, so the images are captured as screenshots at a fixed view and georeferenced afterwards. Google Earth imagery is © Google and its data providers. It is not redistributed in this repository, so capture your own snapshots.

## 1. Capture protocol

1. Open Google Earth web at this camera position (top-down, north up, heading 0, tilt 0):
   `https://earth.google.com/web/@5.45118973,100.18285,3.4243048a,1105.37373818d,35y,0h,0t,0r`
2. Optional, to reduce label clutter: press Ctrl+B for basemap settings and choose **Clean**. Some place labels may stay on; mask them (see `LABELS` in `run_vegetation_line.py`).
3. Turn on **Historical Imagery** (clock icon in the toolbar). Use the arrows next to the date to step through the dates.
4. Wait about 5 seconds for each date to load fully, then take a screenshot of the browser window at 800 × 600 pixels. Do not pan or zoom between dates.
5. Save each one as `GE_YYYY-MM-DD.jpg`, using the date shown in the toolbar.
6. Check the bottom bar ("Data attribution ... older-date"). Some dates are mosaics, where the beach still comes from an older image. Skip a date if the beach looks identical to an earlier one. Skip cloudy dates too.

In our run, 19 dates were usable: 2002-06-21, 2005-01-23, 2010-02-20, 2011-02-15, 2014-02-21, 2015-02-08, 2015-10-02, 2016-01-12, 2016-02-13, 2016-12-29, 2017-08-19, 2017-09-07, 2018-07-04, 2019-09-01, 2019-12-23, 2020-01-22, 2022-01-13, 2024-02-08 and 2026-01-13.

If your screenshots are a different size, change `centre_px` in `georeference()` and the `valid` and `LABELS` boxes.

## 2. Georeferencing

`georeference()` fits the metres per pixel and the offset of the fixed view. It maximises the overlap between the water mask of one snapshot and the water mask of a Sentinel-2 scene of a similar date. We used the snapshot of 13 Jan 2026 and Sentinel-2 of 5 Jan 2026, which gave 1.155 m per pixel and an overlap (IoU) of 0.95. All snapshots share the same camera, so one fit serves every date.

Alignment between dates was checked by matching the end of the jetty. It agreed within 3 to 4 m from 2010 onward. The 2002 and 2005 images are less certain, and the jetty itself differs in 2005.

## 3. Vegetation line

Along each transect, starting near the typical waterline and moving landward, the vegetation line is placed where image brightness drops below the midpoint between sand and vegetation brightness and stays below it for 5 m. Positions are in metres from the transect origin. They use the same reference frame as the waterline analysis, so dry beach width is simply waterline minus vegetation line.

Run `run_vegetation_line.py` after editing the paths at the top. With our snapshots it reproduces `data/vegetation_line_by_transect.csv` to within 0.05 m.
