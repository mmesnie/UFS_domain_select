#!/usr/bin/env python3

#
# TODO / FIXME
#
# - resolve g_extent vs. g_lambert_extent. Do we need both??
# - fix filering for GFS (currently filters based on lon/lat spans)
# - clean up regional vs. global, esp. with settting write compoment when in global view
# - make sure GFS plots look ok, based on args we generated
#

#
# This script has two uses. First, it automates the creation of the YAML 
# configuration that specifies the compute grid and write component. The write
# component will be shown in red and the compute grid in blue.  Usage is -
#
g_help="""
******************************************************************************
*         Zoom: zoom grid                                                    *
*          Pan: shift grid                                                   *
*     Spacebar: set write component (view taken from Lambert Conformal plot) *
*        Key g: toggle global view                                           *
*        Key y: output yaml config                                           *
*        Key -: decrement Gnomonic compute grid by 1%                        *
*        Key +: increment Gnomonic compute grid by 1%                        *
*        Key R: set resolution (3km, 13km, 25km, AUTO)                       *
*        Key 0: restore default settings                                     *
*        Key h: show this help                                               *
* Middle Click: center grid                                                  *
******************************************************************************
"""

#
# The second use is to graph a grib file (500 MB geopotential and vorticity), using
# the -f <FILENAME> option. This can be a GFS grib file (initial or boundary conditions)
# or the output from a SRW forecast (i.e., the grib files in "postprd" directory). The
# center latitude and longitude can be specified via the --cen_lon and _cen_lat options,
# and the lower left cornet longitude and latitude via --crn_lon and --crn_lat. 
#
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--cen_lon", help="center longitude", required=False)
parser.add_argument("--cen_lat", help="center latitude", required=False)
parser.add_argument("--crn_lon", help="lower left corner longitude", required=False)
parser.add_argument("--crn_lat", help="lower left corner latitude", required=False)
parser.add_argument("--compute", "-c", help="show compute grid (Gnomonic)", required=False, action="store_true")
parser.add_argument("--file", "-f", help="grib file to plot", required=False)
parser.add_argument("--close", "-x", help="close after saving plot of grib file", required=False, action="store_true")
args = parser.parse_args()

#
# Imports
#
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.backend_bases import MouseButton
import time
import math
import io
import PIL
import pygrib
import scipy
import numpy as np
import cartopy
import matplotlib

g_output_filename = '/home/mmesnie/ufs-srweather-app-v2.2.0/ush/config.yaml'

#
# Static user config
#
g_dt_atmos = 36
g_blocksize = 40
g_layout_x = 3
g_layout_y = 3
g_write_groups = 1
g_write_tasks_per_group = 3
g_date ='20250403'
g_cycle ='00'
g_fcst_len_hrs = 6
g_lbc_spec_intvl_hrs = 6
g_extrn_mdl_source_basedir_ics = f"/home/mmesnie/DATA-2.2.0/input_model_data/FV3GFS/grib2/{g_date}{g_cycle}"
g_extrn_mdl_source_basedir_lbcs = f"/home/mmesnie/DATA-2.2.0/input_model_data/FV3GFS/grib2/{g_date}{g_cycle}"

#
# Forecast domain defaults (g_res of -1 is "auto" mode)
#
g_cen_lon_dflt=-59.5;   g_cen_lat_dflt=-51.7; g_crn_lon_dflt=-61.98;  g_crn_lat_dflt=-52.81 # Falkland Islands
g_cen_lon_dflt=-97.5;   g_cen_lat_dflt=38.5;  g_crn_lon_dflt=-122.72; g_crn_lat_dflt=21.14  # CONUS
g_cen_lon_dflt=-141.87; g_cen_lat_dflt=40.48; g_crn_lon_dflt=-160.29; g_crn_lat_dflt=16.64  # Eastern Pacific 
g_cen_lon_dflt=-127.68; g_cen_lat_dflt=45.72; g_crn_lon_dflt=-132.86; g_crn_lat_dflt=41.77  # Oregon coast


#
# Internal globals (don't modify)
#
g_lambert_extent = None
g_crn_lon = None
g_crn_lat = None
g_res_dflt=-1 #  Auto mode (picks largest resolution)
g_compute_grid_dflt = 0.05 # Fraction that compute grid is larger than the write compoment
if args.cen_lon:
    g_cen_lon_dflt = float(args.cen_lon)
if args.cen_lat:
    g_cen_lat_dflt = float(args.cen_lat)
if args.crn_lon:
    g_crn_lon_dflt = float(args.crn_lon)
if args.crn_lat:
    g_crn_lat_dflt = float(args.crn_lat)
g_dflt = "regional" # regional or global
g_cen_lon = g_cen_lon_dflt
g_cen_lat = g_cen_lat_dflt
g_res = g_res_dflt
g_view = {}
g_axis = {}
g_proj = {}
g_projs = {}
g_extent = {}
g_mode = None

#
# Function definitions
#

def cen_lat_adjust(cen_lat):
    if cen_lat > 80:
        print("MAX CENTER LATITUDE is 80")
        cen_lat = 80
    if cen_lat < -80:
        print("MIN CENTER LATITUDE is -80")
        cen_lat = -80
    if cen_lat == 0:
        # LambertConformal doesn't work with center latitudes of 0, so adjust slightly
        print("CENTER LATITUDE CANNOT BE ZERO - setting to 0.001 ")
        cen_lat = 0.001

    return cen_lat

def on_button_press(event):
    global g_cen_lon, g_cen_lat

    if event.inaxes and event.button is MouseButton.MIDDLE:
        g_cen_lon, g_cen_lat = ccrs.PlateCarree().transform_point(event.xdata, event.ydata, 
                                                                  event.inaxes.projection)
        g_cen_lat = cen_lat_adjust(g_cen_lat)
        plots_draw("center") #2

def fmt_tuple(tuple_in):
    return tuple([str(round(x,2)) if isinstance(x, float) else x for x in tuple_in])

def get_dims(index, color):
    extent = g_axis[index].get_extent()
    x1, x2, y1, y2 = extent

    if color == "red":
        xspan = abs(x2-x1)
        x1 -= xspan*g_compute_grid/2
        x2 += xspan*g_compute_grid/2
        yspan = abs(y2-y1)
        y1 -= yspan*g_compute_grid/2
        y2 += yspan*g_compute_grid/2
        extent = x1, x2, y1, y2

    xc = (x1+x2)/2
    yc = (y1+y2)/2
    xspan = abs(x2-x1)
    yspan = abs(y2-y1)
    cen_lon, cen_lat = ccrs.PlateCarree().transform_point(xc, yc, g_axis[index].projection)
    lwr_lon, lwr_lat = ccrs.PlateCarree().transform_point(x1, y1, g_axis[index].projection)
    upr_lon, upr_lat = ccrs.PlateCarree().transform_point(x2, y2, g_axis[index].projection)

    return cen_lon, cen_lat, lwr_lon, lwr_lat, extent

def pick_max_delta(xspan, yspan):
    if (xspan/25000 < 50) or (yspan/25000<50):
        if (xspan/13000 < 50) or (yspan/13000<50):
            return 3000
        else:
            return 13000
    else:
        return 25000

def output_config(index):
    index = "LambertConformal"
    res = g_res
    # computational grid
    cen_lon_r, cen_lat_r, lwr_lon_r, lwr_lat_r, extent_r = get_dims(index, "red")
    x1, x2, y1, y2 = extent_r
    xspan_r = abs(x2-x1)
    yspan_r = abs(y2-y1)
    max_delta = pick_max_delta(xspan_r, yspan_r)
    if g_res == -1:
        res = max_delta
        print(f"             resolution: automatically set to {res} meters")
    elif (g_res>max_delta):
        res = max_delta
        print(f"             resolution: forced to max ({res} meters)")
    percent = round(100*(1+g_compute_grid))
    print(f"      compute grid size: {percent}% of write component")
    nx_r = int(xspan_r/res)
    ny_r = int(yspan_r/res)

    print(f"  gnomonic compute grid: nx {nx_r} ny {ny_r}")
    # write component grid
    cen_lon_b, cen_lat_b, lwr_lon_b, lwr_lat_b, extent_b = get_dims(index, "blue")
    x1, x2, y1, y2 = extent_b
    xspan_b = abs(x2-x1)
    yspan_b = abs(y2-y1)
    nx_b = int(xspan_b/res)
    ny_b = int(yspan_b/res)

    print(f"lambert write component: nx {nx_b} ny {ny_b}")
    print(f"       center longitude: {g_cen_lon}")
    print(f"        center latitude: {g_cen_lat}")
    print(f"     lambert corner lon: {round(lwr_lon_b, 2)}")
    print(f"     lambert corner lat: {round(lwr_lat_b, 2)}")

    config_text = f"""
metadata:
  description: >-
    Automatically generated via UFS_domain_select.py
user:
  RUN_ENVIR: community
  MACHINE: linux
  ACCOUNT: an_account
workflow:
  EXPT_SUBDIR: /home/mmesnie/expt-2.2.0/test_community
  USE_CRON_TO_RELAUNCH: false
  CCPP_PHYS_SUITE: FV3_GFS_v16
  DATE_FIRST_CYCL: '{g_date}{g_cycle}'
  DATE_LAST_CYCL: '{g_date}{g_cycle}'
  FCST_LEN_HRS: {g_fcst_len_hrs}
  PREEXISTING_DIR_METHOD: rename
  VERBOSE: true
  COMPILER: gnu
  GRID_GEN_METHOD: "ESGgrid"
task_make_grid:
  ESGgrid_LON_CTR: {g_cen_lon}
  ESGgrid_LAT_CTR: {g_cen_lat}
  ESGgrid_NX: {nx_r} 
  ESGgrid_NY: {ny_r}
  ESGgrid_DELX: {res}
  ESGgrid_DELY: {res} 
  ESGgrid_PAZI: 0.0
  ESGgrid_WIDE_HALO_WIDTH: 6
task_get_extrn_ics:
  EXTRN_MDL_NAME_ICS: FV3GFS
  FV3GFS_FILE_FMT_ICS: grib2
  EXTRN_MDL_SOURCE_BASEDIR_ICS: {g_extrn_mdl_source_basedir_ics}
task_get_extrn_lbcs:
  EXTRN_MDL_NAME_LBCS: FV3GFS
  LBC_SPEC_INTVL_HRS: {g_lbc_spec_intvl_hrs} 
  FV3GFS_FILE_FMT_LBCS: grib2
  EXTRN_MDL_SOURCE_BASEDIR_LBCS: {g_extrn_mdl_source_basedir_lbcs}
task_run_fcst:
  QUILTING: true
  DT_ATMOS: {g_dt_atmos}
  LAYOUT_X: {g_layout_x} 
  LAYOUT_Y: {g_layout_y}
  BLOCKSIZE: {g_blocksize}
  WRTCMP_write_groups: {g_write_groups}
  WRTCMP_write_tasks_per_group: {g_write_tasks_per_group}
  WRTCMP_output_grid: "lambert_conformal"
  WRTCMP_cen_lon: {g_cen_lon}
  WRTCMP_cen_lat: {g_cen_lat} 
  WRTCMP_stdlat1: {g_cen_lat} 
  WRTCMP_stdlat2: {g_cen_lat} 
  WRTCMP_nx: {nx_b}
  WRTCMP_ny: {ny_b}
  WRTCMP_lon_lwr_left: {round(lwr_lon_b, 2)} 
  WRTCMP_lat_lwr_left: {round(lwr_lat_b, 2)}
  WRTCMP_dx: {res}
  WRTCMP_dy: {res}
task_run_post:
  POST_OUTPUT_DOMAIN_NAME: 'mesnier'
"""
    with open(g_output_filename, "w") as file:
        file.write(config_text)
    print(f"       YAML file output: {g_output_filename}")
    print(f"                    cmd: export DATE={g_date} CYCLE={g_cycle} LEN={g_fcst_len_hrs} LBC={g_lbc_spec_intvl_hrs}; time ./forecast")

def on_key_press(event):
    global g_cen_lon, g_cen_lat
    global g_compute_grid
    global g_res
    global g_extent

    ax = event.inaxes
    index = get_index(ax)

    if event.key == 'h':
        show_help()
    elif event.key == '0':
        print(f"restored default settings")
        plots_draw("init") #3
    elif event.key == 'y':
        index = "LambertConformal"
        output_config("LambertConformal")
    elif event.key == 'R':
        if g_res == 3000:
            g_res = 13000
        elif g_res == 13000:
            g_res = 25000
        elif g_res == 25000:
            g_res = -1 
        elif g_res == -1:
            g_res = 3000
        if (g_res == -1):
            print(f"resolution set to AUTO")
        else:
            print(f"resolution set to {int(g_res/1000)}km")
    elif event.key == '+':
        g_compute_grid += .01
        percent = round(100*(1+g_compute_grid))
        print(f"set compute grid to {percent}% of write component")
        plots_draw("set") #4
    elif event.key == '-':
        g_compute_grid -= .01
        percent = round(100*(1+g_compute_grid))
        print(f"set compute grid to {percent}% of write component")
        plots_draw("set") #5
    elif event.key == ' ':
        plots_draw("set") #6
    elif event.key == 'g':
        if index == "Gnomonic":
            print("global view disabled for Gnomonic")
        elif index:
            print(f"toggling global view for {index}...")
            if g_view[index] == "regional":
                g_extent[index] = g_axis[index].get_extent()
                g_view[index] = "global"
                g_axis[index].set_global()
            elif g_view[index] == "global": 
                g_view[index] = "regional"
                g_axis[index].set_extent(g_extent[index], crs=g_proj[index])
    elif event.key == 'q':
        exit(0)

def get_index(ax):
    for p in g_projs:
        if g_axis[p] == ax:
            return p

def projs_create(mode):
    global g_proj
    global g_projs
    global g_cen_lon, g_cen_lat
    global g_view
    global g_dim_x, g_dim_y
    global g_compute_grid
    global g_res

    #print(f"projs_create: mode {mode} g_cen_lon {g_cen_lon} g_cen_lat {g_cen_lat}")

    if mode == "init":
        g_cen_lon = g_cen_lon_dflt
        g_cen_lat = g_cen_lat_dflt
        g_compute_grid = g_compute_grid_dflt
        g_res = g_res_dflt

    g_cen_lat = cen_lat_adjust(g_cen_lat)

    g_proj = {}

    # Class: Cylindrical, Conic, Planar/Azimuthal
    # Case: Tangent, Secant
    # Aspect: Normal, Transverse, Oblique

    # Cylindrical
    g_proj["PlateCarree"] = ccrs.PlateCarree(central_longitude=g_cen_lon)
    g_proj["Mercator"] = ccrs.Mercator(central_longitude=g_cen_lon)
    g_proj["Miller"] = ccrs.Miller(central_longitude=g_cen_lon)

    # Conic
    g_proj["EquidistantConic"] = ccrs.EquidistantConic(central_longitude=g_cen_lon, central_latitude=g_cen_lat)
    g_proj["AlbersEqualArea"] = ccrs.AlbersEqualArea(central_longitude=g_cen_lon, central_latitude=g_cen_lat)
    if g_cen_lat>0:
        lat_cutoff = -30
    else:
        lat_cutoff = 30
    g_proj["LambertConformal"] = ccrs.LambertConformal(central_longitude=g_cen_lon, central_latitude=g_cen_lat,
                                                       standard_parallels=(g_cen_lat, g_cen_lat), cutoff=lat_cutoff)

    # Planar/Azimuthal
    g_proj["Stereographic"] = ccrs.Stereographic(central_longitude=g_cen_lon, central_latitude=g_cen_lat)
    g_proj["Orthographic"] = ccrs.Orthographic(central_longitude=g_cen_lon, central_latitude=g_cen_lat)
    g_proj["Gnomonic"] = ccrs.Gnomonic(central_longitude=g_cen_lon, central_latitude=g_cen_lat)

    # Pseudo-cylindrical
    g_proj["Robinson"] = ccrs.Robinson(central_longitude=g_cen_lon)
    g_proj["Mollweide"] = ccrs.Mollweide(central_longitude=g_cen_lon)
    g_proj["Sinusoidal"] = ccrs.Sinusoidal(central_longitude=g_cen_lon)

    g_proj["InterruptedGoodeHomolosine"] = ccrs.InterruptedGoodeHomolosine(central_longitude=g_cen_lon)

    #g_proj["RotatedPole"] = ccrs.RotatedPole(pole_latitude=g_cen_lat, pole_longitude=g_cen_lon)
    #g_proj["RotatedPole"] = ccrs.RotatedPole()

    g_projs = [ "PlateCarree", "Mercator", "Miller",
                "EquidistantConic", "LambertConformal", "AlbersEqualArea",
                "Stereographic", "Gnomonic", "Orthographic",
                "Robinson", "Mollweide", "Sinusoidal", "InterruptedGoodeHomolosine" ]

    g_projs = [ "PlateCarree", "Orthographic", "LambertConformal", "Gnomonic" ]
    g_projs = [ "LambertConformal", "Gnomonic", "Mercator" ]

    if args.compute:
        g_projs = [ "LambertConformal", "Gnomonic" ]
        g_dim_x = 2; g_dim_y = 1
    else:
        g_projs = [ "LambertConformal" ]
        g_dim_x = 1; g_dim_y = 1

    if mode == "init":
        for p in g_projs:
            g_view[p] = "regional"

def plots_draw(mode):
    global g_proj
    global g_projs
    global g_axis
    global g_view
    global g_cen_lon, g_cen_lat
    global g_dim_x, g_dim_y
    global g_xdata_span, g_ydata_span
    global g_mode
    global g_lambert_extent
    global g_crn_lon, g_crn_lat

    # Save mode
    if not mode == "init":
        for p in g_projs:
            g_mode = g_axis[p].get_navigate_mode()

    # Save new extent info
    if mode == "set":
        g_cen_lon, g_cen_lat, g_crn_lon_new, g_crn_lat_new, g_lambert_extent_new = get_dims('LambertConformal', 'blue')
        g_cen_lat = cen_lat_adjust(g_cen_lat)
        print(f"plots_draw: set: lwr_lon {g_crn_lon_new} lwr_lat {g_crn_lat_new}")

    # Remove any old plots
    if g_projs:
        for p in g_projs:
            g_axis[p].remove()

    # Create/recreate projections (recreate needed when the center changes)
    projs_create(mode)

    j = 1; g_axis = {}
    for p in g_projs:
        g_axis[p] = g_fig.add_subplot(g_dim_x, g_dim_y, j, projection=g_proj[p])
        g_axis[p].margins(x=0.0, y=0.0)
        g_axis[p].add_feature(cfeature.COASTLINE)
        g_axis[p].add_feature(cfeature.BORDERS)
        g_axis[p].add_feature(cfeature.STATES)
        g_axis[p].gridlines()

        # Set axis view (global or regional/zoomed)
        if g_view[p] == "global":
            g_axis[p].set_title(f"{p} global")
            g_axis[p].set_global()
        else:
            g_view[p] = "regional"
            g_axis[p].set_title(f"{p}")

        # Restore mode
        g_axis[p].set_navigate_mode(g_mode)
    
        j = j + 1

    if mode == "init":
        g_crn_lon = g_crn_lon_dflt
        g_crn_lat = g_crn_lat_dflt
        xc, yc = g_proj["LambertConformal"].transform_point(g_cen_lon, g_cen_lat, ccrs.Geodetic())
        xll, yll = g_proj["LambertConformal"].transform_point(g_crn_lon, g_crn_lat, ccrs.Geodetic())
        xlr = xc+(xc-xll); ylr = yc-(yc-yll)
        xul = xll; yul = yc+(yc-yll)
        xur = xlr; yur = yul
        g_lambert_extent = (xll, xlr, yll, yul)
        print(f"plots_draw: init: g_cen_lon {g_cen_lon} g_cen_lat {g_cen_lat} g_crn_lon {g_crn_lon} g_crn_lat {g_crn_lat}")
    elif mode == "center":
        print(f"plots_draw: center: reusing g_lambert_extent = {fmt_tuple(g_lambert_extent)}");
        print(f"plots_draw: center: g_cen_lon {g_cen_lon} g_cen_lat {g_cen_lat} g_crn_lon {g_crn_lon} g_crn_lat {g_crn_lat}")
    elif mode == "set":
        g_crn_lon = g_crn_lon_new
        g_crn_lat = g_crn_lat_new
        xc, yc = g_proj["LambertConformal"].transform_point(g_cen_lon, g_cen_lat, ccrs.Geodetic())
        xll, yll = g_proj["LambertConformal"].transform_point(g_crn_lon, g_crn_lat, ccrs.Geodetic())
        xlr = xc+(xc-xll); ylr = yc-(yc-yll)
        xul = xll; yul = yc+(yc-yll)
        xur = xlr; yur = yul
        g_lambert_extent = (xll, xlr, yll, yul)
        print(f"plots_draw: set: g_cen_lon {g_cen_lon} g_cen_lat {g_cen_lat} g_crn_lon {g_crn_lon} g_crn_lat {g_crn_lat}")

    if True or not args.file:

        # Set LambertConformal extent
        g_axis["LambertConformal"].set_extent(g_lambert_extent, crs=g_proj["LambertConformal"])

        # The write component (shown in blue) is taken from the LambertConformal
        # extents and, as such, will be a perfect rectangle when plotted on the
        # LambertConformal plot. All other projects will show the transformed
        # rectangle.
        for p_tgt in g_projs:
            draw_box_xy_data(p_tgt, "LambertConformal", 'blue') # write component

        # The compute component (shown in red) is taken from the Gnomonic
        # extents and, as such, will be a perfect rectangle when plotted on the
        # Gnomonic plot.  All other projects will show the transformed
        # rectangle.
        if args.compute:
            for p_tgt in g_projs:
                draw_box_xy_data(p_tgt, "Gnomonic", 'red') # compute component

    if mode == "center":
        x1, x2, y1, y2 = g_axis['LambertConformal'].get_extent()
        g_crn_lon, g_crn_lat = ccrs.PlateCarree().transform_point(x1, y1, g_proj['LambertConformal'])
        print(f"plots_draw: center: new lambert corner is lon {g_crn_lon} lat {g_crn_lat}")

    if args.file:
        plot_grib()
        ram = io.BytesIO()
        for p in g_projs:
            g_axis[p].set_title(p)
        g_fig.suptitle(g_title)
        plt.savefig(ram, format="png", bbox_inches="tight", dpi=150)
        ram.seek(0)
        im = PIL.Image.open(ram)
        im2 = im.convert("RGB")
        im2.save(args.file + ".png", format="PNG")
        if args.close:
            exit(0)

    plt.show() # end plots_draw()

def draw_box_xy_data(tgt_index, src_index, color):
    cen_lon, cen_lat, lwr_lon, lwr_lat, extent = get_dims(src_index, color)
    #cen_lon, cen_lat, lwr_lon, lwr_lat, extent = get_dims(tgt_index, color)
    x1, x2, y1, y2 = extent 
    if color == "blue":
        g_axis[tgt_index].plot(*(cen_lon, cen_lat), transform=ccrs.Geodetic(), 
                               marker='*', ms=10, color='red')
    xs = []
    ys = []
    #left
    for i in range(0, 33):
        xs.append(x1)
        ys.append(y1+i*(y2-y1)/32)
    #top
    for i in range(0, 33):
        xs.append(x1+i*(x2-x1)/32)
        ys.append(y2)
    #right
    for i in range(0, 33):
        xs.append(x2)
        ys.append(y2-i*(y2-y1)/32)
    #bottom
    for i in range(0, 33):
        xs.append(x2-i*(x2-x1)/32)
        ys.append(y1)

    g_axis[tgt_index].plot(xs, ys, transform=g_axis[src_index].projection, 
                           color=color, linewidth=2, alpha=1.0, linestyle='solid')

def plot_grib():
    global g_title

    g_data = pygrib.open(args.file)

    zmsg = g_data.select(name="Geopotential height", level=500)[0]
    zdata = zmsg.values
    lats, lons = zmsg.latlons()

    vdata = g_data.select(name="Absolute vorticity", level=500)[0].values

    sio = io.StringIO()
    print(zmsg, file=sio)
    fcst = sio.getvalue().split(':')[6]
    time = sio.getvalue().split(':')[7]
    g_title = fcst + " " + time

    # Filter out by lat/lon (need to generalize this using lambert center and corner)
    if False:
        lat_span = 30
        lon_span = 30
        lat_min, lat_max = round(g_cen_lat-lat_span/2), \
                           round(g_cen_lat+lat_span/2)
        dim0 = (lat_max-lat_min)*4+1
        lon_min, lon_max = round(360+g_cen_lon-lon_span/2), \
                           round(360+g_cen_lon+lon_span/2)
        dim1 = (lon_max-lon_min)*4+1
        mask = (lats >= lat_min) & (lats <= lat_max) & (lons >= lon_min) & (lons <= lon_max)
        filtered_lats = lats[mask]
        filtered_lons = lons[mask]
        filtered_zdata = zdata[mask]
        filtered_vdata = vdata[mask]
        lats = np.array(filtered_lats)
        lons = np.array(filtered_lons)
        zdata = np.array(filtered_zdata)
        vdata = np.array(filtered_vdata)
        lats = np.reshape(lats, (dim0, dim1))
        lons = np.reshape(lons, (dim0, dim1))
        zdata = np.reshape(zdata, (dim0, dim1))
        vdata = np.reshape(vdata, (dim0, dim1))

    # Geopotential height
    z500 = zdata * 0.1
    z500 = scipy.ndimage.gaussian_filter(z500, 6.89)
    for p in g_projs:
        contours = g_axis[p].contour(lons, lats, z500, 
                                     np.arange(0, 900, 6), 
                                     colors="blue", linewidths=1, 
                                     transform=cartopy.crs.PlateCarree(), alpha=0.25)
    plt.clabel(contours, np.arange(0, 900, 6), inline_spacing=1, fmt="%d", fontsize=8)

    # Vorticity
    vort500 = vdata * 100000
    vort500 = scipy.ndimage.gaussian_filter(vort500, 1.7225)
    #vort500[vort500 > 1000] = 0  # Mask out undefined values on domain edge
    vortlevs = [16, 20, 24, 28, 32, 36, 40]
    colorlist = ["yellow", "gold", "goldenrod", "orange", "orangered", "red" ]
    cm = matplotlib.colors.ListedColormap(colorlist)
    norm = matplotlib.colors.BoundaryNorm(vortlevs, cm.N)
    for p in g_projs:
        cs1_a = g_axis[p].pcolormesh(lons, lats,
                                     vort500, transform=cartopy.crs.PlateCarree(),
                                     cmap=cm, norm=norm)
    cs1_a.cmap.set_under("none")   # under 16
    cs1_a.cmap.set_over("darkred") # above 40

def show_help():
    print(g_help)

#
# Main
#

plt.rcParams["figure.raise_window"] = False
g_fig = plt.figure(figsize=(10, 10))
g_fig.canvas.mpl_connect('button_press_event', on_button_press)
g_fig.canvas.mpl_connect('key_press_event', on_key_press)
show_help()
plots_draw("init") #1
