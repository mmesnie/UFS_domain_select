#!/usr/bin/env python3

g_menu = True
g_debug = False

#
# TODO
#
# - limit domain by square nautical miles??
#   we can potentially eliminate any "global" checks, which may catch
#   the same thing
# - fix "-" and "+" when LambertConformal isn't enabled
# ...
# - update Makefile to download latest cycle and integrate into GUI
# - make sure g_compute_grid (gnomonic) is correct and add option to show it
# - write status to screen (e.g., YAML written, res selected, ...)
# - test w/ GFS files (currently filters based on lon/lat spans)
# - understand and integrate regional_latlon
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
parser.add_argument("--file", "-f", help="grib file to plot", required=False)
parser.add_argument("--proj", "-p", help="LambertConformal or RotatedPole (required with -f)", required=False)
parser.add_argument("--close", "-x", help="close after saving plot of grib file", required=False, action="store_true")
args = parser.parse_args()

#
# Internal defaults
#
import os
HOME = f"{os.environ['HOME']}"
UFS_DOMAIN_SELECT_HOME=os.path.dirname(os.path.abspath(__file__))

g_res_dflt=-1                                             # 3000, 13000, 25000, or -1 (auto)
g_yaml_file = f"{UFS_DOMAIN_SELECT_HOME}/build/ufs-srweather-app-v2.2.0/ush/config.yaml" # Location of YAML output
g_compute_grid_dflt = 0.1                                 # 5% larger than write component grid

g_index_dflt = 'InterruptedGoodeHomolosine'
g_index_dflt = 'RotatedPole'
g_index_dflt = 'LambertConformal'
g_index_dflt = 'Orthographic'
g_index_dflt = 'PlateCarree'
g_index_dflt = 'Mercator'

#
# Default YAML values
#
g_dt_atmos = 36
g_blocksize = 40
g_layout_x = 3
g_layout_y = 3
g_write_groups = 1
g_write_tasks_per_group = 3
g_date ='20190615'
g_cycle ='18'
g_fcst_len_hrs = 1
g_lbc_spec_intvl_hrs = 1
g_extrn_mdl_source_basedir_ics = f"{UFS_DOMAIN_SELECT_HOME}/build/DATA-2.2.0/input_model_data/FV3GFS/grib2/{g_date}{g_cycle}"
g_extrn_mdl_source_basedir_lbcs = f"{UFS_DOMAIN_SELECT_HOME}/build/DATA-2.2.0/input_model_data/FV3GFS/grib2/{g_date}{g_cycle}"

#g_cen_lon_dflt=-59.5;   g_cen_lat_dflt=-51.7; g_crn_lon_dflt=-61.98;  g_crn_lat_dflt=-52.81 # Falkland Islands
g_cen_lon_dflt=-127.68; g_cen_lat_dflt=45.72; g_crn_lon_dflt=-132.86; g_crn_lat_dflt=41.77  # Oregon coast
#g_cen_lon_dflt=-97.5;   g_cen_lat_dflt=38.5;  g_crn_lon_dflt=-122.72; g_crn_lat_dflt=21.14  # CONUS
#g_cen_lon_dflt=-61.13;   g_cen_lat_dflt=10.65;  g_crn_lon_dflt=-61.98; g_crn_lat_dflt=9.85  # CONUS
#g_cen_lon_dflt=-141.87; g_cen_lat_dflt=40.48; g_crn_lon_dflt=-160.29; g_crn_lat_dflt=16.64  # Eastern Pacific 

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
from matplotlib.widgets import RadioButtons
from matplotlib.widgets import CheckButtons

#
# Function definitions
#
def init_dflts():
    global g_cen_lon, g_cen_lat
    global g_crn_lon, g_crn_lat
    global g_compute_grid
    global g_res
    global g_index

    if args.cen_lon:
        g_cen_lon = float(args.cen_lon)
    else:
        g_cen_lon = g_cen_lon_dflt
    if args.cen_lat:
        g_cen_lat = float(args.cen_lat)
    else:
        g_cen_lat = g_cen_lat_dflt
    if args.crn_lon:
        g_crn_lon = float(args.crn_lon)
    else:
        g_crn_lon = g_crn_lon_dflt
    if args.crn_lat:
        g_crn_lat = float(args.crn_lat)
    else:
        g_crn_lat = g_crn_lat_dflt
    if args.proj:
        g_index = args.proj
    else:
        g_index = g_index_dflt

    g_compute_grid = g_compute_grid_dflt
    g_res = g_res_dflt

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
    global g_crn_lon, g_crn_lat

    if not (event.inaxes and event.button is MouseButton.MIDDLE):
        return

    index = get_index(event.inaxes)

    print(f"centering projection ({index})")

    if g_view[index] == "regional" and g_global[index] == g_axis[index].get_extent():
        print(f"EXTENT IS REGIONAL BUT GLOBAL: {g_axis[index].get_extent()}")
        #g_axis[index].set_global()
        force_global = True
    else:
        force_global = False

    if g_view[index] == "global":
        try:
            g_axis[index].set_extent(g_extent[index], crs=g_proj[index])
            print(f"on_button_press: A: set extent for {index}: {g_extent[index]}")
        except:
            print(f"on_button_press: A: failed to set extent for {index}: {g_extent[index]} -- setting global")
            g_axis[index].set_global()
        restore = True
    else:
        restore = False

    g_cen_lon, g_cen_lat = ccrs.PlateCarree().transform_point(event.xdata, event.ydata,
                                                              event.inaxes.projection)
    x1, x2, y1, y2 = g_axis[index].get_extent()
    if no_cen_lat(index):
        yc = event.ydata
        _, _, lower, upper = g_global[index]
        if yc-(y2-y1)/2>=lower and yc+(y2-y1)/2<=upper:
            g_extent[index] = -(x2-x1)/2, +(x2-x1)/2, yc-(y2-y1)/2, yc+(y2-y1)/2
        else:
            if yc-(y2-y1)/2<lower:
                print(f"HIT BOTTOM EDGE: {index}")
            if yc+(y2-y1)/2>upper:
                print(f"HIT TOP EDGE: {index}")
            print(f"g_global{index} = {g_global[index]}")
            g_extent[index] = -(x2-x1)/2, +(x2-x1)/2, y1, y2
    else:
        g_extent[index] = -(x2-x1)/2, +(x2-x1)/2, -(y2-y1)/2, +(y2-y1)/2

    g_axis[index].remove()
    g_proj[index] = proj_create(index)
    g_axis[index] = g_fig.add_subplot(g_dim_x, g_dim_y, g_loc[index], 
                                      projection=g_proj[index])
    g_axis[index].margins(x=0.0, y=0.0)
    g_axis[index].add_feature(cfeature.COASTLINE)
    g_axis[index].add_feature(cfeature.BORDERS)
    g_axis[index].add_feature(cfeature.STATES)
    g_axis[index].gridlines()
    g_axis[index].set_title(index + " (centered)")
    g_axis[index].plot(*(g_cen_lon, g_cen_lat), transform=ccrs.Geodetic(), 
                         marker='*', ms=20, color='green')
    if force_global:
        g_axis[index].set_global()
        g_axis[index].set_title(index + " (centered forced global)")
    else:
        try:
            g_axis[index].set_extent(g_extent[index], crs=g_proj[index])
            #print(f"on_button_press: B: set extent for {index}: {g_extent[index]}")
        except:
            print(f"on_button_press: B: failed to set extent for {index}: {g_extent[index]} -- setting global")
            g_axis[index].set_global()

    #g_axis[index].callbacks.connect('xlim_changed', on_xlim_changed)
    #g_axis[index].callbacks.connect('ylim_changed', on_ylim_changed)

    if restore:
        g_axis[index].set_global()
        g_axis[index].set_title(index + " (global center restored)")
        debug(f"restore: global extent is {g_axis[index].get_extent()}")

    plt.draw()

    # Now overlay the actual center
    x1, x2, y1, y2 = g_axis[index].get_extent()
    cen_lon, cen_lat = ccrs.PlateCarree().transform_point((x1+x2)/2, (y1+y2)/2,
                                                           g_axis[index].projection)
    g_axis[index].plot(*(cen_lon, cen_lat), transform=ccrs.Geodetic(), 
                         marker='*', ms=10, color='red')

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
    res = g_res

    # computational grid
    cen_lon_r, cen_lat_r, lwr_lon_r, lwr_lat_r, extent_r = get_dims(index, "red")
    x1, x2, y1, y2 = extent_r
    xspan_r = abs(x2-x1)
    yspan_r = abs(y2-y1)
    if g_index == "RotatedPole":
        max_delta = pick_max_delta(xspan_r*1852*60, yspan_r*1852*60)
    else:
        max_delta = pick_max_delta(xspan_r, yspan_r)
    if g_res == -1:
        res = max_delta
        print(f"             resolution: automatically set to {res} meters")
    elif (g_res>max_delta):
        res = max_delta
        print(f"             resolution: forced to max ({res} meters)")
    percent = round(100*(1+g_compute_grid))
    print(f"      compute grid size: {percent}% of write component")
    if g_index == "RotatedPole":
        nx_r = int(xspan_r*1852*60/res)
        ny_r = int(yspan_r*1852*60/res)
    else:
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

    config_text = f"""
metadata:
  description: >-
    Automatically generated via UFS_domain_select.py
user:
  RUN_ENVIR: community
  MACHINE: linux
  ACCOUNT: an_account
workflow:
  EXPT_SUBDIR: {UFS_DOMAIN_SELECT_HOME}/build/expt-2.2.0/test_community
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
task_run_post:
  POST_OUTPUT_DOMAIN_NAME: 'mesnier'
task_run_fcst:
  QUILTING: true
  DT_ATMOS: {g_dt_atmos}
  LAYOUT_X: {g_layout_x} 
  LAYOUT_Y: {g_layout_y}
  BLOCKSIZE: {g_blocksize}
  WRTCMP_write_groups: {g_write_groups}
  WRTCMP_write_tasks_per_group: {g_write_tasks_per_group}"""

    if g_index == "LambertConformal":
        config_text_wrtcmp = f"""
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
"""

    elif g_index == "RotatedPole":
        config_text_wrtcmp = f"""
  WRTCMP_output_grid: "rotated_latlon"
  WRTCMP_cen_lon: {g_cen_lon}
  WRTCMP_cen_lat: {g_cen_lat} 
  WRTCMP_lon_lwr_left: {x1} 
  WRTCMP_lat_lwr_left: {y1} 
  WRTCMP_lon_upr_rght: {x2} 
  WRTCMP_lat_upr_rght: {y2}
  WRTCMP_dlon: {res/(1852*60)} 
  WRTCMP_dlat: {res/(1852*60)} 
"""

    print(config_text+config_text_wrtcmp)

    with open(g_yaml_file, "w") as file:
        file.write(config_text)
        file.write(config_text_wrtcmp)
    print(f"       YAML file output: {g_yaml_file}")
    print(f"                    cmd: export DATE={g_date} CYCLE={g_cycle} LEN={g_fcst_len_hrs} LBC={g_lbc_spec_intvl_hrs}; time ./forecast")

def on_key_press(event):
    global g_compute_grid
    global g_res
    global g_extent
    global g_index

    g_index = g_index_dflt
    if event.inaxes:
        g_index = get_index(event.inaxes)

    if event.key == 'h':
        show_help()
    elif event.key == '0':
        print(f"restoring default settings")
        plots_draw("init")
    elif event.key == 'y':
        if (g_index == "LambertConformal" or g_index == "RotatedPole"):
            output_config(g_index)
        else:
            print("choose from LambertConformal or RotatedPole")
            return
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
        plots_draw("set")
    elif event.key == '-':
        g_compute_grid -= .01
        percent = round(100*(1+g_compute_grid))
        print(f"set compute grid to {percent}% of write component")
        plots_draw("set")
    elif event.key == ' ':
        print(f"setting source projection ({g_index})")
        if (g_axis[g_index].get_extent() == g_global[g_index]):
            print(f"SORRY - ATTEMPT TO SET WITH GLOBAL EXTENT: {g_global[g_index]}")
            return
        plots_draw("set")
    elif event.key == 'g':
        print(f"toggling global view (source {g_index})")
        if g_view[g_index] == "regional":
            g_extent[g_index] = g_axis[g_index].get_extent()
            g_view[g_index] = "global"
            g_axis[g_index].set_global()
            g_axis[g_index].set_title(g_index + " (global toggle)")
        elif g_view[g_index] == "global": 
            g_view[g_index] = "regional"
            g_axis[g_index].set_extent(g_extent[g_index], crs=g_proj[g_index])
            g_axis[g_index].set_title(g_index + " (regional toggle)")

    elif event.key == 'q':
        exit(0)

def get_index(ax):
    for p in g_projs:
        if g_axis[p] == ax:
            return p

def proj_create(index):
    match index:

        # Cylindrical
        case "PlateCarree":
            return ccrs.PlateCarree(central_longitude=g_cen_lon)
        case "Mercator":
            return ccrs.Mercator(central_longitude=g_cen_lon)
        case "Miller":
            return ccrs.Miller(central_longitude=g_cen_lon)

        # Pseudo Cylindrical
        case "Sinusoidal":
            return ccrs.Sinusoidal(central_longitude=g_cen_lon)
        case "Mollweide":
            return ccrs.Mollweide(central_longitude=g_cen_lon)
        case "Robinson":
            return ccrs.Robinson(central_longitude=g_cen_lon)
        case "InterruptedGoodeHomolosine":
            return ccrs.InterruptedGoodeHomolosine(central_longitude=g_cen_lon)

        case "EquidistantConic":
            return ccrs.EquidistantConic(central_longitude=g_cen_lon, central_latitude=g_cen_lat)
        case "AlbersEqualArea":
            return ccrs.AlbersEqualArea(central_longitude=g_cen_lon, central_latitude=g_cen_lat)
        case "LambertConformal":
            if g_cen_lat>0:
                lat_cutoff = -30
            else:
                lat_cutoff = 30
            return ccrs.LambertConformal(central_longitude=g_cen_lon, 
                                         central_latitude=g_cen_lat,
                                         standard_parallels=(g_cen_lat, g_cen_lat), 
                                         cutoff=lat_cutoff)
        case "Stereographic":
            return ccrs.Stereographic(central_longitude=g_cen_lon, central_latitude=g_cen_lat)
        case "Orthographic":
            return ccrs.Orthographic(central_longitude=g_cen_lon, central_latitude=g_cen_lat)
        case "Gnomonic":
            return ccrs.Gnomonic(central_longitude=g_cen_lon, central_latitude=g_cen_lat)
        case "RotatedPole":
            return ccrs.RotatedPole(pole_latitude=90-g_cen_lat, pole_longitude=g_cen_lon-180)

def no_cen_lat(index):
    if index == "Mercator" or index == "Miller" or index == "PlateCarree" \
            or index == "Robinson" or index == "Sinusoidal" \
            or index == "Mollweide" or index == "InterruptedGoodeHomolosine":
            return True
    else:
            return False

def projs_create(mode):
    global g_proj
    global g_projs
    global g_plotted
    global g_view
    global g_extent
    global g_cen_lat
    global g_dim_x, g_dim_y
    global g_color
    global g_enabled
    global g_global

    if mode == "init":
        g_extent = {}
        g_view = {}
        g_global = {}

    g_proj = {}
    g_plotted = {}
    g_color = {}

    g_cen_lat = cen_lat_adjust(g_cen_lat)

    # Class: Cylindrical, Conic, Planar/Azimuthal
    # Case: Tangent, Secant
    # Aspect: Normal, Transverse, Oblique

    # Cylindrical
    g_proj["PlateCarree"] = proj_create("PlateCarree")
    g_proj["Mercator"] = proj_create("Mercator")
    g_proj["Miller"] = proj_create("Miller")

    # Conic
    g_proj["EquidistantConic"] = proj_create("EquidistantConic") 
    g_proj["AlbersEqualArea"] = proj_create("AlbersEqualArea")
    g_proj["LambertConformal"] = proj_create("LambertConformal")

    # Planar/Azimuthal
    g_proj["Stereographic"] = proj_create("Stereographic")
    g_proj["Orthographic"] = proj_create("Orthographic")
    g_proj["Gnomonic"] = proj_create("Gnomonic")

    # Pseudo-cylindrical
    g_proj["Robinson"] = proj_create("Robinson")
    g_proj["Mollweide"] = proj_create("Mollweide")
    g_proj["Sinusoidal"] = proj_create("Sinusoidal")
    g_proj["InterruptedGoodeHomolosine"] = proj_create("InterruptedGoodeHomolosine")

    # Other
    g_proj["RotatedPole"] = proj_create("RotatedPole")

    if mode == "init":
        g_enabled = {}
        for p in g_proj:
            g_enabled[p] = False
        if args.file:
            g_enabled[g_index_dflt] = True
        else:
            g_enabled['Mercator'] = True
            #g_enabled['PlateCarree'] = True
            #g_enabled['Miller'] = True
            #g_enabled['LambertConformal'] = True
            #g_enabled['RotatedPole'] = True
            #g_enabled['InterruptedGoodeHomolosine'] = True
            #g_enabled['Gnomonic'] = True
            #g_enabled['Orthographic'] = True

    g_projs = []

    count = 0
    for p in g_proj:
        if g_enabled[p]:
            g_projs.append(p)
            count += 1

    match count:
        case 1:
            g_dim_x = 1; g_dim_y = 1
        case 2:
            g_dim_x = 1; g_dim_y = 2
        case 3 | 4:
            g_dim_x = 2; g_dim_y = 2
        case 5 | 6:
            g_dim_x = 2; g_dim_y = 3
        case 7 | 8:
            g_dim_x = 3; g_dim_y = 3
        case 9:
            g_dim_x = 3; g_dim_y = 3
        case 10 | 11 | 12:
            g_dim_x = 4; g_dim_y = 3
        case 13 | 14:
            g_dim_x = 4; g_dim_y = 4

    # Give each projection a color (brown is the default).
    for p in g_projs:
        g_color[p] = "brown"
        g_plotted[p] = False
    g_color["LambertConformal"] = "blue"
    g_color["RotatedPole"] = "purple"
    g_color["Mercator"] = "yellow"
    g_color["Gnomonic"] = "green"
    g_color["Orthographic"] = "orange"

    if mode == "init":
        for p in g_proj:
            g_view[p] = "regional"

def plots_remove():
    global g_plotted

    if (g_menu):
        g_axis['menu1'].remove()
        g_axis['menu2'].remove()
    if g_projs:
        for p in g_projs:
            if g_plotted[p]:
                g_axis[p].remove()
                g_plotted[p] = False

def find_extent(tx, ty):
    min_x = tx[0]; max_x = tx[0]
    min_y = ty[0]; max_y = ty[0]
    for x in tx:
        if x < min_x: min_x = x
        if x > max_x: max_x = x
    for y in ty:
        if y < min_y: min_y = y
        if y > max_y: max_y = y
    return (min_x, max_x, min_y, max_y)

def plots_draw(mode):
    global g_proj
    global g_projs
    global g_axis
    global g_view
    global g_cen_lon, g_cen_lat
    global g_crn_lon, g_crn_lat
    global g_dim_x, g_dim_y
    global g_xdata_span, g_ydata_span
    global g_compute_grid
    global g_res
    global g_index
    global g_enabled
    global g_check1
    global g_check2
    global g_loc

    #
    # Notes on global variables:
    #
    # g_index is the index of the source axis that we're trransforming from.
    # g_projs is the list of projections, only one of which is g_index.
    #

    # Change all target axes to regional view. Any target
    # axis in global view will be returned to global view
    # at the end of this procedure.
    if not mode == "init":
        restore_global = {}
        for p in g_projs:
            if p == g_index:
                g_extent[p] = g_axis[p].get_extent()
            elif (not p == g_index) and g_plotted[p] and g_enabled[p]:
                try:
                    g_axis[p].set_extent(g_extent[p], crs=g_proj[p])
                except:
                    print(f"*** WARN: FAILED TO SET EXTENT ({p}) ***")
                if g_view[p] == "global":
                    restore_global[p] = True
                g_view[p] = "regional"

    if mode == "init":
        init_dflts()
    else:
        if mode == "set":
            g_cen_lon, g_cen_lat, g_crn_lon, g_crn_lat, (x1, x2, y1, y2) = get_dims(g_index, 'blue')
        elif mode == "center":
            _, _, _, _, (x1, x2, y1, y2) = get_dims(g_index, 'blue')

        if no_cen_lat(g_index):
            xc, yc = g_proj[g_index].transform_point(g_cen_lon, g_cen_lat, ccrs.Geodetic())
            # This calculation doesn't work for InterruptedGoodeHomolosine
            if g_index == "InterruptedGoodeHomolosine":
                print(f"SPECIAL CASE FOR InterruptedGoodeHomolosine (set)")
                new_extent = (xc-abs(x2-x1)/2, xc+abs(x2-x1)/2, yc-abs(y2-y1)/2, yc+abs(y2-y1)/2)
            else:
                new_extent = (-abs(x2-x1)/2, abs(x2-x1)/2, yc-abs(y2-y1)/2, yc+abs(y2-y1)/2)
        else:
            new_extent = (-abs(x2-x1)/2, abs(x2-x1)/2,   -abs(y2-y1)/2,    abs(y2-y1)/2)

    # Remove old plots (the try will fail on the first "init", as g_projs hasn't
    # yet been defined via projs_create()
    try:
        plots_remove()
    except:
        debug("NO PLOTS TO REMOVE")

    # Create projections
    projs_create(mode)

    # Create plots
    g_axis = {}
    g_loc = {}
    j = 1

    debug(f"CREATING PLOTS w/ INDEX {g_index}")

    for p in g_projs:

        g_axis[p] = g_fig.add_subplot(g_dim_x, g_dim_y, j, projection=g_proj[p])
        g_axis[p].margins(x=0.0, y=0.0)
        g_axis[p].add_feature(cfeature.COASTLINE)
        g_axis[p].add_feature(cfeature.BORDERS)
        g_axis[p].add_feature(cfeature.STATES)
        g_axis[p].gridlines()
        #g_axis[p].callbacks.connect('xlim_changed', on_xlim_changed)
        #g_axis[p].callbacks.connect('ylim_changed', on_ylim_changed)
        g_loc[p] = j

        # Set axis view (global or regional/zoomed)
        debug(f"p is {p}, mode is {mode}, g_view[{p}] is {g_view[p]}")

        if ((p != g_index) and g_view[p] == "global"):
            g_axis[p].set_title(f"{p} ({mode} global)")
            g_axis[p].set_global()

        else:
            g_view[p] = "regional"
            if p == g_index:
                g_axis[p].set_title(f"{p} (regional index)")
            else:
                g_axis[p].set_title(f"{p} (regional non-index)")

        g_plotted[p] = True

        j = j + 1

    if mode == "init":
        xc, yc = g_proj[g_index].transform_point(g_cen_lon, g_cen_lat, ccrs.Geodetic())
        xll, yll = g_proj[g_index].transform_point(g_crn_lon, g_crn_lat, ccrs.Geodetic())
        xlr = xc+(xc-xll); ylr = yc-(yc-yll)
        xul = xll; yul = yc+(yc-yll)
        xur = xlr; yur = yul
        g_extent[g_index] = (xll, xlr, yll, yul)
    else:
        g_extent[g_index] = new_extent
        debug("CASE 2")

    #print(f"A: set g_extent to {g_extent[g_index]}")

    if args.file:
        print(f"PLOTTING GRIB FILE {args.file}")
        plot_grib()
        print(f"DONE PLOTTING GRIB FILE {args.file}")
        for p in g_projs:
            if g_plotted[p]:
                g_axis[p].set_title(p)
        g_fig.suptitle(g_title)
        if False:
            ram = io.BytesIO()
            print(f"SAVING GRIB PLOT TO {args.file}")
            plt.savefig(ram, format="png", bbox_inches="tight", dpi=150)
            ram.seek(0)
            im = PIL.Image.open(ram)
            im2 = im.convert("RGB")
            im2.save(args.file + ".png", format="PNG")
            print(f"DONE SAVING GRIB PLOT TO {args.file}")
        if args.close:
            exit(0)

    if True:

        assert(g_view[g_index] == "regional")
        assert(mode=="set" or mode=="init")

        try:
            if not args.file:
                print(f"GRIB FILE NOT SPECIFIED: setting default extent for {g_index}")
                g_axis[g_index].set_extent(g_extent[g_index], crs=g_proj[g_index])
                debug(f"A:   set extent to {fmt_tuple(g_extent[g_index])}")
            else:
                print(f"GRIB FILE SPECIFIED: not setting extent for {g_index}")
            g_axis[g_index].set_title(f"{g_index} (STILL REGIONAL after mode {mode})")
        except:
            print(f"*** CASE 2: FAILED TO SET EXTENT ({g_index}): {g_extent[g_index]} ***")

        # Outline the extent of the index's axis
        x, y = create_box_xy_data(g_index, g_color[g_index])
        g_axis[g_index].plot(x, y, color=g_color[g_index], linewidth=2, alpha=1.0, linestyle='dashed')
        g_axis[g_index].plot(*(g_cen_lon, g_cen_lat), transform=ccrs.Geodetic(), 
                               marker='*', ms=10, color='purple')

        targets = [t for t in g_projs if not t == g_index]
        for p in targets:

            style = 'solid' 

            # Draw the transformed index's extent on the target's axis
            pts = g_proj[p].transform_points(g_proj[g_index], np.array(x), np.array(y))
            tx = pts[:, 0]; ty = pts[:, 1]
            g_axis[p].plot(tx, ty, color=g_color[g_index], linewidth=2, alpha=1.0, linestyle=style)

            # Outline the extent of the target's axis
            min_x, max_x, min_y, max_y = find_extent(tx, ty)
            xp, yp = create_box_xy((min_x, max_x, min_y, max_y))
            g_axis[p].plot(xp, yp, color=g_color[p], linewidth=2, alpha=1.0, linestyle=style)

            # Draw the transformed targets's extent on the index's axis
            pts = g_proj[g_index].transform_points(g_proj[p], np.array(xp), np.array(yp))
            tx = pts[:, 0]; ty = pts[:, 1]
            g_axis[g_index].plot(tx, ty, color=g_color[p], linewidth=2, alpha=1.0, linestyle=style)

    if mode == "center":
        x1, x2, y1, y2 = g_axis[g_index].get_extent()
        g_crn_lon, g_crn_lat = ccrs.PlateCarree().transform_point(x1, y1, g_proj[g_index])

    # Set/restore view
    if mode == "init":
        for p in g_projs:
            if g_plotted[p]:
                if p == "Orthographic":
                    print(f"special case: initializing global view for Orthographic")
                    g_extent[p] = g_axis[p].get_extent()
                    g_view[p] = "global"
                    g_axis[p].set_global()
                    g_axis[p].set_title(p + " (global init)")
    else:
        for p in restore_global:
            if g_enabled[p]:
                debug(f"restoring global view for {p} ***")
                g_extent[p] =  g_axis[p].get_extent()
                g_axis[p].set_global()
                g_view[p] = "global"
                g_axis[p].set_title(p + " (global restore)")

    # Keep track of what the global extent is
    for p in g_projs:
            debug(f"recording global extent for {p} (current view {g_view[p]})")
            if g_view[p] == "global":
                g_global[p] = g_axis[p].get_extent()
            else:
                g_extent[p] = g_axis[p].get_extent()
                g_axis[p].set_global()
                g_global[p] = g_axis[p].get_extent()
                try:
                    g_axis[p].set_extent(g_extent[p], crs=g_proj[p])
                except:
                    print(f"FAILED TO SET EXTENT FOR {p}: {g_extent[p]}")

    # Check buttons
    if (g_menu):
        g_axis['menu1'] = g_fig.add_axes([0.25, 0.0, 0.1, 0.1], frameon=False)
        g_axis['menu2'] = g_fig.add_axes([0.65, 0.0, 0.1, 0.1], frameon=False)

    proj1 = []
    proj2 = []
    status1 = []
    status2 = []
    count = 1
    for p in g_proj:
        if count<len(g_proj)/2:
            if g_enabled[p]:
                status1.append(True)
            else:
                status1.append(False)
            proj1.append(p)
        else:
            if g_enabled[p]:
                status2.append(True)
            else:
                status2.append(False)
            proj2.append(p)
        count += 1
    if (g_menu):
        g_check1 = CheckButtons(g_axis['menu1'], proj1, status1)
        g_check1.on_clicked(checkfunc)
        g_check2 = CheckButtons(g_axis['menu2'], proj2, status2)
        g_check2.on_clicked(checkfunc)

    plt.show()

def checkfunc(label):
    global g_enabled
    global g_index
    global g_view

    if not g_enabled[label]:
        g_enabled[label] =  True
    else:
        if len(g_projs) == 1:
            debug(f"checkfunc: cannot remove last projection")
            plots_draw("set")
            return
        elif label == g_index:
            g_enabled[label] =  False
            g_index = g_projs[0]
            debug(f"checkfunc: changed index to {g_index}")
            if label == g_index:
                g_index = g_projs[1]
                debug(f"checkfunc: changed index again to {g_index}")
            g_view[g_index] = "regional"
            g_axis[g_index].set_extent(g_extent[g_index], crs=g_proj[g_index])
            debug(f"checkfunc: forcing {g_index} to regional")
        else:
            debug(f"checkfunc: forcing {g_index} to regional")
            g_view[g_index] = "regional"
            g_axis[g_index].set_extent(g_extent[g_index], crs=g_proj[g_index])
        g_enabled[label] = False

    plots_draw("set")

def create_box_xy(extent):
    x1, x2, y1, y2 = extent 
    xs = []
    ys = []
    steps = 64
    #left
    for i in range(0, steps+1):
        xs.append(x1)
        ys.append(y1+i*(y2-y1)/steps)
    #top
    for i in range(0, steps+1):
        xs.append(x1+i*(x2-x1)/steps)
        ys.append(y2)
    #right
    for i in range(0, steps+1):
        xs.append(x2)
        ys.append(y2-i*(y2-y1)/steps)
    #bottom
    for i in range(0, steps+1):
        xs.append(x2-i*(x2-x1)/steps)
        ys.append(y1)

    return xs, ys

def create_box_xy_data(src_index, color):
    cen_lon, cen_lat, lwr_lon, lwr_lat, extent = get_dims(src_index, color)
    x1, x2, y1, y2 = extent 
    xs = []
    ys = []

    #print(f"cen_lon {cen_lon} cen_lat {cen_lat} lwr_lon {lwr_lon} lwr_lat {lwr_lat}")

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

    return xs, ys

def draw_box_xy_data(tgt_index, src_index, color):
    cen_lon, cen_lat, lwr_lon, lwr_lat, extent = get_dims(src_index, color)
    x1, x2, y1, y2 = extent 
    if color == "blue":
        g_axis[tgt_index].plot(*(cen_lon, cen_lat), transform=ccrs.Geodetic(), 
                               marker='*', ms=10, color='orange')
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

    return xs, ys

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
        if g_plotted[p]:
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
        if g_plotted[p]:
            cs1_a = g_axis[p].pcolormesh(lons, lats,
                                         vort500, transform=cartopy.crs.PlateCarree(),
                                         cmap=cm, norm=norm)
    cs1_a.cmap.set_under("none")   # under 16
    cs1_a.cmap.set_over("darkred") # above 40

def show_help():
    print(g_help)

#def on_resize(event):
#    print(f"Figure resized: width={event.width}, height={event.height}")

#def on_xlim_changed(axes):
#    print("X-axis limits changed:", axes.get_xlim())

#def on_ylim_changed(axes):
#    print("Y-axis limits changed:", axes.get_ylim())

def debug(format):
    if g_debug:
        print(format)

#
# Main
#

#if args.file and not (args.cen_lon and args.cen_lat and args.proj):
#    print("--file requires --cen_lon, --cen_lat, and --proj")
#    exit(1)
#else:

if args.file:
    g_data = pygrib.open(args.file)
    msg = g_data[1]
    debug(f"params: {msg.projparams}")
    match msg.projparams['proj']:
        case "lcc":
            g_index_dflt = 'LambertConformal'
            g_cen_lon_dflt = msg.projparams['lon_0']
            g_cen_lat_dflt = msg.projparams['lat_0']
            print(f"lon_0: {msg.projparams['lon_0']}")
            print(f"lat_0: {msg.projparams['lat_0']}")
        case "ob_tran":
            g_index_dflt = 'RotatedPole'
            g_cen_lon_dflt = msg.projparams['lon_0']
            g_cen_lat_dflt = 90 - msg.projparams['o_lat_p']
        case _:
            print(f"unsupported projection {msg.projparams['proj']}")
            exit(1)
    print(f" proj: {msg.projparams['proj']} --> {g_index_dflt}")

plt.rcParams["figure.raise_window"] = False
g_fig = plt.figure(figsize=(10, 10))
g_fig.canvas.mpl_connect('button_press_event', on_button_press)
g_fig.canvas.mpl_connect('key_press_event', on_key_press)
#g_fig.canvas.mpl_connect('resize_event', on_resize)
show_help()
plots_draw("init")
