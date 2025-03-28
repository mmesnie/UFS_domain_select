#!/usr/bin/env python3

#
# TODO / FIXME
#
# - double star
# - messed up blue/red when we go global
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
*        Key R: set resolution (3km, 13km, 25km)                             *
*        Key 0: reset all settings                                           *
*        Key h: show this help                                               *
* Middle Click: center grid                                                  *
******************************************************************************
"""

#
# The second use is to graph a grib file (500 MB geopotential and vorticity), using
# the -f <FILENAME> option. This can be a GFS grib file (initial or boundary conditions)
# or the output from a SRW forecast (i.e., the grib files in "postprd" directory). The
# center latitude and longitude can be specified via the --cen_lon and _cen_lat options,
# and the longitude and latitude extents via --lon_span and --lat_span. Run with -h for
# more details.
#
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--cen_lon", help="center longitude", required=False)
parser.add_argument("--cen_lat", help="center latitude", required=False)
parser.add_argument("--lon_span", help="longitude span (degrees)", required=False)
parser.add_argument("--lat_span", help="latitude span (degrees)", required=False)
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

#
# User globals (modify as needed)
#

g_cen_lon_default = -97.5; g_cen_lat_default =  38.5; g_lon_span_default = 60; g_lat_span_default = 30; g_res=25000  # CONUS
g_cen_lon_default = -59.5; g_cen_lat_default = -51.8; g_lon_span_default = 10; g_lat_span_default =  5; g_res=13000  # Falkland Islands
g_cen_lon_default = -143;  g_cen_lat_default =  36;   g_lon_span_default = 42; g_lat_span_default = 34; g_res=25000  # Eastern Pacific

g_halo_width = 6
g_dt_atmos = 36
g_blocksize = 40
g_layout_x = 3
g_layout_y = 3
g_write_groups = 1
g_write_tasks_per_group = 3

#
# Internal globals (don't modify)
#

g_compute_grid_default = 0.05
if args.cen_lon:
    g_cen_lon_default = float(args.cen_lon)
if args.cen_lat:
    g_cen_lat_default = float(args.cen_lat)
if args.lon_span:
    g_lon_span_default = float(args.lon_span)
if args.lat_span:
    g_lat_span_default = float(args.lat_span)
g_default_view = "regional" # regional or global
g_cen_lon = g_cen_lon_default
g_cen_lat = g_cen_lat_default
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
        print(f"centering plots at mouse click...")
        g_cen_lon, g_cen_lat = ccrs.PlateCarree().transform_point(event.xdata, event.ydata, 
                                                                  event.inaxes.projection)
        g_cen_lat = cen_lat_adjust(g_cen_lat)
        plots_draw("centering")

def fmt_dims(dims):
    return tuple([str(round(x,2)) if isinstance(x, float) else x for x in dims])

def get_dims(index, color):
    extent = g_axis[index].get_extent()
    x1, x2, y1, y2 = extent

    if color == "red":
        extent = g_axis[index].get_extent()
        xspan = abs(x2-x1)
        x1 -= xspan*(g_compute_grid)/2
        x2 += xspan*(g_compute_grid)/2
        yspan = abs(y2-y1)
        y1 -= yspan*(g_compute_grid)/2
        y2 += yspan*(g_compute_grid)/2
        extent = x1, x2, y1, y2

    nxc = (x1+x2)/2
    nyc = (y1+y2)/2
    xspan = abs(x2-x1)
    yspan = abs(y2-y1)
    cen_lon, cen_lat = ccrs.PlateCarree().transform_point(nxc, nyc, g_axis[index].projection)
    lwr_lon, lwr_lat = ccrs.PlateCarree().transform_point(x1, y1, g_axis[index].projection)
    upr_lon, upr_lat = ccrs.PlateCarree().transform_point(x2, y2, g_axis[index].projection)

    lft_lon, lft_lat = ccrs.PlateCarree().transform_point(x1, nyc, g_axis[index].projection)
    rgt_lon, rgt_lat = ccrs.PlateCarree().transform_point(x2, nyc, g_axis[index].projection)
    top_lon, top_lat = ccrs.PlateCarree().transform_point(nxc, y2, g_axis[index].projection)
    bot_lon, bot_lat = ccrs.PlateCarree().transform_point(nxc, y1, g_axis[index].projection)
    lon_span = abs(lft_lon-rgt_lon)
    lat_span = abs(top_lat-bot_lat)

    return cen_lon, cen_lat, lwr_lon, lwr_lat, upr_lon, upr_lat, xspan, yspan, extent, lon_span, lat_span

def output_yaml(index):
    if index != "LambertConformal":
        print("Write component must be written from LambertConformal grid")
        return
    # computational grid
    cen_lon_r, cen_lat_r, lwr_lon_r, lwr_lat_r, \
        upr_lon_r, upr_lat_r, xspan_r, yspan_r, extent_r, _, _ = get_dims(index, "red")
    # write component grid
    cen_lon_b, cen_lat_b, lwr_lon_b, lwr_lat_b, \
        upr_lon_b, upr_lat_b, xspan_b, yspan_b, extent_b, _, _ = get_dims(index, "blue")
    for res in [ g_res ]:
        nx_r = int(xspan_r/res)
        ny_r = int(yspan_r/res)
        #print(f"\"CUSTOM_XX_{int(res/1000)}km\":")
        print(f"\"CUSTOM_XX\":")
        print(f"  GRID_GEN_METHOD: \"ESGgrid\"")
        print(f"  ESGgrid_LON_CTR: {round(cen_lon_r, 2)}")
        print(f"  ESGgrid_LAT_CTR: {round(cen_lat_r, 2)}")
        print(f"  ESGgrid_NX: {nx_r}")
        print(f"  ESGgrid_NY: {ny_r}")
        print(f"  ESGgrid_DELX: {res}")
        print(f"  ESGgrid_DELY: {res}")
        print(f"  ESGgrid_PAZI: 0.0")
        print(f"  ESGgrid_WIDE_HALO_WIDTH: {g_halo_width}")
        print(f"  DT_ATMOS: {g_dt_atmos}")
        print(f"  LAYOUT_X: {g_layout_x}")
        print(f"  LAYOUT_Y: {g_layout_y}")
        print(f"  BLOCKSIZE: {g_blocksize}")
        nx_b = int(xspan_b/res)
        ny_b = int(yspan_b/res)
        print(f"  QUILTING:")
        print(f"    WRTCMP_write_groups: {g_write_groups}")
        print(f"    WRTCMP_write_tasks_per_group: {g_write_tasks_per_group}")
        print(f"    WRTCMP_output_grid: \"lambert_conformal\"")
        print(f"    WRTCMP_cen_lon: {round(cen_lon_b, 2)}")
        print(f"    WRTCMP_cen_lat: {round(cen_lat_b, 2)}")
        print(f"    WRTCMP_stdlat1: {round(cen_lat_b, 2)}")
        print(f"    WRTCMP_stdlat2: {round(cen_lat_b, 2)}")
        print(f"    WRTCMP_nx: {nx_b}")
        print(f"    WRTCMP_ny: {ny_b}")
        print(f"    WRTCMP_lon_lwr_left: {round(lwr_lon_b, 2)}")
        print(f"    WRTCMP_lat_lwr_left: {round(lwr_lat_b, 2)}")
        print(f"    WRTCMP_dx: {res}")
        print(f"    WRTCMP_dy: {res}")
    return

def set_write_component(index):
    global g_cen_lon, g_cen_lat

    if index != "LambertConformal":
        print("Write component must be set from the LambertConformal grid")
        return
    g_view[index] = "regional"

    g_cen_lon, g_cen_lat, lwr_lon, lwr_lat, upr_lon, upr_lat, xspan, yspan, extent, _, _ = get_dims(index, 'blue')
    g_cen_lat = cen_lat_adjust(g_cen_lat)
    plots_draw("setting write grid") #3
    p_src = index

    return

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
        print(f"resetting to default settings...")
        plots_draw("init") #4
    elif event.key == 'y':
        index = "LambertConformal"
        output_yaml("LambertConformal")
        x1, x2, y1, y2 = g_axis[index].get_extent()
        left_lon, left_lat = ccrs.PlateCarree().transform_point(x1, (y1+y2)/2, g_axis[index].projection)
        right_lon, right_lat = ccrs.PlateCarree().transform_point(x2, (y1+y2)/2, g_axis[index].projection)
        top_lon, top_lat = ccrs.PlateCarree().transform_point((x1+x2)/2, y2, g_axis[index].projection)
        bottom_lon, bottom_lat = ccrs.PlateCarree().transform_point((x1+x2)/2, y1, g_axis[index].projection)
        lon_span = abs(right_lon - left_lon)
        lat_span = abs(top_lat - bottom_lat)
        print(f"    UFS_domain_select.py args: "  \
              f"--cen_lon {int(g_cen_lon)} " \
              f"--cen_lat {int(g_cen_lat)} " \
              f"--lon_span {int(lon_span)} " \
              f"--lat_span {int(lat_span)}")
    elif event.key == 'R':
        if g_res == 3000:
            g_res = 13000
        elif g_res == 13000:
            g_res = 25000
        else:
            g_res = 3000
        print(f"resolution set to {int(g_res/1000)}km")
    elif event.key == '+':
        g_compute_grid += .01
        print(f"set compute grid component to {g_compute_grid}")
        set_write_component("LambertConformal")
    elif event.key == '-':
        g_compute_grid -= .01
        print(f"set compute grid component to {g_compute_grid}")
        set_write_component("LambertConformal")
    elif event.key == ' ':
        print(f"setting write component...")
        set_write_component("LambertConformal")
        print(f"write component set ok. 'y' to output yaml")
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

    #plt.show() # end on_key_press


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
    global g_lambert_xspan, g_lambert_yspan

    #print(f"projs_create: mode {mode} g_cen_lon {g_cen_lon} g_cen_lat {g_cen_lat}")

    if mode == "init":
        g_cen_lon = g_cen_lon_default
        g_cen_lat = g_cen_lat_default
        g_compute_grid = g_compute_grid_default

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

    if False and args.file:
        g_projs = [ "LambertConformal" ]
        g_dim_x = 1; g_dim_y = 1
    else:
        g_projs = [ "LambertConformal", "Gnomonic", "Mercator", "Orthographic" ]
        g_dim_x = 2; g_dim_y = 2

    if mode == "init":
        for p in g_projs:

            if p == "LambertConformal":
                lon_span = g_lon_span_default
                lat_span = g_lat_span_default
                nxc, nyc = g_proj[p].transform_point(g_cen_lon, g_cen_lat, ccrs.Geodetic())
                nxl, nyl = g_proj[p].transform_point(g_cen_lon-lon_span/2, g_cen_lat, ccrs.Geodetic())
                nxr, nyr = g_proj[p].transform_point(g_cen_lon+lon_span/2, g_cen_lat, ccrs.Geodetic())
                nxt, nyt = g_proj[p].transform_point(g_cen_lon, g_cen_lat+lat_span/2, ccrs.Geodetic())
                nxb, nyb = g_proj[p].transform_point(g_cen_lon, g_cen_lat-lat_span/2, ccrs.Geodetic())
                g_lambert_xspan = abs(nxr-nxl)
                g_lambert_yspan = abs(nyt-nyb)

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
    global g_lambert_xspan, g_lambert_yspan

    # Save mode
    if not mode == "init":
        for p in g_projs:
            g_mode = g_axis[p].get_navigate_mode()

    # Remove any old plots
    if g_projs:
        x1, x2, y1, y2 = g_axis["LambertConformal"].get_extent()
        g_lambert_xspan = abs(x2-x1)
        g_lambert_yspan = abs(y2-y1)
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

    if not args.file:

        # Set extent of Lambert grid and draw write components in all
        # projections. Draw compute grid in just the Gnomonic.
        g_axis["LambertConformal"].set_extent((-g_lambert_xspan/2, g_lambert_xspan/2, 
                                               -g_lambert_yspan/2, g_lambert_yspan/2), 
                                               crs=g_proj["LambertConformal"])

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
        for p_tgt in g_projs:
            draw_box_xy_data(p_tgt, "Gnomonic", 'red') # compute component

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
    cen_lon, cen_lat, lwr_lon, lwr_lat, upr_lon, upr_lat, xspan, yspan, extent, _, _ = \
                get_dims(src_index, color)
    x1, x2, y1, y2 = extent 
    # Star inside of circle using geodetic coords (just a test)
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

    # Filter out by lat/lon
    if args.cen_lon and args.cen_lat and args.lon_span and args.lat_span: 
        lat_min, lat_max = round(g_cen_lat-g_lat_span_default/2), \
                           round(g_cen_lat+g_lat_span_default/2)
        dim0 = (lat_max-lat_min)*4+1
        lon_min, lon_max = round(360+g_cen_lon-g_lon_span_default/2), \
                           round(360+g_cen_lon+g_lon_span_default/2)
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
plots_draw("init")
