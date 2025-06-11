#!/usr/bin/env python3

#
# TODO
#
# Bugs/tweaks/cleanup
# - FIXME: clean up plot+save.  Why is contour under plot and colormesh under save? 
# ...
# - persist projection selection when we change domain
# - make checked Orthographic global, by default
# - re-organize script comments, arguments and such
# - move plots_draw and plots_remove into Class definitions??
# - fix figure names (like "regional non-index")
# - add menu for predefined UFS regions (to be parsed from the file ;-)
# - things get sluggish when we remove LambertConformal?
#
# Validation
# - recreate various predefined grids (should be identical!)
#   whatever we do to recreate needs to be automated by reading
#   in the predefined file!
#
#                             DFLT  REPR
#   L RRFS_CONUS_25km           OK
#   L RRFS_CONUScompact_25km    OK
#   L RRFS_CONUS_13km           OK
#   L RRFS_CONUScompact_13km    OK
#   L RRFS_CONUS_3km            X - OOM
#   L RRFS_CONUScompact_3km     X - OOM
#   L SUBCONUS_Ind_3km          OK
#   L RRFS_AK_13km              OK - but contours wrap
#   L RRFS_AK_3km               X - OOM 
#   L WoFS_3km                  OK
#   R CONUS_25km_GFDLgrid       OK
#   R CONUS_3km_GFDLgrid        X - OOM 
#   L GSD_HRRR_25km             OK
#   R AQM_NA_13km               X - OOM 
#   R RRFS_NA_13km              X - OOM 
#   R RRFS_NA_3km
#   R RRFS_NA_25km              OK - but odd projection
#
# - make sure uds.compute_grid (gnomonic) is correct and add option to show it
# - test w/ GFS files (currently filters based on lon/lat spans)
#
# Features
# - make it so that all projections go to Orthographic
# - filter big GFS by region
# - update Makefile to download latest cycle and integrate into GUI
# - write status to screen (e.g., YAML written, res selected, ...)?
# - understand and integrate regional_latlon
#

import argparse
import os
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
import yaml

###########
# Globals #
###########

g_debug = False

HOME = f"{os.environ['HOME']}"
UFS_DOMAIN_SELECT_HOME=os.path.dirname(os.path.abspath(__file__))

g_res_dflt=-1                                             # 3000, 13000, 25000, or -1 (auto)
g_yaml_file = f"{UFS_DOMAIN_SELECT_HOME}/build/ufs-srweather-app-v2.2.0/ush/config.yaml" # Location of YAML output
g_compute_grid_dflt = 0.1                                 # 5% larger than write component grid


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

#####################
# Command-line args #
#####################

parser = argparse.ArgumentParser()
parser.add_argument("--cen_lon", help="center longitude", required=False)
parser.add_argument("--cen_lat", help="center latitude", required=False)
parser.add_argument("--crn_lon", help="lower left corner longitude", required=False)
parser.add_argument("--crn_lat", help="lower left corner latitude", required=False)
parser.add_argument("--file", "-f", help="grib file to plot", required=False)
parser.add_argument("--proj", "-p", help="LambertConformal or RotatedPole (required with -f)", required=False)
parser.add_argument("--close", "-x", help="close after saving plot of grib file", required=False, action="store_true")
g_args = parser.parse_args()

#####################
# Class definitions #
#####################

class grib():
    global g_index_dflt, g_cen_lon_dflt, g_cen_lat_dflt

    def init(self, file):

        debug(f"grib: loading file {file}")

        # Open file and get forecast time
        self.data = pygrib.open(file)

        msg = self.data[1]
        print(f"grib: params: {msg.projparams}")
        match msg.projparams['proj']:
            case "lcc":
                g_index_dflt = 'LambertConformal'
                g_cen_lon_dflt = msg.projparams['lon_0']
                g_cen_lat_dflt = msg.projparams['lat_0']
                print(f"grib: lon_0: {msg.projparams['lon_0']}")
                print(f"grib: lat_0: {msg.projparams['lat_0']}")
            case "ob_tran":
                g_index_dflt = 'RotatedPole'
                g_cen_lon_dflt = msg.projparams['lon_0']
                g_cen_lat_dflt = 90 - msg.projparams['o_lat_p']
            case _:
                print(f"grib: unknown projection {msg.projparams['proj']} -- using default")
        print(f"grib: proj: {msg.projparams['proj']} --> {g_index_dflt}")

        zmsg = self.data.select(name="Geopotential height", level=500)[0]
        zdata = zmsg.values
        lats, lons = zmsg.latlons()
        vdata = self.data.select(name="Absolute vorticity", level=500)[0].values
        sio = io.StringIO()
        print(zmsg, file=sio)
        fcst = sio.getvalue().split(':')[6]
        time = sio.getvalue().split(':')[7]
        self.title = fcst + " " + time

        # Get geopotential height
        z500 = zdata * 0.1
        z500 = scipy.ndimage.gaussian_filter(z500, 6.89)

        if True:
            for p in self.uds.projs:
             if self.uds.plotted[p]:
                contours = self.uds.axis[p].contour(lons, lats, z500,
                                                    np.arange(0, 900, 6),
                                                    colors="blue", linewidths=1,
                                                    transform=cartopy.crs.PlateCarree(), alpha=0.25)
            plt.clabel(contours, np.arange(0, 900, 6), inline_spacing=1, fmt="%d", fontsize=8)

        # Get vorticity
        self.vort500 = vdata * 100000
        self.vort500 = scipy.ndimage.gaussian_filter(self.vort500, 1.7225)
        if False:
            self.vort500[self.vort500 > 1000] = 0  # Mask out undefined values on domain edge
        vortlevs = [16, 20, 24, 28, 32, 36, 40]
        colorlist = ["yellow", "gold", "goldenrod", "orange", "orangered", "red" ]
        self.cm = matplotlib.colors.ListedColormap(colorlist)
        self.norm = matplotlib.colors.BoundaryNorm(vortlevs, self.cm.N)

        self.lons = lons
        self.lats = lats
        self.initialized = True

        debug(f"grib: done loading file {file}")

    def plot(self, file, parent):
        if not self.initialized:
            self.init(file)
        debug("grib: plotting")
        for p in self.uds.projs:
            if self.uds.plotted[p]:
                cs1_a = self.uds.axis[p].pcolormesh(self.lons, self.lats,
                                                    self.vort500, transform=cartopy.crs.PlateCarree(),
                                                    cmap=self.cm, norm=self.norm)
        cs1_a.cmap.set_under("none")   # under 16
        cs1_a.cmap.set_over("darkred") # above 40
        debug("grib: done plotting")

    def __init__(self, uds):
        self.uds = uds
        self.initialized = False

class ufs_domain_select:

    def set_dflts(self):
        s = self
        if g_args.cen_lon:
            s.cen_lon = float(g_args.cen_lon)
        else:
            s.cen_lon = g_cen_lon_dflt
        if g_args.cen_lat:
            s.cen_lat = float(g_args.cen_lat)
        else:
            s.cen_lat = g_cen_lat_dflt
        if g_args.crn_lon:
            s.crn_lon = float(g_args.crn_lon)
        else:
            s.crn_lon = g_crn_lon_dflt
        if g_args.crn_lat:
            s.crn_lat = float(g_args.crn_lat)
        else:
            s.crn_lat = g_crn_lat_dflt
        if g_args.proj:
            s.index = g_args.proj
        else:
            s.index = g_index_dflt
        s.compute_grid = g_compute_grid_dflt
        s.res = g_res_dflt

    def projs_create(self, mode):
        s = self
        if mode == "init":
            s.enabled = {}
            s.view = {}
            s.extent = {}
        s.globe = {}
        s.proj = {}
        s.plotted = {}
        s.color = {}
        s.projs = []
        s.axis = {}
        s.loc = {}
        s.cen_lat = cen_lat_adjust(s.cen_lat)

        # Cylindrical
        s.proj["PlateCarree"] = proj_create(s, "PlateCarree")
        s.proj["Mercator"] = proj_create(s, "Mercator")
        s.proj["Miller"] = proj_create(s, "Miller")
    
        # Conic
        s.proj["EquidistantConic"] = proj_create(s, "EquidistantConic") 
        s.proj["AlbersEqualArea"] = proj_create(s, "AlbersEqualArea")
        s.proj["LambertConformal"] = proj_create(s, "LambertConformal")
    
        # Planar/Azimuthal
        s.proj["Stereographic"] = proj_create(s, "Stereographic")
        s.proj["Orthographic"] = proj_create(s, "Orthographic")
        s.proj["Gnomonic"] = proj_create(s, "Gnomonic")
    
        # Pseudo-cylindrical
        s.proj["Robinson"] = proj_create(s, "Robinson")
        s.proj["Mollweide"] = proj_create(s, "Mollweide")
        s.proj["Sinusoidal"] = proj_create(s, "Sinusoidal")
        s.proj["InterruptedGoodeHomolosine"] = proj_create(s, "InterruptedGoodeHomolosine")
    
        # Other
        s.proj["RotatedPole"] = proj_create(s, "RotatedPole")

        if mode == "init":
            for p in s.proj:
                s.enabled[p] = False
            s.enabled[g_index_dflt] = True
            if not g_args.file:
                #s.enabled['PlateCarree'] = True
                s.enabled['Mercator'] = True
                #s.enabled['Miller'] = True
                #s.enabled['EquidistantConic'] = True
                #s.enabled['AlbersEqualArea'] = True
                s.enabled['LambertConformal'] = True
                #s.enabled['Stereographic'] = True
                s.enabled['Orthographic'] = True
                s.enabled['Gnomonic'] = True
                #s.enabled['Robinson'] = True
                #s.enabled['Mollweide'] = True
                #s.enabled['Sinusoidal'] = True
                #s.enabled['InterruptedGoodeHomolosine'] = True
                s.enabled['RotatedPole'] = True

        count = 0
        for p in s.proj:
            if s.enabled[p]:
                s.projs.append(p)
                count += 1

        match count:
            case 1:
                s.dim_x = 1; s.dim_y = 1
            case 2:
                s.dim_x = 1; s.dim_y = 2
            case 3 | 4:
                s.dim_x = 2; s.dim_y = 2
            case 5 | 6:
                s.dim_x = 2; s.dim_y = 3
            case 7 | 8:
                s.dim_x = 3; s.dim_y = 3
            case 9:
                s.dim_x = 3; s.dim_y = 3
            case 10 | 11 | 12:
                s.dim_x = 4; s.dim_y = 3
            case 13 | 14:
                s.dim_x = 4; s.dim_y = 4
    
        # Give each projection a color (brown is the default).
        for p in s.projs:
            s.color[p] = "brown"
            s.plotted[p] = False
        s.color["LambertConformal"] = "blue"
        s.color["Gnomonic"] = "red"
        s.color["Orthographic"] = "black"
        s.color["RotatedPole"] = "purple"
        s.color["Mercator"] = "yellow"
    
        if mode == "init":
            for p in s.proj:
                s.view[p] = "regional"

    def on_button_press(self, event):
        if not (event.inaxes and event.button is MouseButton.MIDDLE):
            return
    
        index = get_index(self, event.inaxes)
    
        print(f"on_button_press: centering projection ({index})")
    
        if self.view[index] == "regional" and self.globe[index] == self.axis[index].get_extent():
            print(f"EXTENT IS REGIONAL BUT GLOBAL: {self.axis[index].get_extent()}")
            #self.axis[index].set_global()
            force_global = True
        else:
            force_global = False
    
        if self.view[index] == "global":
            try:
                self.axis[index].set_extent(self.extent[index], crs=self.proj[index])
                debug(f"on_button_press: A: set extent for {index}: {self.extent[index]}")
            except:
                debug(f"on_button_press: A: failed to set extent for {index}: {self.extent[index]} -- setting global")
                self.axis[index].set_global()
            restore = True
        else:
            restore = False
    
        self.cen_lon, self.cen_lat = ccrs.PlateCarree().transform_point(event.xdata, event.ydata,
                                                                  event.inaxes.projection)
        x1, x2, y1, y2 = self.axis[index].get_extent()
        if no_cen_lat(index):
            yc = event.ydata
            _, _, lower, upper = self.globe[index]
            if yc-(y2-y1)/2>=lower and yc+(y2-y1)/2<=upper:
                self.extent[index] = -(x2-x1)/2, +(x2-x1)/2, yc-(y2-y1)/2, yc+(y2-y1)/2
            else:
                if yc-(y2-y1)/2<lower:
                    print(f"HIT BOTTOM EDGE: {index}")
                if yc+(y2-y1)/2>upper:
                    print(f"HIT TOP EDGE: {index}")
                print(f"self.globe{index} = {self.globe[index]}")
                self.extent[index] = -(x2-x1)/2, +(x2-x1)/2, y1, y2
        else:
            self.extent[index] = -(x2-x1)/2, +(x2-x1)/2, -(y2-y1)/2, +(y2-y1)/2
    
        self.axis[index].remove()
        self.proj[index] = proj_create(self, index)
        self.axis[index] = self.fig.add_subplot(self.dim_x, self.dim_y, self.loc[index], 
                                          projection=self.proj[index])
        self.axis[index].margins(x=0.0, y=0.0)
        self.axis[index].add_feature(cfeature.COASTLINE)
        self.axis[index].add_feature(cfeature.BORDERS)
        self.axis[index].add_feature(cfeature.STATES)
        self.axis[index].gridlines()
        self.axis[index].set_title(index + " (centered)")
        self.axis[index].plot(*(self.cen_lon, self.cen_lat), transform=ccrs.Geodetic(), 
                             marker='*', ms=20, color='green')
        if force_global:
            self.axis[index].set_global()
            self.axis[index].set_title(index + " (centered forced global)")
        else:
            try:
                self.axis[index].set_extent(self.extent[index], crs=self.proj[index])
                debug(f"on_button_press: B: set extent for {index}: {self.extent[index]}")
            except:
                debug(f"on_button_press: B: failed to set extent for {index}: {self.extent[index]} -- setting global")
                self.axis[index].set_global()
    
        #self.axis[index].callbacks.connect('xlim_changed', on_xlim_changed)
        #self.axis[index].callbacks.connect('ylim_changed', on_ylim_changed)
    
        if restore:
            self.axis[index].set_global()
            self.axis[index].set_title(index + " (global center restored)")
            debug(f"restore: global extent is {self.axis[index].get_extent()}")
    
        # Now overlay the actual center
        x1, x2, y1, y2 = self.axis[index].get_extent()
        cen_lon, cen_lat = ccrs.PlateCarree().transform_point((x1+x2)/2, (y1+y2)/2,
                                                               self.axis[index].projection)
        self.axis[index].plot(*(cen_lon, cen_lat), transform=ccrs.Geodetic(), 
                             marker='*', ms=10, color='red')
        plt.draw()

    def on_key_press(self, event):
        self.index = g_index_dflt
        if event.inaxes:
            self.index = get_index(self, event.inaxes)
    
        if event.key == 'h':
            show_help()
        elif event.key == '0':
            print(f"on_key_press: restoring default settings")
            plots_draw(self, "init")
        elif event.key == 'y':
            if (self.index == "LambertConformal" or self.index == "RotatedPole"):
                output_config(self, self.index)
            else:
                print("choose from LambertConformal or RotatedPole")
                return
        elif event.key == 'R':
            if self.res == 3000:
                self.res = 13000
            elif self.res == 13000:
                self.res = 25000
            elif self.res == 25000:
                self.res = -1 
            elif self.res == -1:
                self.res = 3000
            if (self.res == -1):
                print(f"on_key_press: resolution set to AUTO")
            else:
                print(f"on_key_press: resolution set to {int(self.res/1000)}km")
        elif event.key == '+':
            self.compute_grid += .01
            percent = round(100*(1+self.compute_grid))
            print(f"on_key_press: set compute grid to {percent}% of write component")
            plots_draw(self, "set")
        elif event.key == '-':
            self.compute_grid -= .01
            percent = round(100*(1+self.compute_grid))
            print(f"on_key_press: set compute grid to {percent}% of write component")
            plots_draw(self, "set")
        elif event.key == ' ':
            if not event.inaxes:
                print(f"on_key_press: hover over an axis to select source projection)")
                return
            print(f"on_key_press: setting source projection ({self.index})")
            if (self.axis[self.index].get_extent() == self.globe[self.index]):
                print(f"on_key_press: failed attempt to set with global extent: {self.globe[self.index]}")
                return
            plots_draw(self, "set")
        elif event.key == 'g':
            print(f"toggling global view (source {self.index})")
            if self.view[self.index] == "regional":
                self.extent[self.index] = self.axis[self.index].get_extent()
                self.view[self.index] = "global"
                self.axis[self.index].set_global()
                self.axis[self.index].set_title(self.index + " (global toggle)")
            elif self.view[self.index] == "global": 
                self.view[self.index] = "regional"
                self.axis[self.index].set_extent(self.extent[self.index], crs=self.proj[self.index])
                self.axis[self.index].set_title(self.index + " (regional toggle)")
        elif event.key == 'q':
            exit(0)

    def checkfunc(self, label):

        # get index (i) of projection in self.proj
        i = 0
        for p in self.proj:
            if p == label:
                break
            i = i + 1

        if len(self.projs) == 1:
            true_count = sum(self.check1.get_status())
            if not self.check1.get_status()[i] and true_count == 0:
                print(f"CANNOT DISABLE LAST PROJECTION ({label})")
                self.check1.set_active(i)
                return
            elif self.check1.get_status()[i] and true_count == 1:
                print(f"RE-ENABLING LAST PROJECTION ({label})")
                # the above if statement will call back to this point
                # to re-enabled the button, so just return
                return

        if not self.enabled[label]:
            self.enabled[label] =  True
        else:
            if label == self.index:
                self.enabled[label] =  False
                self.index = self.projs[0]
                debug(f"checkfunc: changed index to {self.index}")
                if label == self.index:
                    self.index = self.projs[1]
                    debug(f"checkfunc: changed index again to {self.index}")
                self.view[self.index] = "regional"
                self.axis[self.index].set_extent(self.extent[self.index], crs=self.proj[self.index])
                debug(f"checkfunc: forcing {self.index} to regional")
            else:
                debug(f"checkfunc: forcing {self.index} to regional")
                self.view[self.index] = "regional"
                self.axis[self.index].set_extent(self.extent[self.index], crs=self.proj[self.index])
            self.enabled[label] = False

        plots_draw(self, "set")

    def __init__(self):
        plt.rcParams["figure.raise_window"] = False
        self.fig = plt.figure(figsize=(8, 5))
        self.fig_control = plt.figure(figsize=(5.5, 5))
        self.fig.canvas.mpl_connect('button_press_event', self.on_button_press)
        self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)
        self.grib = grib(self)

########################
# Function definitions #
########################

def debug(format):
    if g_debug:
        print(format)

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

def fmt_tuple(tuple_in):
    return tuple([str(round(x,2)) if isinstance(x, float) else x for x in tuple_in])

def get_dims(uds, index, color):
    extent = uds.axis[index].get_extent()
    x1, x2, y1, y2 = extent

    if color == "red":
        xspan = abs(x2-x1)
        x1 -= xspan*uds.compute_grid/2
        x2 += xspan*uds.compute_grid/2
        yspan = abs(y2-y1)
        y1 -= yspan*uds.compute_grid/2
        y2 += yspan*uds.compute_grid/2
        extent = x1, x2, y1, y2

    xc = (x1+x2)/2
    yc = (y1+y2)/2
    xspan = abs(x2-x1)
    yspan = abs(y2-y1)
    cen_lon, cen_lat = ccrs.PlateCarree().transform_point(xc, yc, uds.axis[index].projection)
    lwr_lon, lwr_lat = ccrs.PlateCarree().transform_point(x1, y1, uds.axis[index].projection)
    upr_lon, upr_lat = ccrs.PlateCarree().transform_point(x2, y2, uds.axis[index].projection)

    return cen_lon, cen_lat, lwr_lon, lwr_lat, extent

def pick_max_delta(xspan, yspan):
    if (xspan/25000 < 50) or (yspan/25000<50):
        if (xspan/13000 < 50) or (yspan/13000<50):
            return 3000
        else:
            return 13000
    else:
        return 25000

def output_config(uds, index):
    res = uds.res

    # computational grid
    cen_lon_r, cen_lat_r, lwr_lon_r, lwr_lat_r, extent_r = get_dims(uds, index, "red")
    x1, x2, y1, y2 = extent_r
    xspan_r = abs(x2-x1)
    yspan_r = abs(y2-y1)
    if uds.index == "RotatedPole":
        max_delta = pick_max_delta(xspan_r*1852*60, yspan_r*1852*60)
    else:
        max_delta = pick_max_delta(xspan_r, yspan_r)
    if uds.res == -1:
        res = max_delta
        print(f"             resolution: automatically set to {res} meters")
    elif (uds.res>max_delta):
        res = max_delta
        print(f"             resolution: forced to max ({res} meters)")
    percent = round(100*(1+uds.compute_grid))
    print(f"      compute grid size: {percent}% of write component")
    if uds.index == "RotatedPole":
        nx_r = int(xspan_r*1852*60/res)
        ny_r = int(yspan_r*1852*60/res)
    else:
        nx_r = int(xspan_r/res)
        ny_r = int(yspan_r/res)

    print(f"  gnomonic compute grid: nx {nx_r} ny {ny_r}")

    # write component grid
    cen_lon_b, cen_lat_b, lwr_lon_b, lwr_lat_b, extent_b = get_dims(uds, index, "blue")
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
  ESGgrid_LON_CTR: {uds.cen_lon}
  ESGgrid_LAT_CTR: {uds.cen_lat}
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

    if uds.index == "LambertConformal":
        config_text_wrtcmp = f"""
  WRTCMP_output_grid: "lambert_conformal"
  WRTCMP_cen_lon: {uds.cen_lon}
  WRTCMP_cen_lat: {uds.cen_lat} 
  WRTCMP_stdlat1: {uds.cen_lat} 
  WRTCMP_stdlat2: {uds.cen_lat} 
  WRTCMP_nx: {nx_b}
  WRTCMP_ny: {ny_b}
  WRTCMP_lon_lwr_left: {round(lwr_lon_b, 2)} 
  WRTCMP_lat_lwr_left: {round(lwr_lat_b, 2)}
  WRTCMP_dx: {res}
  WRTCMP_dy: {res}
"""

    elif uds.index == "RotatedPole":
        config_text_wrtcmp = f"""
  WRTCMP_output_grid: "rotated_latlon"
  WRTCMP_cen_lon: {uds.cen_lon}
  WRTCMP_cen_lat: {uds.cen_lat} 
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

def get_index(uds, ax):
    for p in uds.projs:
        if uds.axis[p] == ax:
            return p

def proj_create(self, index):
    match index:

        # Cylindrical
        case "PlateCarree":
            return ccrs.PlateCarree(central_longitude=self.cen_lon)
        case "Mercator":
            return ccrs.Mercator(central_longitude=self.cen_lon)
        case "Miller":
            return ccrs.Miller(central_longitude=self.cen_lon)

        # Pseudo Cylindrical
        case "Sinusoidal":
            return ccrs.Sinusoidal(central_longitude=self.cen_lon)
        case "Mollweide":
            return ccrs.Mollweide(central_longitude=self.cen_lon)
        case "Robinson":
            return ccrs.Robinson(central_longitude=self.cen_lon)
        case "InterruptedGoodeHomolosine":
            return ccrs.InterruptedGoodeHomolosine(central_longitude=self.cen_lon)

        case "EquidistantConic":
            return ccrs.EquidistantConic(central_longitude=self.cen_lon, central_latitude=self.cen_lat)
        case "AlbersEqualArea":
            return ccrs.AlbersEqualArea(central_longitude=self.cen_lon, central_latitude=self.cen_lat)
        case "LambertConformal":
            if self.cen_lat>0:
                lat_cutoff = -30
            else:
                lat_cutoff = 30
            return ccrs.LambertConformal(central_longitude=self.cen_lon, 
                                         central_latitude=self.cen_lat,
                                         standard_parallels=(self.cen_lat, self.cen_lat), 
                                         cutoff=lat_cutoff)
        case "Stereographic":
            return ccrs.Stereographic(central_longitude=self.cen_lon, central_latitude=self.cen_lat)
        case "Orthographic":
            return ccrs.Orthographic(central_longitude=self.cen_lon, central_latitude=self.cen_lat)
        case "Gnomonic":
            return ccrs.Gnomonic(central_longitude=self.cen_lon, central_latitude=self.cen_lat)
        case "RotatedPole":
            return ccrs.RotatedPole(pole_latitude=90-self.cen_lat, pole_longitude=self.cen_lon-180)

def no_cen_lat(index):
    if index == "Mercator" or index == "Miller" or index == "PlateCarree" \
            or index == "Robinson" or index == "Sinusoidal" \
            or index == "Mollweide" or index == "InterruptedGoodeHomolosine":
            return True
    else:
            return False

def plots_remove(uds):
    if uds.projs:
        for p in uds.projs:
            if uds.plotted[p]:
                print(f"plots_remove: removing {p} axis")
                uds.axis[p].remove()
                uds.plotted[p] = False

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

def plots_draw(uds, mode):

    #
    # Notes on variables:
    #
    # uds.index is the index of the source axis that we're trransforming from.
    # uds.projs is the list of projections, only one of which is uds.index.
    #

    # Change all target axes to regional view. Any target
    # axis in global view will be returned to global view
    # at the end of this procedure.
    if not mode == "init":
        restore_global = {}
        for p in uds.projs:
            if p == uds.index:
                uds.extent[p] = uds.axis[p].get_extent()
            elif (not p == uds.index) and uds.plotted[p] and uds.enabled[p]:
                try:
                    uds.axis[p].set_extent(uds.extent[p], crs=uds.proj[p])
                except:
                    print(f"*** WARN: FAILED TO SET EXTENT ({p}) ***")
                if uds.view[p] == "global":
                    restore_global[p] = True
                uds.view[p] = "regional"

    if mode == "init":
        uds.set_dflts()
    else:
        if mode == "set":
            uds.cen_lon, uds.cen_lat, uds.crn_lon, uds.crn_lat, (x1, x2, y1, y2) = get_dims(uds, uds.index, 'blue')
        elif mode == "center":
            _, _, _, _, (x1, x2, y1, y2) = get_dims(uds, uds.index, 'blue')

        if no_cen_lat(uds.index):
            xc, yc = uds.proj[uds.index].transform_point(uds.cen_lon, uds.cen_lat, ccrs.Geodetic())
            # This calculation doesn't work for InterruptedGoodeHomolosine
            if uds.index == "InterruptedGoodeHomolosine":
                print(f"plots_draw: special case for interruptedGoodeHomolosine (set)")
                new_extent = (xc-abs(x2-x1)/2, xc+abs(x2-x1)/2, yc-abs(y2-y1)/2, yc+abs(y2-y1)/2)
            else:
                new_extent = (-abs(x2-x1)/2, abs(x2-x1)/2, yc-abs(y2-y1)/2, yc+abs(y2-y1)/2)
        else:
            new_extent = (-abs(x2-x1)/2, abs(x2-x1)/2,   -abs(y2-y1)/2,    abs(y2-y1)/2)

    # Remove old plots (the try will fail on the first "init", as uds.projs hasn't
    # yet been defined via projs_create()
    try:
        plots_remove(uds)
    except:
        debug("plots_draw: no plots to remove")

    # Create projections
    uds.projs_create(mode)

    # Create plots
    j = 1

    debug(f"plots_draw: creating plots W/ index {uds.index}")

    for p in uds.projs:

        uds.axis[p] = uds.fig.add_subplot(uds.dim_x, uds.dim_y, j, projection=uds.proj[p])
        uds.axis[p].margins(x=0.0, y=0.0)
        uds.axis[p].add_feature(cfeature.COASTLINE)
        uds.axis[p].add_feature(cfeature.BORDERS)
        uds.axis[p].add_feature(cfeature.STATES)
        uds.axis[p].gridlines()
        uds.loc[p] = j

        # Set axis view (global or regional/zoomed)
        debug(f"plots_draw: p is {p}, mode is {mode}, uds.view[{p}] is {uds.view[p]}")

        if ((p != uds.index) and uds.view[p] == "global"):
            uds.axis[p].set_title(f"{p} ({mode} global)")
            uds.axis[p].set_global()

        else:
            uds.view[p] = "regional"
            if p == uds.index:
                uds.axis[p].set_title(f"{p} (regional index)")
            else:
                uds.axis[p].set_title(f"{p} (regional non-index)")

        uds.plotted[p] = True

        j = j + 1

    if mode == "init":
        if uds.index == "RotatedPole":
            # Transform the RotatedPole corner to Geodetic
            print("plots_draw: init for RotatedPole")
            uds.crn_lon, uds.crn_lat = ccrs.PlateCarree().transform_point(g_crn_lon_dflt, g_crn_lat_dflt, uds.proj[uds.index])
            print(f"plots_draw: rotated corner: lon {g_crn_lon_dflt} lat {g_crn_lat_dflt}")
            print(f"plots_draw: geodetic corner: lon {uds.crn_lon} lat {uds.crn_lat}")

        xc, yc = uds.proj[uds.index].transform_point(uds.cen_lon, uds.cen_lat, ccrs.Geodetic())
        xll, yll = uds.proj[uds.index].transform_point(uds.crn_lon, uds.crn_lat, ccrs.Geodetic())
        xlr = xc+(xc-xll); ylr = yc-(yc-yll)
        xul = xll; yul = yc+(yc-yll)
        xur = xlr; yur = yul
        uds.extent[uds.index] = (xll, xlr, yll, yul)
    else:
        uds.extent[uds.index] = new_extent
        debug("plots_draw: CASE 2")

    # FIXME: plot and save (before)
    if g_args.file:
        uds.grib.plot(g_args.file, uds)
        for p in uds.projs:
            if uds.plotted[p]:
                uds.axis[p].set_title(p)
        uds.fig.suptitle(uds.grib.title)
        print(f"plots_draw: BEFORE: saving forecast to {g_args.file + '.png'} ***")
        uds.fig.savefig(g_args.file + ".A.png")
        if False and g_args.close:
            exit(0)

    if True:

        assert(uds.view[uds.index] == "regional")
        assert(mode=="set" or mode=="init")

        try:
            if not g_args.file:
                print(f"plots_draw: grib file not specified: setting default extent for {uds.index}")
                uds.axis[uds.index].set_extent(uds.extent[uds.index], crs=uds.proj[uds.index])
                debug(f"plots_draw: A: set extent to {fmt_tuple(uds.extent[uds.index])}")
            else:
                debug(f"plots_draw: grib file specified: not setting extent for {uds.index}")
            uds.axis[uds.index].set_title(f"{uds.index} (STILL REGIONAL after mode {mode})")
        except:
            print(f"*** CASE 2: FAILED TO SET EXTENT ({uds.index}): {uds.extent[uds.index]} ***")

        # Outline the extent of the index's axis
        x, y = create_box_xy_data(uds, uds.index, uds.color[uds.index])
        uds.axis[uds.index].plot(x, y, color=uds.color[uds.index], linewidth=2, alpha=1.0, linestyle='dashed')
        uds.axis[uds.index].plot(*(uds.cen_lon, uds.cen_lat), transform=ccrs.Geodetic(), 
                               marker='*', ms=10, color='purple')

        targets = [t for t in uds.projs if not t == uds.index]
        for p in targets:

            style = 'solid' 

            # Draw the transformed index's extent on the target's axis
            pts = uds.proj[p].transform_points(uds.proj[uds.index], np.array(x), np.array(y))
            tx = pts[:, 0]; ty = pts[:, 1]
            uds.axis[p].plot(tx, ty, color=uds.color[uds.index], linewidth=2, alpha=1.0, linestyle=style)

            # Outline the extent of the target's axis
            min_x, max_x, min_y, max_y = find_extent(tx, ty)
            xp, yp = create_box_xy((min_x, max_x, min_y, max_y))
            uds.axis[p].plot(xp, yp, color=uds.color[p], linewidth=2, alpha=1.0, linestyle=style)

            # Draw the transformed targets's extent on the index's axis
            pts = uds.proj[uds.index].transform_points(uds.proj[p], np.array(xp), np.array(yp))
            tx = pts[:, 0]; ty = pts[:, 1]
            uds.axis[uds.index].plot(tx, ty, color=uds.color[p], linewidth=2, alpha=1.0, linestyle=style)

    if mode == "center":
        x1, x2, y1, y2 = uds.axis[uds.index].get_extent()
        uds.crn_lon, uds.crn_lat = ccrs.PlateCarree().transform_point(x1, y1, uds.proj[uds.index])

    # Set/restore view
    if mode == "init":
        for p in uds.projs:
            if uds.plotted[p]:
                if p == "Orthographic":
                    print(f"plots_draw: special case: initializing global view for Orthographic")
                    uds.extent[p] = uds.axis[p].get_extent()
                    uds.view[p] = "global"
                    uds.axis[p].set_global()
                    uds.axis[p].set_title(p + " (global init)")
    else:
        for p in restore_global:
            if uds.enabled[p]:
                debug(f"plots_draw: restoring global view for {p} ***")
                uds.extent[p] =  uds.axis[p].get_extent()
                uds.axis[p].set_global()
                uds.view[p] = "global"
                uds.axis[p].set_title(p + " (global restore)")

    # Keep track of what the global extent is
    for p in uds.projs:
        debug(f"plots_draw: recording global extent for {p} (current view {uds.view[p]})")
        if uds.view[p] == "global":
            uds.globe[p] = uds.axis[p].get_extent()
        else:
            uds.extent[p] = uds.axis[p].get_extent()
            uds.axis[p].set_global()
            uds.globe[p] = uds.axis[p].get_extent()
            try:
                uds.axis[p].set_extent(uds.extent[p], crs=uds.proj[p])
            except:
                print(f"plots_draw: *** FAILED TO SET EXTENT FOR {p} ***: {uds.extent[p]}")

    # Check buttons
    uds.axis['menu1'] = uds.fig_control.add_axes([0.10, 0.0, 0.10, 1.0], frameon=False)
    uds.axis['menu3'] = uds.fig_control.add_axes([0.50, 0.0, 0.25, 1.0], frameon=False)

    proj1 = []
    proj2 = []
    status1 = []
    count = 1
    for p in uds.proj:
        if uds.enabled[p]:
            status1.append(True)
        else:
            status1.append(False)
        proj1.append(p)
        count += 1

    uds.check1 = CheckButtons(uds.axis['menu1'], proj1, status1)
    uds.check1.on_clicked(uds.checkfunc)
    radio = RadioButtons(uds.axis['menu3'], g_radio_buttons)
    radio.on_clicked(radio_func)

    # FIXME: plot and save (after)
    if mode == "init":
        if g_args.file:
            print(f"plots_draw: AFTER: saving forecast to {g_args.file + '.png'} ***")
            uds.fig.savefig(g_args.file + ".B.png")
            if g_args.close:
                exit(0)
        plt.show()
    else:
        uds.fig_control.canvas.draw()
        uds.fig.canvas.draw()

def radio_func(region):
    global g_cen_lon_dflt
    global g_cen_lat_dflt
    global g_crn_lon_dflt
    global g_crn_lat_dflt
    global g_index_dflt
    print(f"radio_func ({region})")
    if not g_args.file:
        g_cen_lon_dflt = WRTCMP_cen_lon[region]['WRTCMP_cen_lon']
        g_cen_lat_dflt = WRTCMP_cen_lat[region]['WRTCMP_cen_lat']
        g_crn_lon_dflt = WRTCMP_lon_lwr_left[region]['WRTCMP_lon_lwr_left']
        g_crn_lat_dflt = WRTCMP_lat_lwr_left[region]['WRTCMP_lat_lwr_left']
        match WRTCMP_output_grid[region]['WRTCMP_output_grid']:
            case "lambert_conformal":
                print("SELECTING lambert_conformal")
                g_index_dflt = "LambertConformal"
            case "rotated_latlon":
                print("SELECTING rotated_latlon")
                g_index_dflt = "RotatedPole"
    plots_draw(myuds, "init")
    print("radio_func: drawing control canvas")
    myuds.fig_control.canvas.draw()
    print("radio_func: done drawing control canvas")

def create_box_xy(extent):
    x1, x2, y1, y2 = extent 
    xs = []
    ys = []
    steps = 256
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

def create_box_xy_data(uds, src_index, color):
    cen_lon, cen_lat, lwr_lon, lwr_lat, extent = get_dims(uds, src_index, color)
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
    cen_lon, cen_lat, lwr_lon, lwr_lat, extent = get_dims(uds, src_index, color)
    x1, x2, y1, y2 = extent 
    if color == "blue":
        uds.axis[tgt_index].plot(*(cen_lon, cen_lat), transform=ccrs.Geodetic(), 
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

    uds.axis[tgt_index].plot(xs, ys, transform=uds.axis[src_index].projection, 
                           color=color, linewidth=2, alpha=1.0, linestyle='solid')

    return xs, ys

def show_help():
    print(g_help)

########
# Main #
########

with open("/home/mmesnie/UFS_domain_select/build/ufs-srweather-app-v2.2.0/ush/predef_grid_params.yaml", 'r') as file:
    yaml_data = yaml.safe_load(file)
    WRTCMP_cen_lon = {}
    WRTCMP_cen_lat = {}
    WRTCMP_lon_lwr_left = {}
    WRTCMP_lat_lwr_left = {}
    WRTCMP_output_grid = {}
    g_radio_buttons = []
    for region in yaml_data:
        WRTCMP_cen_lon[region] = {}
        WRTCMP_cen_lat[region] = {}
        WRTCMP_lon_lwr_left[region] = {}
        WRTCMP_lat_lwr_left[region] = {}
        WRTCMP_output_grid[region] = {}
        print(f"region is {region} ({yaml_data[region]['QUILTING']['WRTCMP_output_grid']})")
        for key in yaml_data[region]['QUILTING']:
            WRTCMP_output_grid[region][key] = yaml_data[region]['QUILTING']['WRTCMP_output_grid']
            WRTCMP_cen_lon[region][key] = yaml_data[region]['QUILTING']['WRTCMP_cen_lon']
            WRTCMP_cen_lat[region][key] = yaml_data[region]['QUILTING']['WRTCMP_cen_lat']
            WRTCMP_lon_lwr_left[region][key] = yaml_data[region]['QUILTING']['WRTCMP_lon_lwr_left']
            WRTCMP_lat_lwr_left[region][key] = yaml_data[region]['QUILTING']['WRTCMP_lat_lwr_left']
        g_radio_buttons.append(region)

# LambertConformal
#region='RRFS_CONUS_25km'
#region='RRFS_CONUScompact_25km'
#region='RRFS_CONUS_13km'
#region='RRFS_CONUScompact_13km'
#region='RRFS_CONUS_3km'
#region='RRFS_CONUScompact_3km'
#region='SUBCONUS_Ind_3km'
#region='RRFS_AK_13km'
#region='RRFS_AK_3km'
#region='WoFS_3km'
#region='GSD_HRRR_25km'

# RotatedPole
#region='CONUS_25km_GFDLgrid'
#region='CONUS_3km_GFDLgrid'
#region='AQM_NA_13km'
#region='RRFS_NA_13km'
#region='RRFS_NA_3km'

show_help()
myuds = ufs_domain_select()

if True:
    #g_index_dflt = 'RotatedPole'
    #g_cen_lon_dflt=-59.5;   g_cen_lat_dflt=-51.7; g_crn_lon_dflt=-61.98;  g_crn_lat_dflt=-52.81 # Falkland Islands
    #g_cen_lon_dflt=-127.68; g_cen_lat_dflt=45.72; g_crn_lon_dflt=-132.86; g_crn_lat_dflt=41.77  # Oregon coast
    #g_cen_lon_dflt=-141.87; g_cen_lat_dflt=40.48; g_crn_lon_dflt=-160.29; g_crn_lat_dflt=16.64  # Eastern Pacific 
    g_index_dflt = 'LambertConformal'
    g_cen_lon_dflt=-61.13;   g_cen_lat_dflt=10.65;  g_crn_lon_dflt=-61.98; g_crn_lat_dflt=9.85  # Trinidad & Tobago
    plots_draw(myuds, "init")
else:
    #radio_func('RRFS_CONUS_25km') # This will call plots_draw()
    radio_func('SUBCONUS_Ind_3km') # This will call plots_draw()
