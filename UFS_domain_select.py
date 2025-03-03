import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.backend_bases import MouseButton
import time
import math

g_write_component_default = 0.97
g_res = 3000 # choose among 3000, 13000 or 25000

#
# This script automates the creation of the write
# component parameters. 'h' for detailed options.
#
# -The write component grid is shown in blue on all projections.
# -The write component grid will be smaller than the Lambert plot
#  extents, as per g_write_component. E.g., if g_write_component
#  is 0.98, the write component grid will be 98% of the size of the
#  Lambert plot window.
# -The compute grid will be shown in red on the Gnomonic plot. The
#  compute grid represents 100% of the Lambert plot extents and will
#  be 100/g_write_component percent larger than the write component
#  grid.
#

#
# Globals
#

# User globals
g_cen_lon_default = -97.5; g_cen_lat_default =  38.5; g_lon_span_default = 25; g_lat_span_default = 15 # CONUS
#g_cen_lon_default = -55.86; g_cen_lat_default = -63.21; g_lon_span_default = 15; g_lat_span_default = 15 # Antarctica
#g_cen_lon_default = 16.73; g_cen_lat_default = 79.15; g_lon_span_default = 15; g_lat_span_default = 5 # Svalbard 
#g_cen_lon_default = -59.5236; g_cen_lat_default = -51.7963; g_lon_span_default = 15; g_lat_span_default = 15 # FALKLAND ISLANDS
#g_cen_lon_default = 0; g_cen_lat_default =  0; g_lon_span_default = 30; g_lat_span_default = 20
#g_cen_lon_default = 0; g_cen_lat_default = 85; g_lon_span_default = 30; g_lat_span_default = 20

g_halo_width = 6
g_dt_atmos = 36
g_blocksize = 40
g_layout_x = 3
g_layout_y = 3
g_write_groups = 1
g_write_tasks_per_group = 3

g_default_view = "regional" # regional or global

# Internal globals
g_cen_lon = g_cen_lon_default
g_cen_lat = g_cen_lat_default
g_extent = {}
g_view = {}
g_axis = {}
g_proj = {}
g_projs = {}
g_mode = None
g_xdata_span = {}
g_ydata_span = {}

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
        plots_draw("centering") #2
        print("Plots centered. Spacebar to set write component.")

def get_dims(index, color):
    extent = g_axis[index].get_extent()
    x1, x2, y1, y2 = extent

    if color == "blue":
        xspan = abs(x2-x1)
        x1 += xspan*(1-g_write_component)/2
        x2 -= xspan*(1-g_write_component)/2
        yspan = abs(y2-y1)
        y1 += yspan*(1-g_write_component)/2
        y2 -= yspan*(1-g_write_component)/2
        extent = x1, x2, y1, y2

    nxc = (x1+x2)/2
    nyc = (y1+y2)/2
    xspan = abs(x2-x1)
    yspan = abs(y2-y1)
    cen_lon, cen_lat = ccrs.PlateCarree().transform_point(nxc, nyc, g_axis[index].projection)
    lwr_lon, lwr_lat = ccrs.PlateCarree().transform_point(x1, y1, g_axis[index].projection)
    upr_lon, upr_lat = ccrs.PlateCarree().transform_point(x2, y2, g_axis[index].projection)

    return cen_lon, cen_lat, lwr_lon, lwr_lat, upr_lon, upr_lat, xspan, yspan, extent

def set_views():
    for p in g_projs:
        if g_view[p] == "global":
            g_axis[p].set_title(f"{p} global")
            g_axis[p].set_global()
        else:
            g_axis[p].set_title(f"{p} regional")
            g_axis[p].set_extent(g_extent[p], crs=g_axis[p].projection)
            print(f"RESTORED g_extent[{p}]")

def output_yaml(index):
    if index != "LambertConformal":
        print("Write component must be written from LambertConformal grid")
        return
    # computational grid
    cen_lon_r, cen_lat_r, lwr_lon_r, lwr_lat_r, \
        upr_lon_r, upr_lat_r, xspan_r, yspan_r, extent_r = get_dims(index, "red")
    # write component grid
    cen_lon_b, cen_lat_b, lwr_lon_b, lwr_lat_b, \
        upr_lon_b, upr_lat_b, xspan_b, yspan_b, extent_b = get_dims(index, "blue")
    for res in [ g_res ]:
        nx_r = int(xspan_r/res)
        ny_r = int(yspan_r/res)
        print(f"\"CUSTOM_XX_{int(res/1000)}km\":")
        print(f"  GRID_GEN_METHOD: \"ESGgrid\"")
        print(f"  ESGgrid_LON_CTR: {cen_lon_r}")
        print(f"  ESGgrid_LAT_CTR: {cen_lat_r}")
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
        print(f"    WRTCMP_cen_lon: {cen_lon_b}")
        print(f"    WRTCMP_cen_lat: {cen_lat_b}")
        print(f"    WRTCMP_stdlat1: {cen_lat_b}")
        print(f"    WRTCMP_stdlat2: {cen_lat_b}")
        print(f"    WRTCMP_nx: {nx_b}")
        print(f"    WRTCMP_ny: {ny_b}")
        print(f"    WRTCMP_lon_lwr_left: {lwr_lon_b}")
        print(f"    WRTCMP_lat_lwr_left: {lwr_lat_b}")
        print(f"    WRTCMP_dx: {res}")
        print(f"    WRTCMP_dy: {res}")
    return

def set_write_component(index):
    global g_cen_lon, g_cen_lat

    if index != "LambertConformal":
        print("Write component must be set from the LambertConformal grid")
        return
    g_view[index] = "regional"
    g_cen_lon, g_cen_lat, lwr_lon, lwr_lat, upr_lon, upr_lat, xspan, yspan, extent = get_dims(index, 'blue')
    if False:
        print(f"cen_lon {g_cen_lon} cen_lat {g_cen_lat}")
        print(f"lwr_lon {lwr_lon} lwr_lat {lwr_lat}")
        print(f"upr_lon {upr_lon} upr_lat {upr_lat}")
        print(f"xspan {xspan} yspan {yspan}")
    g_cen_lat = cen_lat_adjust(g_cen_lat)
    plots_draw("setting write grid") #3
    p_src = index
    for p_tgt in g_projs:
        draw_box_xy_data(p_tgt, p_src, 'blue') # write component grid
    draw_box_xy_data("Gnomonic", "Gnomonic", 'red') # write component grid

    return

def on_key_press(event):
    global g_cen_lon, g_cen_lat
    global g_write_component
    global g_res

    ax = event.inaxes
    index = get_index(ax)

    if event.key == 'h':
        show_help()
    elif event.key == '0':
        plots_draw("init") #4
    elif event.key == 'v':
        set_views()
    elif event.key == 'y':
        output_yaml("LambertConformal")
    elif event.key == 'R':
        if g_res == 3000:
            g_res = 13000
        elif g_res == 13000:
            g_res = 25000
        else:
            g_res = 3000
        print(f"Resolution set to {int(g_res/1000)}km")
    elif event.key == '+':
        g_write_component += .01
        print(f"Set write component to {int(100*g_write_component)}% of compute grid")
        set_write_component("LambertConformal")
    elif event.key == '-':
        g_write_component -= .01
        print(f"Set write component to {int(100*g_write_component)}% of compute grid")
        set_write_component("LambertConformal")
    elif event.key == ' ':
        print(f"Write component set. \'y\' to output yaml")
        set_write_component("LambertConformal")
    elif event.key == 'g':
        print("DOING REGIONAL")
        if g_view[index] == "regional":
            _, _, _, _, _, _, _, _, g_extent[index] = get_dims(index, "red")
            print(f"SAVED g_extent[{index}]")
            g_view[index] = "global"
            set_views()
        else:
            print(f"ALREADY GLOBAL")
    elif event.key == 'r':
        print("DOING REGIONAL")
        if g_view[index] == "global":
            g_view[index] = "regional"
            set_views()
        else:
            print(f"ALREADY REGIONAL")
    elif event.key == 'q':
        exit(0)

def get_index(ax):
    for p in g_projs:
        if g_axis[p] == ax:
            return p

def projs_create(mode):
    global g_proj
    global g_projs
    global g_cen_lat
    global g_cen_lon
    global g_view, g_extent
    global g_dim_x, g_dim_y
    global g_write_component

    #print(f"projs_create: mode {mode} g_cen_lon {g_cen_lon} g_cen_lat {g_cen_lat}")

    if mode == "init":
        g_cen_lon = g_cen_lon_default
        g_cen_lat = g_cen_lat_default
        g_write_component = g_write_component_default

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
    g_projs = [ "LambertConformal", "Gnomonic", "Orthographic" ]
    g_dim_x = 3; g_dim_y = 1

    if mode == "init":
        for p in g_projs:
            if p == "Orthographic":
                g_view[p] = "global"
            else:
                g_view[p] = g_default_view
            print("WARN: lon/lat spans should use a transform")
            if p == "PlateCarree":
                g_xdata_span[p] = g_lon_span_default
                g_ydata_span[p] = g_lat_span_default
            else:
                g_xdata_span[p] = g_lon_span_default * 60 * 1852
                g_ydata_span[p] = g_lat_span_default * 60 * 1852

            x_center, y_center = g_proj[p].transform_point(g_cen_lon, g_cen_lat, ccrs.Geodetic())
            g_extent[p] = (x_center-g_xdata_span[p], x_center+g_xdata_span[p], 
                           y_center-g_ydata_span[p], y_center+g_ydata_span[p])

def plots_draw(mode):
    global g_proj
    global g_projs
    global g_axis
    global g_view
    global g_cen_lon, g_cen_lat
    global g_dim_x, g_dim_y
    global g_xdata_span, g_ydata_span
    global g_mode

    #print(f"plots_draw: {mode}")

    # Save extent info
    if not mode == "init":
        for p in g_projs:
            x1, x2, y1, y2 = g_axis[p].get_extent()
            g_extent[p] = (x1, x2, y1, y2)
            g_xdata_span[p] = abs(x2-x1)/2
            g_ydata_span[p] = abs(y2-y1)/2
            g_mode = g_axis[p].get_navigate_mode()
            if mode == "setting write grid" and not p == "Orthographic":
                g_extent[p] = g_extent["LambertConformal"]
                g_xdata_span[p] = g_xdata_span["LambertConformal"]
                g_ydata_span[p] = g_ydata_span["LambertConformal"]

    # Remove any old plots
    if g_projs:
        for p in g_projs:
            g_axis[p].remove()

    # Create/recreate projections (recreate needed when the center changes)
    projs_create(mode)

    j = 1; g_axis = {}
    for p in g_projs:

        g_axis[p] = fig.add_subplot(g_dim_x, g_dim_y, j, projection=g_proj[p])
        g_axis[p].add_feature(cfeature.COASTLINE)
        g_axis[p].add_feature(cfeature.BORDERS)
        g_axis[p].add_feature(cfeature.STATES)
        g_axis[p].gridlines()

        # Plot central longitude and latitude
        g_axis[p].plot(*(g_cen_lon, g_cen_lat), transform=ccrs.Geodetic(), marker='*', ms=5, color='red')

        # Set axis view (global or regional/zoomed)
        if g_view[p] == "global":
            g_axis[p].set_title(f"{p} global")
            g_axis[p].set_global()
        else:
            g_view[p] = "regional"
            g_axis[p].set_title(f"{p} regional")
            x_center, y_center = g_proj[p].transform_point(g_cen_lon, g_cen_lat, ccrs.Geodetic())
            g_extent[p] = (x_center-g_xdata_span[p], x_center+g_xdata_span[p], 
                           y_center-g_ydata_span[p], y_center+g_ydata_span[p])

            # Adjust for potential whole-globe boundary condition
            #x1, x2, y1, y2 = g_extent[p]
            #if abs(x1) == abs(x2) and abs(y1) == abs(y2):
            #    print("plots_draw: ADJUSTING EXTENT FOR POTENTIAL BOUNDARY CONDITION")
            #    x1, x2, y1, y2 = g_extent[p]
            #    g_extent[p] = (x1+1, x2-1, y1+1, y2-1)

            g_axis[p].set_extent(g_extent[p], crs=g_proj[p])

        # Restore mode
        g_axis[p].set_navigate_mode(g_mode)
    
        j = j + 1
    if mode == "init":
        set_write_component("LambertConformal")

def draw_box_geodetic(ax):
    x_degrees = 10
    y_degrees = 10
    for p in g_projs:
        lons = []
        lats = []
        #left
        for i in range(0, 33):
            lons.append(g_cen_lon-x_degrees)
            lats.append(g_cen_lat-y_degrees+i*y_degrees/16)
        #top
        for i in range(0, 33):
            lons.append(g_cen_lon-x_degrees+i*x_degrees/16)
            lats.append(g_cen_lat+y_degrees)
        #right
        for i in range(0, 33):
            lons.append(g_cen_lon+x_degrees)
            lats.append(g_cen_lat+y_degrees-i*y_degrees/16)
        #bottom
        for i in range(0, 33):
            lons.append(g_cen_lon+x_degrees-i*x_degrees/16)
            lats.append(g_cen_lat-y_degrees)
    ax.plot(lons, lats, transform=ccrs.Geodetic(), color='green', linewidth=5, alpha=1.0)

def draw_box_xy_data(tgt_index, src_index, color):
    cen_lon, cen_lat, lwr_lon, lwr_lat, upr_lon, upr_lat, xspan, yspan, extent = get_dims(src_index, color)
    x1, x2, y1, y2 = extent 
    # Center circle usng data coords
    g_axis[tgt_index].plot(*((x1+x2)/2, (y1+y2)/2), transform=g_axis[src_index].projection, marker='o', ms=15, color='orange')
    # Star inside of circle using geodetic coords (just a test)
    g_axis[tgt_index].plot(*(cen_lon, cen_lat), transform=ccrs.Geodetic(), marker='*', ms=10, color='blue')
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

    g_axis[tgt_index].plot(xs, ys, transform=g_axis[src_index].projection, color=color, linewidth=2, alpha=1.0, linestyle='solid')

def show_help():
    print("**************************************")
    print("*         Zoom: zoom grid            *")                  
    print("*          Pan: shift grid           *")
    print("*     Spacebar: set write component  *")
    print("*        Key g: set global view      *")
    print("*        Key r: set regional view    *")
    print("*        Key y: output yaml config   *")
    print("*        Key 0: reset grid to CONUS  *")
    print("*        Key h: show this help       *")
    print("* Middle Click: center grid          *")
    print("**************************************")

#
# Main
#

plt.rcParams["figure.raise_window"] = False

fig = plt.figure()
fig.canvas.mpl_connect('button_press_event', on_button_press)
fig.canvas.mpl_connect('key_press_event', on_key_press)
show_help()

plots_draw("init") #1

while True:
    plt.pause(.5)
