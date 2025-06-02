#
# SCRAP
#

#g_fig.canvas.mpl_connect('resize_event', on_resize)
#def on_resize(event):
#    print(f"Figure resized: width={event.width}, height={event.height}")
#def on_xlim_changed(axes):
#    print("X-axis limits changed:", axes.get_xlim())
#def on_ylim_changed(axes):
#    print("Y-axis limits changed:", axes.get_ylim())
#uds.axis[p].callbacks.connect('xlim_changed', on_xlim_changed)
#uds.axis[p].callbacks.connect('ylim_changed', on_ylim_changed)

# Old filtering code
#    if False:
#        lat_span = 30
#        lon_span = 30
#        lat_min, lat_max = round(uds.cen_lat-lat_span/2), \
#                           round(uds.cen_lat+lat_span/2)
#        dim0 = (lat_max-lat_min)*4+1
#        lon_min, lon_max = round(360+uds.cen_lon-lon_span/2), \
#                           round(360+uds.cen_lon+lon_span/2)
#        dim1 = (lon_max-lon_min)*4+1
#        mask = (lats >= lat_min) & (lats <= lat_max) & (lons >= lon_min) & (lons <= lon_max)
#        filtered_lats = lats[mask]
#        filtered_lons = lons[mask]
#        filtered_zdata = zdata[mask]
#        filtered_vdata = vdata[mask]
#        lats = np.array(filtered_lats)
#        lons = np.array(filtered_lons)
#        zdata = np.array(filtered_zdata)
#        vdata = np.array(filtered_vdata)
#        lats = np.reshape(lats, (dim0, dim1))
#        lons = np.reshape(lons, (dim0, dim1))
#        zdata = np.reshape(zdata, (dim0, dim1))
#        vdata = np.reshape(vdata, (dim0, dim1))
