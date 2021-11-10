#!/bin/env python


from netCDF4 import Dataset         
import matplotlib.pyplot as plt     
from matplotlib.cm import get_cmap  
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim, cartopy_ylim, latlon_coords, ALL_TIMES)   # Import all the necessary functions from wrf-python module to deal with data, map projection, lats & lons, plot bounds etc.
import cartopy.crs as ccrs          
import cartopy.feature as cfeature  
from cartopy.feature import ShapelyFeature  
from cartopy.io.shapereader import Reader   
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER  
import os                           
import math
import numpy as np


wrf_out_path = os.path.join('d1')


wrf_rd = Dataset(wrf_out_path)


slp = getvar(wrf_rd, "slp", timeidx=ALL_TIMES)
wspd = getvar(wrf_rd, "wspd10", timeidx=ALL_TIMES,units="kt")



lats, lons = latlon_coords(slp)

#Func to return index of lat long values
def ret_index(var,lat_1,lon_1,lat_2,lon_2):
    v_lats,v_lons = latlon_coords(var)
    a = abs(v_lats-lat_1)+abs(v_lons-lon_1)
    i,j = np.unravel_index(a.argmin(),a.shape)
    b = abs(v_lats-lat_2)+abs(v_lons-lon_2)
    k,l = np.unravel_index(b.argmin(),b.shape)
    return i,j,k,l




cart_proj = get_cartopy(slp)


smooth_slp = smooth2d(slp, 3, cenweight=8)
ax = plt.axes(projection=ccrs.Mercator())

shapefile_path = os.path.join("shape_files", "world", "world_shape_modified.shp")


world_shape_file_features = ShapelyFeature(Reader(shapefile_path).geometries(), crs=ccrs.PlateCarree(), facecolor='none')

ax.add_feature(world_shape_file_features, linewidth=0.7, edgecolor='black')




ax.set_extent([50, 95, 0, 30], crs=ccrs.PlateCarree())



ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)



gl = ax.gridlines(draw_labels=True, color='gray', alpha=0.5, linestyle='--')


gl.top_labels = False
gl.right_labels = False


gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER



lon_start=55
lon_end=74
lat_start=10
lat_end=28
#bound
#call_indexfunc
lon_li,lat_li,lon_ui,lat_ui=ret_index(smooth_slp,lat_start,lon_start,lat_end,lon_end)

#outerbound
slp_bound=smooth_slp.sel(south_north=slice(lon_li,lon_ui),west_east=slice(lat_li,lat_ui))
wspd_n=wspd.sel(south_north=slice(lon_li,lon_ui),west_east=slice(lat_li,lat_ui))




ilon_start=70
ilon_end=72
ilat_start=10
ilat_end=11
start_t=4
end_t=len(smooth_slp.Time.values)
bound_incr=0.5
dist_travelled=100
app_counter=1




min_a=[]
latt_a=[]
lonn_a=[]
max_wind_list=[]

#########################Func to find distance bet two lat/lon
def find_dist(lat1,lon1,lat2,lon2):
    r=6373.0
    lat1=math.radians(lat1)
    lon1=math.radians(lon1)
    lat2=math.radians(lat2)
    lon2=math.radians(lon2)
    dlon=lon2-lon1
    dlat=lat2-lat1
    a=math.sin(dlat/2)**2+math.cos(lat1)*math.cos(lat2)*math.sin(dlon/2)**2
    c=2*math.atan2(math.sqrt(a),math.sqrt(1-a))
    distance=r*c
    return distance




for idx_time, ele_time in enumerate(slp_bound.Time.values):
    fcst_time = ele_time.astype('datetime64[h]')
    slp_t = slp_bound[idx_time,:,:]
    wspd_t=wspd[idx_time]
    
    
    
    slp_t=slp_t.where((slp_t.XLAT>=ilat_start)&(slp_t.XLAT<=ilat_end)&(slp_t.XLONG>= ilon_start) & (slp_t.XLONG<=ilon_end))
    max_wind=wspd_t.max()
    mini=slp_t.min()


    latt=slp_t.where(slp_t!=slp_t.min().data,slp_t.XLAT).min()

    lonn=slp_t.where(slp_t!=slp_t.min().data,slp_t.XLONG).min()
    
      
    
    if(idx_time!=0):
        arr_l=len(min_a)
        
        
        if(find_dist(latt_a[arr_l-1],lonn_a[arr_l-1],latt,lonn)<(dist_travelled*app_counter)):
            
            app_counter=1
            max_wind_list.append(max_wind)
            min_a.append(mini)
            latt_a.append(latt)
            lonn_a.append(lonn)
        else:
            app_counter=app_counter+1
    else:
        
        app_counter=1
        max_wind_list.append(max_wind)
        min_a.append(mini)
        latt_a.append(latt)
        lonn_a.append(lonn)  

    
    if(ilat_start<=slp_t.XLAT.min()):
        ilat_start=slp_t.XLAT.min()
    else:
        ilat_start=ilat_start-bound_incr

    if(ilat_end>=slp_t.XLAT.max()):
        ilat_end=slp_t.XLAT.max()
    else:
        ilat_end=ilat_end+bound_incr


    if(ilon_start<=slp_t.XLONG.min()):
        ilon_start=slp_t.XLONG.min()
    else:
        ilon_start=ilon_start-bound_incr

    if(ilon_end>=slp_t.XLONG.max()):
        ilon_end=slp_t.XLONG.max()
    else:
        ilon_end=ilon_end+bound_incr

  
###############################################################
print(*max_wind_list)
clr="blue"
cs=0
scs=0
vscs=0
escs=0
sucs=0
for l in range(len(latt_a)-1):
    if(l>=start_t and l<=end_t):
        if(max_wind_list[l]>=34 and max_wind_list[l-1]<34 and l>0 and cs==0):
            scs=0
            cs=1
            vscs=0
            escs=0
            sucs=0
            print("####################")
            print(max_wind_list[l])
            print(l)
            plt.plot(lonn_a[l],latt_a[l],'o',transform=ccrs.PlateCarree())
            plt.text(lonn_a[l],latt_a[l],'CS',color="red",transform=ccrs.PlateCarree())

        if(max_wind_list[l]>=48 and max_wind_list[l-1]<48 and l!=0 and scs==0):
            scs=1
            cs=0
            vscs=0
            escs=0
            sucs=0
            plt.plot(lonn_a[l],latt_a[l],'o',transform=ccrs.PlateCarree())
            plt.text(lonn_a[l],latt_a[l],'SCS',color="red",transform=ccrs.PlateCarree())

        if(max_wind_list[l]>=64 and max_wind_list[l-1]<64 and l!=0 and vscs==0):
            vscs=1
            escs=0
            scs=0
            cs=0
            sucs=0
            plt.plot(lonn_a[l],latt_a[l],'o',transform=ccrs.PlateCarree())
            plt.text(lonn_a[l],latt_a[l],'VSCS',color="red",transform=ccrs.PlateCarree())

        if(max_wind_list[l]>=90 and max_wind_list[l-1]<90 and l!=0 and escs==0):
            escs=1
            scs=0
            cs=0
            vscs=0
            sucs=0
            plt.plot(lonn_a[l],latt_a[l],'o',transform=ccrs.PlateCarree())
            plt.text(lonn_a[l],latt_a[l],'ESCS',color="red",transform=ccrs.PlateCarree())
        if(max_wind_list[l]>=120 and max_wind_list[l-1]<120 and l!=0 and sucs==0):
            escs=0
            scs=0
            cs=0
            vscs=0
            sucs=1
            plt.plot(lonn_a[l],latt_a[l],'o',transform=ccrs.PlateCarree())
            plt.text(lonn_a[l],latt_a[l],'Super CS',color="red",transform=ccrs.PlateCarree())
        plt.plot([lonn_a[l],lonn_a[l+1]],[latt_a[l],latt_a[l+1]],color=clr,transform=ccrs.PlateCarree())



#################################################

ax.set_title("Cyclone Track", fontdict={'fontsize': 20, 'fontweight': 'medium'})
plt.legend(loc="lower right")

plt.tight_layout()





plt.savefig('Track.png')            # Relative or current path


plt.show()
