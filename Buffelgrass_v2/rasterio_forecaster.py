import rasterio as ras
from datetime import date
from datetime import timedelta
import numpy as np
import matplotlib.pyplot as plt
#%%

def buffer_ppt_vals(data):
    #data = list(data)
    events = 0
    low_ppt_days = 4
    for ppt in data:
        if ppt == -9999:
            return -9999
        if ppt >= 6.35 and low_ppt_days >= 3:
            events += 1
            low_ppt_days = 0
        elif ppt == 0:
            low_ppt_days += 1
        elif ppt < 6.35 and low_ppt_days != 0:
            low_ppt_days += 1
        else:
            low_ppt_days = 0
    if events > 4:
        events = 4
    return events


today = date(2023, 10, 22)
i = 29
raster_stack = np.zeros((137,144,30))
while i >= 0:
    delay = timedelta(days = i)
    lookback_day = today - delay
    lookback_day_string = lookback_day.strftime('%Y%m%d')
    #/geo-data/climate_data/prism/prism_data/prism_ppt
    path1 = 'prism-precip-files/PRISM_ppt_early_4kmD2_'+lookback_day_string+'_bil.tif'
    path2 = 'prism-precip-files/PRISM_ppt_provisional_4kmD2_'+lookback_day_string+'_bil.tif'
    try: dataset = ras.open(path1)
    except: 
        dataset = ras.open(path2)
        print('used path2')
    raster_stack[:,:,i] = dataset.read(1)[310:447,241:385]
    i -= 1
       

forecast = np.apply_along_axis(buffer_ppt_vals, 2, raster_stack)
forecast = np.ma.masked_values(forecast, -9999)
dt_tf = dataset.transform
cropped_transform = ras.Affine(dt_tf[0], dt_tf[1], dt_tf[2]+241*dt_tf[0],
                               dt_tf[3], dt_tf[4], dt_tf[5]+310*dt_tf[4])


new_tif = ras.open(
    'test.tif',
    'w',
    height = forecast.shape[0],
    width = forecast.shape[1],
    count = 1,
    dtype = 'int32',
    crs = dataset.crs,
    transform = cropped_transform)

new_tif.write(forecast, 1)
new_tif.close()
#%%
from_http = ras.open('C:/Users/twmat/USA_NPN/Buffelgrass_v2/test.tif')
fig, ax = plt.subplots()
im = ax.imshow(from_http.read(1), cmap='magma', vmin=-2, vmax=5)