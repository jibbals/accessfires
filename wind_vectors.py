# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 13:32:37 2019
    Look at wind vector fields in 3D
@author: jgreensl
"""

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.ticker as tick
import numpy as np
from datetime import datetime

from utilities import fio, plotting



# Read cloud and wind data

slv, ro1, th1, th2 = fio.read_waroona(dtime= datetime(2016,1,6,9),extent=plotting._extents_['waroona'], add_winds=True)

topog, = slv.extract(['topog'])
u,v = ro1.extract(['u','v'])
w   = th1.extract(['upward_air_velocity'])[0]
zth = th1.extract('z_th')[0] 
# show 1 in 10 lats and lons, and just show specific heights
skip = (slice(None,None,10), slice(None,None,10), np.array([10,28,40,50]))

lon = u.coord('longitude').points
lat = u.coord('latitude').points
lh = u.coord('level_height').points # approximate level height
mlon2,mlat2 = np.meshgrid(lon,lat)
mlon3,mlat3,mlh3 = np.meshgrid(lon,lat,lh)
#for thing in [u,v,w,zth,lon,lat,lh,mlon,mlat,mlh]:
#    print(np.shape(thing))

# Need to rotate so z dim is last:
u3 = np.moveaxis(u[0].data, 0,2)
v3 = np.moveaxis(v[0].data, 0,2)
w3 = np.moveaxis(w[0].data, 0,2)
zth3=np.moveaxis(zth.data,0,2)
mlon3s, mlat3s, mlh3s = mlon3[skip],mlat3[skip],mlh3[skip]
u3s = u3[skip]
v3s = v3[skip]
w3s = w3[skip]
zth3s = zth3[skip]

# w3 will form the colour map
# flatten and normalize
#w3n = (w3s.ravel() - w3s.min()) / w3s.ptp()
w3sr = w3s.ravel()
# repeat for each body line and two head lines
w3src = np.concatenate((w3sr, np.repeat(w3sr,2)))

plt.close() 
fig=plt.figure(figsize=[10,10])
#ax=plt.subplot(1,2,1, projection='3d')
ax = fig.gca(projection='3d')

# Plot topography surface
ax.plot_surface(mlon2, mlat2, topog.data, cmap='terrain', 
                linewidth=0, antialiased=True, vmin=-150,vmax=550)

# Show wind vectors overhead, coloured by vertical velocity
quiv= plt.quiver(mlon3s,mlat3s,zth3s, u3s, v3s, w3s,
                 cmap=plotting._cmaps_['verticalvelocity'],
                 norm=col.SymLogNorm(0.5),
                 length=0.025,normalize=True)
# this seems to set the colours!
quiv.set_array(w3src)
cb=plt.colorbar(format=tick.ScalarFormatter())
cb.set_label('vertical motion (m/s)')

# lower camera elevation
ax.view_init(elev=14)
plt.savefig('figures/waroona/wind_vectors1.png')
#plt.close()

