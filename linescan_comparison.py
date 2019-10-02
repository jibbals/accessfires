# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 12:16:43 2019
    Compare model output to sirivan line scans
@author: jgreensl
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from glob import glob
from datetime import datetime
from os.path import basename

from utilities import fio, utils, plotting

# sirivan linescan figures extent
linescan_extent =  [149.48, 150.04, -32.18, -31.85]

# loop over linescans
files = glob('data/linescans/SIRIVAN_*.jpg')
files.sort()

for file in files:
    fname=basename(file)
    dtime = datetime.strptime(fname,"SIRIVAN_%Y%m%d %H%M.jpg")
    print(dtime)

#show linescans on row 1
img1 = mpimg.imread('data/linescans/SIRIVAN_20170211 1525.jpg')
plt.close()
fig = plt.figure(figsize=(10,14))
ax1 = fig.add_axes([0,.5,1,.5], frameon=False)
plt.imshow(img1)

ax1.axes.get_yaxis().set_visible(False)
ax1.axes.get_xaxis().set_visible(False)

#show fire plan on row 2 (matching linescan times)
# can show either googlemap (roads, rivers), or satellite map (true colour)
google=True
if google:
    _,ax2,mproj = plotting.map_google(linescan_extent, fig=fig, zoom=12,
                                      subplot_extent=[0,0.05,1,0.45],
                                      draw_gridlines=False)
else:
    _,ax2,mproj,dproj = plotting.map_satellite(linescan_extent, fig=fig,
                                               subplot_extent=[0,0.05,1,0.45])
fio.save_fig('sirivan_run1',_sn_,"linescan_%s"%dstr,plt,dpi=400)
