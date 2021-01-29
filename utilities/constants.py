#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 10:04:08 2019

  globals used in analyses

@author: jesse
"""


#Water content plus ice content demarking cloud (from toivanen2019) .1g/kg
# "liquid water mixing ratio is the ratio ofthe mass of liquid water to the mass of dry air in a given volume of air. Likewise, the ice mixing ratio is the ratio of the mass of frozen water to the mass of dry air"
cloud_threshold = 0.01 # g/kg cloud boundary of liquid + ice /kg air

# dictionaries for named locations/extents
latlons, extents={},{}

## Waroona extents
# local (first day of burn, escarp covered)
extents['waroona']    = [115.76,116.2, -33.05,-32.7] 
# full fire area + a little bit
extents['waroonaf']    = [115.6,116.21, -33.2,-32.75] 
# synoptic view
extents['waroonas']   = [112,120,-34.5,-31]
# zoomed in view
extents['waroonaz']    = [115.88, 116.19, -32.92,-32.83] # zoom in on fire
## Nests centre: -32.9, 116.1
## Nests resolution: 0.036 384x384, 0.01, 0.0028 (~300m)


latlons['waroona']    = -32.84, 115.93  # suburb centre: -32.8430, 115.8526
latlons['hamel']      = -32.8725, 115.92 # just south of waroona
latlons['yarloop']    = -32.96, 115.90  # suburb centre: -32.9534, 115.9124
latlons['wagerup']    = -32.92, 115.91  # wagerup
latlons['AWS_wagerup']    = -32.92, 115.91  # AWS at wagerup, 40 m asl
latlons['perth']      = -31.9505, 115.8605
latlons['fire_waroona'] = -32.89, 116.17
latlons['fire_waroona_upwind'] = -32.89 -0.004, 116.17+0.009 # ~ 1km from fire

## Extra locs
latlons['sydney']    = -33.8688, 151.2093
latlons['brisbane']  = -27.4698, 153.0251
latlons['canberra']  = -35.2809, 149.1300
latlons['melbourne'] = -37.8136, 144.9631

# two PyroCB
latlons['pyrocb_waroona1'] = -32.87,116.1 # ~4pm first day
latlons['pyrocb_waroona2'] = 0,0 # 1100-1400 second day

# Sir Ivan locations
extents['sirivan']    = [149.2, 150.4, -32.4, -31.6]
extents['sirivan_linescan']   = [149.48, 150.04, -32.18, -31.85]
extents['sirivanz']   = [149.4, 150.19, -32.2, -31.8]
extents['sirivans']   = [147,154, -34, -29] # synoptic
extents['sirivan_pcb']= [149.675, 150.1, -32.08,-31.9]

latlons['dunedoo']    = -32.019, 149.39
latlons['borambil']   = -32.0367, 150.00
latlons['uarbry']     = -32.047280, 149.71
latlons['coolah']     = -31.8234,149.722
latlons['sirivan']    = latlons['uarbry'] # no idea where sir ivan is..
latlons['cassillis']      = -32.01, 150.0
latlons['leadville'] = -32.0383, 149.5779
latlons['merotherie'] = -32.1586, 149.5696
latlons['turill'] = -32.1692, 149.8360
latlons['fire_sirivan'] = -32.05, 149.59
latlons['fire_sirivan_upwind'] = -32.01, 149.5
# one pyrocb

#latlons['pyrocb_sirivan'] = 0,0 # 0530UTC=XXXX local time
# The pycb over Sir Ivan was around 0530 UTC on Sunday 12 Feb (or a bit after).
# The exact location is a bit tricky, because the pycb would be downstream of 
# the actual fire front, but around the location of Uarbry 
# (which was burned over) is probably a good place to start. 
# You'll probably have to move the cross section around a bit to see what 
# lat/longs get the most interesting picture.

extents['NYE'] = [149.2,150.05, -36.5, -35.85]


extents['KI'] = [136.56,137,-36,-35.66]
extents['KIz'] = [136.56,137,-36,-35.66]
latlons['parndana'] = -35.79, 137.262
latlons['fire_KI'] = -35.72, 136.92
