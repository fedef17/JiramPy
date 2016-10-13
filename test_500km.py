#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os.path
import matplotlib.pyplot as pl
import math as mt
import spect_base_module as sbm
import jirfu
from subprocess import call
import pickle
import scipy.io as io
from scipy.interpolate import PchipInterpolator as spline
#from mayavi import mlab
from mpl_toolkits.basemap import Basemap

nomi = ['line','sample','pc_lon','pc_lat','pg_lon','pg_lat','incid','emiss','pc_lon_500','pc_lat_500','pg_lon_500','pg_lat_500','incid_500','emiss_500']

tipi2 = 2*'i4,'+11*'f8,'+'f8'
geo_ = np.empty(1, dtype=tipi2)
geo_.dtype.names = nomi


cart = '/home/fede/Scrivania/Jiram/DATA/GEO_500/'

cubo = io.readsav(cart+'GEO_500.sav')
geo500_cubo = cubo.geo500
geo_cubo = cubo.geo

n_lin = 11
n_sam = 256

print(np.shape(geo_cubo))
print(geo_cubo[0,:,:])
geo500 = geo_
for i in range(n_lin):
    for j in range(n_sam):
        geol=geo_
        el = 0
        geol[0][el] = i
        el+=1
        geol[0][el] = j
        el+=1
        for l in range(6):
            geol[0][el] = geo_cubo[l,i,j]
            el+=1
        for l in range(6):
            geol[0][el] = geo500_cubo[l,i,j]
            el+=1
        geo500 = np.append(geo500,geol)

geo500 = geo500[1:]

geo500 = geo500.view(np.recarray)
print(geo500.pc_lon)

geo500 = geo500[(geo500.emiss < 75)]

jirfu.stereoplot(geo500.pc_lon,geo500.pc_lat,geo500.emiss,cart+'emiss.eps',title='Emission angle',polo='N',show=True)

fig = pl.figure(figsize=(8, 6), dpi=150)
pl.title('TEST 500 km')
polo = 'N'

if polo == 'N':
    map = Basemap(projection='npstere',boundinglat=60,lon_0=180,resolution='l')
    map.drawparallels(np.arange(60,90,10))
else:
    map = Basemap(projection='spstere',boundinglat=-60,lon_0=180,resolution='l')
    map.drawparallels(np.arange(-80,-50,10))

aur_lon_0,aur_lat,aur_lon,aur_theta = jirfu.leggi_map_aur(polo)

map.drawmeridians(np.arange(0,360,30),labels=[1,1,0,1],fontsize=10)

x, y = map(geo500.pc_lon,geo500.pc_lat)
sca = map.scatter(x,y,color = 'green',s=5,label='GEO 0 km')

x, y = map(geo500.pc_lon_500,geo500.pc_lat_500)
sca = map.scatter(x,y,color = 'red',s=5,label='GEO 500 km')

pl.legend()

x, y = map(360-aur_lon,aur_lat)
#map.scatter(x,y,color = 'white',edgecolors='black',s=15,marker = 'o')
x = np.append(x,x[0])
y = np.append(y,y[0])
pl.plot(x,y,color='black',linewidth=2.0,linestyle='--')
pl.show()
fig.savefig(cart + 'TEST_500.eps', format='eps', dpi=150)
pl.close()

