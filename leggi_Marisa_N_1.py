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
from mayavi import mlab
from mpl_toolkits.basemap import Basemap

cart = '/home/fede/Scrivania/Jiram/DATA/TEST_marisa/'
cart2 = '/home/fede/Scrivania/Jiram/DATA/TEST_marisa/Spe_prova/'

pixtot = pickle.load(open(cart+'pix_nadir_N.pic','r'))
pixs = pixtot.view(np.recarray)

cond = (pixs.emiss_angle < 80) & (pixs.ind_h3p > 0.015)
pl.scatter(pixs[cond].pc_lon,pixs[cond].pc_lat,c=pixs[cond].ind_h3p)
pl.colorbar()
pl.show()
pl.close()
#pl.plot(pixs[cond].wl[13], pixs[cond].spe[13])
print(pixs[cond].ind_h3p[13])
#pl.show()
print('ok')

print(len(pixs[cond]))
n = len(pixs[cond])

# for i in range(n):
#     nome = cart2+'Spe_'+'{0:0=4d}'.format(i)+'.dat'
#     print(nome)
#     err = jirfu.findspi(wls,spe)
#     jirfu.write_obs_JIR(pixs[cond].wl[i],pixs[cond].spe[i],nome,comment='emiss < 80 and ind_h3p > 0.015')
#     # call(cart+'./fomi_mipas')
#     # nomeout = cart+'output__mipas.dat'
#     # alt_fomi, cr_fomi = sbm.leggioutfomi(nomeout)

cartres = '/home/fede/Scrivania/Jiram/DATA/TEST_marisa/Res_prova/'

ii,col, err_col = jirfu.read_res_jir_3(cartres+'CD-H3p.dat')
ii,temp, err_temp = jirfu.read_res_jir_3(cartres+'VT-H3p.dat')
ii,chi = jirfu.read_res_jir_4(cartres+'chisq.dat')
ii,off = jirfu.read_res_jir_4(cartres+'offset_ok.dat')

cond2 = ((chi > 4) | (col < 0.0) | (off < 0.0))
print(len(col))
print(type(col))
print(len(col[cond2]))
col[cond2] = float('NaN')
temp[cond2] = float('NaN')
chi[cond2] = float('NaN')
off[cond2] = float('NaN')

print(pixs[cond].emiss_angle)

col_c = col*np.cos(np.pi*pixs[cond].emiss_angle/180.0)

fig = pl.figure(figsize=(8, 6), dpi=150)
pl.title('H3+ column [cm-2]')
pl.scatter(pixs[cond].pc_lon,pixs[cond].pc_lat,c=col_c)
#pl.show()
pl.colorbar(orientation="horizontal")
fig.savefig(cart+'col_h3p.eps', format='eps', dpi=150)
pl.close()

fig = pl.figure(figsize=(8, 6), dpi=150)
pl.title('H3+ temp [K]')
pl.scatter(pixs[cond].pc_lon,pixs[cond].pc_lat,c=temp)
#pl.show()
pl.colorbar(orientation="horizontal")
fig.savefig(cart+'temp_h3p.eps', format='eps', dpi=150)
pl.close()

fig = pl.figure(figsize=(8, 6), dpi=150)
pl.title('chi squared')
pl.scatter(pixs[cond].pc_lon,pixs[cond].pc_lat,c=chi)
#pl.show()
pl.colorbar(orientation="horizontal")
fig.savefig(cart+'chi.eps', format='eps', dpi=150)
pl.close()

fig = pl.figure(figsize=(8, 6), dpi=150)
pl.title('Offset [W/m2/nm/sr]')
pl.scatter(pixs[cond].pc_lon,pixs[cond].pc_lat,c=off)
#pl.show()
pl.colorbar(orientation="horizontal")
fig.savefig(cart+'offset.eps', format='eps', dpi=150)
pl.close()

fig = pl.figure(figsize=(8, 6), dpi=150)
map = Basemap(projection='npstere',boundinglat=60,lon_0=270,resolution='l')
map.drawmeridians(np.arange(0,360,30))
map.drawparallels(np.arange(60,90,10))
x, y = map(pixs[cond].pg_lon,pixs[cond].pg_lat)
sca = map.scatter(x,y,c = col_c,vmin=0.2e12,vmax=2.9e12)
pl.colorbar(orientation="horizontal")

fig2 = pl.figure(figsize=(8, 6), dpi=150)
map = Basemap(projection='npstere',boundinglat=60,lon_0=270,resolution='l')
map.drawmeridians(np.arange(0,360,30))
map.drawparallels(np.arange(60,90,10))
x, y = map(pixs[cond].pg_lon,pixs[cond].pg_lat)
sca = map.scatter(x,y,c = temp,vmin=0.2e12,vmax=2.9e12)
pl.colorbar(orientation="horizontal")

fig3 = pl.figure(figsize=(8, 6), dpi=150)
map = Basemap(projection='npstere',boundinglat=60,lon_0=270,resolution='l')
map.drawmeridians(np.arange(0,360,30))
map.drawparallels(np.arange(60,90,10))
x, y = map(pixs[cond].pg_lon,pixs[cond].pg_lat)
sca = map.scatter(x,y,c = pixs[cond].emiss_angle,vmin=0.2e12,vmax=2.9e12)
pl.colorbar(orientation="horizontal")
pl.show()