#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os.path
import matplotlib.pyplot as pl
import math as mt
import spect_base_module as sbm
from subprocess import call
import pickle
import scipy.io as io
from scipy.interpolate import PchipInterpolator as spline


lin1 = 20
lin2 = 55
sam1 = 50
sam2 = 240

cart = '/home/fede/Scrivania/Jiram/ANALISI/TEST_1/'

cubo = io.readsav(cart+'cubo_test.sav')
data = cubo.spe

print(type(data))

wew = io.readsav(cart+'wls.sav')
wl = wew.wls

data[(data < 0)] = np.nan

fig = pl.figure(figsize=(8, 6), dpi=150)
pl.title('Test data')
pl.imshow(data[171,20:55,50:240], cmap='jet', interpolation="nearest")
pl.colorbar(orientation="horizontal")
#pl.show()

fig.savefig(cart+'DATA.eps', format='eps', dpi=150)
pl.close()

fig = pl.figure(figsize=(8, 6), dpi=150)
pl.title('Spettro giorno')
pl.plot(wl, data[:,50,128])
#pl.show()
pl.ylabel('Radiance (W m$^{-2}$ nm$^{-1}$ sr$^{-1}$)')
pl.xlabel('Wavelength (nm)')
fig.savefig(cart+'SPE_tot.eps', format='eps', dpi=150)
pl.close()


#sys.exit()

# mappa = np.zeros([lin2-lin1+1,sam2-sam1+1])
# mappa_err = np.zeros([lin2-lin1+1,sam2-sam1+1])
#
# for el1, el2, t, et in zip(i,j,temp,err_t):
#     el11=el1-lin1
#     el22=el2-sam1
#     mappa[el11,el22]=t
#     mappa_err[el11,el22]=et

i,j,temp, err_t = sbm.read_res_jir(cart+'CD-H3p.dat')

#print(i,j,temp,err_t)


#lin1 = min(i)
#lin2 = max(i)
#sam1 = min(j)
#sam2 = max(j)

mappa = np.zeros([lin2-lin1+1,sam2-sam1+1])
mappa_err = np.zeros([lin2-lin1+1,sam2-sam1+1])

mappa[(mappa == 0.0)] = np.nan

for el1, el2, t, et in zip(i,j,temp,err_t):
    el11=el1-lin1
    el22=el2-sam1
    mappa[el11,el22]=t
    mappa_err[el11,el22]=et


mappa[(mappa < 0.0)] = np.nan
mappa[(mappa > 1e13)] = np.nan


i,j,temp, err_t = sbm.read_res_jir(cart+'VT-H3p.dat')

#print(i,j,temp,err_t)


# lin1 = min(i)
# lin2 = max(i)
# sam1 = min(j)
# sam2 = max(j)

mappa_T = np.zeros([lin2-lin1+1,sam2-sam1+1])
mappa_err_T = np.zeros([lin2-lin1+1,sam2-sam1+1])

mappa_T[(mappa == 0.0)] = np.nan

for el1, el2, t, et in zip(i,j,temp,err_t):
    el11=el1-lin1
    el22=el2-sam1
    mappa_T[el11,el22]=t
    mappa_err_T[el11,el22]=et


#mappa[(mappa < 0.0)] = np.nan


i,j,temp = sbm.read_res_jir_2(cart+'chisq.dat')

#print(i,j,temp,err_t)


# lin1 = min(i)
# lin2 = max(i)
# sam1 = min(j)
# sam2 = max(j)

mappa_chi = np.zeros([lin2-lin1+1,sam2-sam1+1])

for el1, el2, t in zip(i,j,temp):
    el11=el1-lin1
    el22=el2-sam1
    mappa_chi[el11,el22]=t


mappa[(mappa < 0.0)] = np.nan
mappa[(mappa > 1e13)] = np.nan
mappa[(mappa < 1.5e12)] = np.nan
mappa[(mappa_chi > 3)] = np.nan

#mappa_T[(mappa_err_T > 100)] = np.nan
#mappa_T[(mappa_T > 1300)] = np.nan
#mappa_T[(mappa < 1.5e12)] = np.nan
#mappa_T[(mappa_chi > 3)] = np.nan

mappa_T[(np.isnan(mappa))] = np.nan


fig = pl.figure(figsize=(8, 6), dpi=150)
pl.title('H3+ column')
pl.imshow(mappa, cmap='jet', interpolation="nearest")
pl.colorbar(orientation="horizontal")

#pl.show()
fig.savefig(cart+'H3p.eps', format='eps', dpi=150)
#pl.close()


fig2 = pl.figure(figsize=(8, 6), dpi=150)
pl.title('H3+ temperature')
pl.imshow(mappa_T, cmap='jet', interpolation="nearest")
pl.colorbar(orientation="horizontal")

#pl.show()
fig2.savefig(cart+'Temp_H3p.eps', format='eps', dpi=150)
pl.close()


# fig3 = pl.figure(figsize=(8, 6), dpi=150)
# pl.title('H3+ temperature')
# pl.imshow(mappa_chi, cmap='jet')
# pl.colorbar()
#
# pl.show()
# fig3.savefig(cart+'Temp_H3p.eps', format='eps', dpi=150)
# pl.close()