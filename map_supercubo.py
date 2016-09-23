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


#cart = '/home/fede/Scrivania/Jiram/DATA/supercubo_ale/Res/'
cart = '/home/fede/Scrivania/Jiram/DATA/supercubo_ale_2/Res/'
#fun=pickle.load(open(cart+'funzmap.pyc','r'))

# ii,col, err_col = sbm.read_res_jir_3(cart+'CD-H3p.dat')
# ii,temp, err_temp = sbm.read_res_jir_3(cart+'VT-H3p.dat')
# ii,chi = sbm.read_res_jir_4(cart+'chisq.dat')
ii,jj,col, err_col = sbm.read_res_jir(cart+'CD-H3p.dat')
ii,jj,temp, err_temp = sbm.read_res_jir(cart+'VT-H3p.dat')
ii,jj,chi = sbm.read_res_jir_2(cart+'chisq.dat')
ii,jj,off = sbm.read_res_jir_2(cart+'offset_ok.dat')

cols = np.zeros((101,101))
err_cols = np.zeros((101,101))
temps = np.zeros((101,101))
err_temps = np.zeros((101,101))
chis = np.zeros((101,101))

# for l in ii:
#     print(fun[l-1001])
#     cols[fun[l-1001][0],fun[l-1001][1]] = col[l-1001]
#     err_cols[fun[l-1001][0],fun[l-1001][1]] = err_col[l-1001]
#     temps[fun[l-1001][0],fun[l-1001][1]] = temp[l-1001]
#     err_temps[fun[l-1001][0],fun[l-1001][1]] = err_temp[l-1001]
#     chis[fun[l-1001][0],fun[l-1001][1]] = chi[l-1001]
for i,j,co,te,erco,erte,chi in zip(ii,jj,col,temp,err_col,err_temp,chi):
    i=i-100
    j=j-100
    cols[i,j]=co
    temps[i,j]=te
    err_temps[i,j]=erte
    err_cols[i,j]=erco
    chis[i,j]=chi

cond = (cols < 0.0) | (cols > 1e13) | (chis > 5)
cols[cond]=float('nan')
temps[cond]=float('nan')
chis[cond]=float('nan')
print(cond)

fig = pl.figure(figsize=(8, 6), dpi=150)
pl.title('H3+ column [cm-2]')
pl.imshow(cols, cmap='jet', interpolation="nearest")
pl.colorbar(orientation="horizontal")
fig.savefig(cart+'col_h3p.eps', format='eps', dpi=150)
pl.close()
#pl.show()
fig = pl.figure(figsize=(8, 6), dpi=150)
pl.title('H3+ temperature [K]')
pl.imshow(temps, cmap='jet', interpolation="nearest")
pl.colorbar(orientation="horizontal")
fig.savefig(cart+'temp_h3p.eps', format='eps', dpi=150)
pl.close()
#pl.show()
pl.imshow(chis, cmap='jet', interpolation="nearest")
pl.colorbar()
#pl.show()
