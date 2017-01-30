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
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as colors

cart = '/home/fede/Scrivania/Jiram/DATA/JM0003_limb/'
#fi = cart + 'aur_JM0003_limb2000.pic'
fi = cart + 'aur_JM0003_mapelimb.pic'
pixlimb = pickle.load(open(fi,'r'))
pi = pixlimb.view(np.recarray)
col2 = pi.ind_h3p
#condN=(~np.isnan(col2)) & (pi.flag_surf == 0) & (pi.tg_lon < 220) & (pi.tg_lon > 160)
condN=(~np.isnan(col2))
y2 = pi[condN].tg_alt
x2 = pi[condN].tg_lat
col=col2[condN]
ny = 20
nx = 30
xgri, ygri = np.mgrid[0:90:nx*1j, 0:1500:ny*1j]
steplo = xgri[1,0]-xgri[0,0]
stepla = ygri[0,1]-ygri[0,0]
cols = -np.ones((nx,ny))
num = -np.ones((nx,ny))
for i in range(nx):
    for j in range(ny):
        glo = xgri[i,j]
        gla = ygri[i,j]
        cond = (x2-glo >= -steplo/2) & (x2-glo < steplo/2) & (y2-gla >= -stepla/2) & (y2-gla < stepla/2)
        if len(col[cond]) > 0:
            cols[i,j] = np.max(col[cond])
            num[i,j] = len(col[cond])
cols[cols<-0.5]=float('nan')
#pl.contourf(xgri, ygri, cols)
#pl.colorbar()
#pl.scatter(pi.tg_lat,pi.tg_alt,c=pi.ind_h3p,edgecolor='none',vmin=0.0,vmax=0.08,s=5)

# ci = pl.hist2d(pi.ind_h3p[condN],pi.tg_alt[condN],bins=[10,15], range=[[0,0.2],[0,1500]])
# pl.close()
# pl.pcolor(ci[1][:-1],ci[2][:-1],ci[0].T, norm=colors.LogNorm())
# pl.colorbar()
# pl.show()

from mpl_toolkits.mplot3d import Axes3D

pis = pi[(pi.tg_lat > 55) & (pi.tg_alt < 1500) & (pi.tg_lon > 150) & (pi.tg_lon < 250)]
# lat, indh3p = np.meshgrid(pis.tg_lat,pis.ind_h3p)
# lat, lon = np.meshgrid(pis.tg_lat,pis.tg_lon)
# lat, alt = np.meshgrid(pis.tg_lat,pis.tg_alt)

#pdf_3d = PdfPages(cart+'3D_plot.pdf')

fig = pl.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(lat,lon,alt,c=indh3p)

for lat,lon,alt,ind in zip(pis.tg_lat,pis.tg_lon,pis.tg_alt,pis.ind_h3p):
    #pl.scatter(lat,alt,c=ind,edgecolor='none',vmin=0.0,vmax=0.2,s=30)
    pl.scatter(lat,lon,alt,c=ind,edgecolor='none',vmin=0.0,vmax=0.2,s=30)

#pl.show()
#sys.exit()

for angle in range(0, 360, 20):
    ax.view_init(15, angle)
    pl.draw()
    pl.savefig(cart+'plot_3d_'+'{0:03d}'.format(angle)+'.png',format='png')
    pdf_3d.savefig()
#ax.colorbar()
#pl.show()
