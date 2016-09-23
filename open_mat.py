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

cart = '/home/fede/Scrivania/Jiram/DATA/supercubo_ale_60/'
cart2 = '/home/fede/Scrivania/Jiram/DATA/supercubo_ale_60/Spe/'

mat1 = io.loadmat(cart+'cub.mat')

wew = io.readsav(cart+'../wls.sav')
wl = wew.wls

rad = mat1['H']
x = mat1['X']
y = mat1['Y']
#em = mat1['EM']

pl.imshow(rad[171,:,:],interpolation="nearest")
#pl.imshow(em,interpolation="nearest")
pl.colorbar()
pl.show()
pl.close()


# for i in range(101):
#     for j in range(101):
#         if(np.isnan(rad[0,i,j])): continue
#         if(np.max(rad[150:200,i,j])>0.008): continue
#         if(rad[171,i,j] > 0.002): pl.plot(wl, rad[:,i,j])
pl.plot(wl, rad[:,39,39])
pl.plot(wl, rad[:,57,14])
pl.show()


# fondo =
#
# linee = np.zeros(4)
# wl1 = []
# for lin in linee:
#     linee = sbm.somma_lin(wl1,wl2,spe,fondo)


sys.exit()



a = np.reshape(rad, (336,-1))

mea = np.zeros(336)
std = np.zeros(336)

n=0
for i in range(101):
    for j in range(101):
        if(np.isnan(rad[0,i,j])): continue
        mea = mea + rad[:,i,j]
        n+=1

mea = mea/n

for i in range(101):
    for j in range(101):
        if(np.isnan(rad[0,i,j])): continue
        std = std + (rad[:,i,j]-mea)**2

std = np.sqrt(std/n)


print(n)

fig = pl.figure(figsize=(8, 6), dpi=150)
# for i in range(101):
#     for j in range(101):
#         if(np.isnan(rad[0,i,j])): continue
#         pl.plot(wl,rad[:,i,j])
pl.plot(wl,mea)
pl.plot(wl,std)
pl.show()

l=0

fun = [1,1]
#fi = open(cart2+'nomi.dat', 'w')
for i in range(101):
    for j in range(101):
        if(np.isnan(rad[0,i,j])): continue
        if(np.mean(rad[:,i,j])<1e-5): continue
        #name = cart2 + 'Spe_'+'{:3d}'.format(100+i)+'_'+'{:3d}'.format(100+j)+'.dat'
        #l+=1
        #fun = np.vstack([fun,[i,j]])
        #fi.write((3*'{:5d}').format(1000+l,i,j)+'\n')
        #name = cart2 + 'Spe_'+'{:4d}'.format(1000+l)+'.dat'
        name = cart2 + 'Spe_'+'{:3d}'.format(i+100)+'_'+'{:3d}'.format(100+j)+'.dat'
        spe = rad[:,i,j]
        sbm.write_obs_JIR(wl,spe,name,comment='Da cubo in '+cart+'cub.mat\n')

#fi.close()
# fun = fun[1:]
# print(fun)
# pickle.dump(fun,open(cart2+'funzmap.pyc','w'))