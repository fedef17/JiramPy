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

cart = '/home/fede/Scrivania/Jiram/DATA/TEST_marisa/'
cart2 = '/home/fede/Scrivania/Jiram/DATA/TEST_marisa/Spe_prova_N_bis/'

pixtot = pickle.load(open(cart+'pix_nadir_N_bis.pic','r'))
pixs = pixtot.view(np.recarray)

cond = (pixs.emiss_angle < 75) & (pixs.ind_h3p > 0.01)

print('ok')

print(len(pixs[cond]))
n = len(pixs[cond])

pixels = pixs[cond]
pickle.dump(pixels,open(cart2+'all_pixs.pic','w'))

for i in range(n):
    nome = cart2+'Spe_'+'{0:0=4d}'.format(i)+'.dat'
    print(nome)
    wl1 = pixs[cond].wl[i]
    spe1 = pixs[cond].spe[i]
    mask = jirfu.findspi(wl1,spe1)
    print(type(pixs[cond].line[i]))
    comm = pixs[cond].cubo+' -> line: {:3d}, sample: {:3d}'.format(pixs[cond].line[i],pixs[cond].sample[i])
    jirfu.write_obs_JIR(wl1,spe1,mask,nome,comment=comm)
    # call(cart+'./fomi_mipas')
    # nomeout = cart+'output__mipas.dat'
    # alt_fomi, cr_fomi = sbm.leggioutfomi(nomeout)

sys.exit()

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

pl.scatter(pixs[cond].pc_lon,pixs[cond].pc_lat,c=col_c)
pl.colorbar()
pl.show()

pl.scatter(pixs[cond].pc_lon,pixs[cond].pc_lat,c=temp)
pl.colorbar()
pl.show()

pl.scatter(pixs[cond].pc_lon,pixs[cond].pc_lat,c=chi)
pl.colorbar()
pl.show()

pl.scatter(pixs[cond].pc_lon,pixs[cond].pc_lat,c=off)
pl.colorbar()
pl.show()