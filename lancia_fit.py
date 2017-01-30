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

#cart = '/home/fede/Scrivania/Jiram/DATA/TEST_marisa/'  # cartella di run

# struttura da cui prendo gli spettri :
picc = '/home/federico/Jiram/DATA/JM0003_SAV_500/nord/aur_JM0003_geo500.pic'

# Dove salvo gli spettri :
#cart2 = '/home/fede/Scrivania/Jiram/DATA/TEST_marisa/Spe_prova_S_bis/'
cart2 = '/home/federico/Jiram/DATA/JM0003_spettri/aur_N_nadir_80/'

pixtot = pickle.load(open(picc,'r'))
pixs = pixtot.view(np.recarray)

# filtro da applicare agli spettri:
cond = (pixs.emiss_angle > 70) & (pixs.emiss_angle < 80) & (pixs.ind_h3p > 0.01)

#############################################################################
print('ok')

print(len(pixs[cond]))
n = len(pixs[cond])

pixels = pixs[cond]
pickle.dump(pixels,open(cart2+'all_pixs.pic','w'))

for i in range(n):
    nome = cart2+'Spe_'+'{0:0=5d}'.format(i)+'.dat'
    #print(nome)
    wl1 = pixs[cond].wl[i]
    spe1 = pixs[cond].spe[i]
    mask = jirfu.findspi(wl1,spe1)
    #print(type(pixs[cond].line[i]))
    comm = pixs[cond].cubo[i]+' -> line: {:3d}, sample: {:3d}'.format(pixs[cond].line[i],pixs[cond].sample[i])
    jirfu.write_obs_JIR(wl1,spe1,mask,nome,comment=comm)
    # call(cart+'./fomi_mipas')
    # nomeout = cart+'output__mipas.dat'
    # alt_fomi, cr_fomi = sbm.leggioutfomi(nomeout)
