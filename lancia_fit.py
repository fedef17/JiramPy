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

######################  SETUP  #########################

# Structure that contains the pixels to be considered. Output of leggi_sav2.py.
picfile = '/home/federico/Jiram/DATA/JM0003_SAV_500/nord/aur_JM0003_geo500.pic'

# Output folder where the single spectra will be written in text format in order to be read by the Fortran code.
cart_out = '/home/federico/Jiram/DATA/JM0003_spettri/aur_N_nadir_80/'

# FILTERING ------>

max_emiss_angle = 70 # Maximum emission angle

thres_h3p = 0.01 # Threshold for the H3+ index, calculated by the jirfu.ind_h3p() function. Select thres_h3p = 0 to include all data.

# The pixels selected (with all geometrical information) will be saved in the pickle file all_pixs.pic inside cart_out.

########################################################

pixtot = pickle.load(open(picfile,'r'))
pixs = pixtot.view(np.recarray)

# filtro da applicare agli spettri:
cond = (pixs.emiss_angle < max_emiss_angle) & (pixs.ind_h3p > thres_h3p)

#############################################################################
print('ok')

print(len(pixs[cond]))
n = len(pixs[cond])

pixels = pixs[cond]
pickle.dump(pixels,open(cart_out+'all_pixs.pic','w'))

for i in range(n):
    nome = cart_out+'Spe_'+'{0:0=5d}'.format(i)+'.dat'

    wl1 = pixs[cond].wl[i]
    spe1 = pixs[cond].spe[i]
    mask = jirfu.findspi(wl1,spe1)

    comm = pixs[cond].cubo[i]+' -> line: {:3d}, sample: {:3d}'.format(pixs[cond].line[i],pixs[cond].sample[i])
    jirfu.write_obs_JIR(wl1,spe1,mask,nome,comment=comm)
