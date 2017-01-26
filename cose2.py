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

# cart = '/home/fede/Scrivania/Dotto/AbstrArt/Titan_workshop/Data/All_data/'
# cubo1 = io.readsav(cart+'PIXs_VIMS_4-5mu_night2.sav')
# cubo2 = io.readsav(cart+'PIXs_VIMS_4-5mu_night_far.sav')

cart2 = '/home/fede/Scrivania/Dotto/AbstrArt/CH4_HCN_climatology/DATA/'
cubo_ch4 = io.readsav(cart2+'PIXs_HCN-CH4-C2H2_season_sza80.sav')

cu = cubo_ch4.comppix
cu = cu[1:]

names = cu.dtype.names
names = names[:-2]

namesm = ['cubo', 'year', 'dist', 'lat', 'sza', 'phang', 'alt', 'wl', 'spe', 'mask']

pixels = []
for piz in cu:
    cose = []
    for name in names:
        if name == 'SPET' or name == 'WL' or name == 'BBL':
            cose.append(piz[name])
        else:
            cose.append(piz[name][0])
    pix = sbm.Pixel(keys=namesm,things=cose)
    pixels.append(pix)

pixels = np.array(pixels)
set_ch4 = sbm.PixelSet(pixels,descr='VIMS: Titan HCN-CH4-C2H2 climatology')

# Definisco le 4 regioni di interesse dello spettro
wl1 = [2950,3150]
wl2 = [3200,3300]
wl3 = [3300,3350]
wl4 = [3350,3500]

set_ch4.integr_tot(key='int_HCN',range=wl1)
set_ch4.integr_tot(key='int_R',range=wl2)
set_ch4.integr_tot(key='int_Q',range=wl3)
set_ch4.integr_tot(key='int_P',range=wl4)

set_ch4 = set_ch4[set_ch4.phang < 130]

setattr(set_ch4,'is_sun',np.array([np.mean(sp[:2])/np.mean(sp[-3:]) for sp in set_ch4.spe]))

alts = np.linspace(400,1100,15)
cond_lo = (set_ch4.alt < 500) & (set_ch4.alt > 400)
set_lo = set_ch4[cond_lo]
#for alt in alts