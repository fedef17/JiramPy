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
#from mpl_toolkits.basemap import Basemap

# import spiceypy.spiceypy as spice
# spice.furnsh('/home/fede/Scrivania/Jiram/DATA/KERNELS_JIRAM/Kernels_jm0003/jm0003.mk')

#cart = '/home/fede/Scrivania/Jiram/DATA/TEST_marisa/Spe_prova_N_bis/'
#cartres = cart + 'Res_CH4/'
#cart = '/home/fede/Scrivania/Jiram/DATA/JM0003_all/aur_S_nadir_70/'

nomi_c = ['pc_lon_1','pc_lon_2','pc_lon_3','pc_lon_4','pc_lat_1','pc_lat_2','pc_lat_3','pc_lat_4',
          'pg_lon_1','pg_lon_2','pg_lon_3','pg_lon_4','pg_lat_1','pg_lat_2','pg_lat_3','pg_lat_4',
          'loc_rad_1','loc_rad_2','loc_rad_3','loc_rad_4']
tipi_c = 19*'f8,'+'f8'
corners_ = np.empty(1, dtype=tipi_c)
corners_.dtype.names = nomi_c

nomi = ['cubo','line','sample',
        'product_id','start_time','stop_time','sc_lat','sc_lon','sc_azi','sc_alt','jup_dist',
        'sol_lat','sol_lon','sol_azi','exposure','inst_mode',
        'pc_lon_1','pc_lon_2','pc_lon_3','pc_lon_4','pc_lat_1','pc_lat_2','pc_lat_3','pc_lat_4',
        'pg_lon_1','pg_lon_2','pg_lon_3','pg_lon_4','pg_lat_1','pg_lat_2','pg_lat_3','pg_lat_4',
        'pc_lon','pc_lat','pg_lon','pg_lat',
        'loc_rad_1','loc_rad_2','loc_rad_3','loc_rad_4','loc_rad','incid_angle',
        'emiss_angle','phase_angle','flag_surf','alt_surf','slant_dist',
        'solar_time','tg_alt','tg_lon','tg_lat','tg_pha','tg_incid_angle',
        'tg_emiss_angle','wl','spe','flags','ind_h3p','fondo']#,'h3p_col','h3p_temp','wl_shift','offset']

tipi2 = '|S50,'+2*'i4,'+3*'|S50,'+9*'f8,'+'|S50,'+ 28 * 'f8,' + 'i4,' + 9*'f8,'+'O,O,O,'+'f8,'+'f8'#5*'f8,'+'f8'
pix_ = np.empty(1, dtype=tipi2)
pix_.dtype.names = nomi


cart = '/home/fede/Scrivania/Jiram/DATA/JM0003_all/aur_S_nadir_70_CH4/'
cartres = cart
cart500 = ''

#aur_lon_0,aur_lat,aur_lon,aur_theta = jirfu.leggi_map_aur(polo)

######################################################################################à

pixs = pickle.load(open(cart+'all_pixs.pic','r'))

# si va a prendere le geometrie per 500 km
tot500 = pickle.load(open(cart500+'all_pixs.pic','r'))

pixs500 = pix_
num = 0
for cubo,i,j in zip(pixs.cubo,pixs.line,pixs.sample):
    pi = pix_
    cond = (tot500.cubo == cubo) & (tot500.line == i) & (tot500.sample == j)
    if len(tot500[cond]) > 0:
        for l in range(len(pi)):
            pi[0][l] = tot500[cond][l]
    else:
        for l in range(3):
            pi[0][l] = pixs[num][l]
        for l in range(3,len(pi)):
            pi[0][l] = np.nan
    pixs500 = np.append(pixs500,pi)
    num+=1

print(len(pixs500))
pickle.dump(pixs500,open(cart+'all_pixs_500.pic','w'))