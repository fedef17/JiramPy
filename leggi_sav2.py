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

########################  SETUP #######################

# Folder where the original .SAV files are located
cart_orig = '/home/fede/Scrivania/Jiram/DATA/TEST_marisa/'

# File that contains the names of the .SAV files to be processed, with their relative path inside cart_orig. If set to None, all .SAV files in cart_orig are processed.
lista = None

# Folder where the output .pic file is located
cart_out = '/home/fede/Scrivania/Jiram/DATA/TEST_marisa/Spe_prova_S_bis/'

# Name of the output file:
nomeout = 'pix_nadir_S_bis.pic'

# FILTERING THE INPUT PIXELS.
thres_h3p = 0.0 # threshold for the H3+ index, calculated by the jirfu.ind_h3p() function. Select thres_h3p = 0 to include all data.

flag_nadir = False # If selected, only nadir measurements are retained.

flag_limb = False # If selected only limb measurements are retained.

max_limb_alt = 5000.0 # Only measurements up to max_limb_alt tangent altitude are retained.

min_lat = -90.0 # Minimum latitude of measurements considered.
max_lat = 90.0 # Maximum latitude of measurements considered.

###########################################################################################

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

pixtot = pix_

if lista_file is None:
    lista_tot = os.listdir(cart_orig)
    lista = [fil for fil in lista_tot if '.sav' in fil or '.SAV' in fil]
else:
    with open(cart_orig+lista_file,'r') as listaf:
        lista = [line.rstrip() for line in listaf.readlines()]

for nomefi in lista:
    print('Processing '+nomefi)
    ind = nomefi.index('.sav')
    if ind < 0:
        ind = nomefi.index('.SAV')
    nome = nomefi[:ind]
    print(nome)

    cubo = io.readsav(cart_orig+nomefi)
    spe = cubo.spe
    geo = cubo.geo
    wls = cubo.wls
    geo_names = cubo.geo_names
    lbls = cubo.lbls

    n_lin = np.shape(spe)[0]
    n_sam = np.shape(spe)[1]

    pixcu = pix_
    found = False
    for i in range(n_lin):
        for j in range(n_sam):
            corners = corners_
            pix = pix_

            if flag_nadir:
                cond = (geo[i,j,28] == 1)
            elif flag_limb:
                cond = (geo[i,j,28] == 0)
            else:
                cond = True

            geo[i,j,32] = 1e-3*geo[i,j,32] # Converting alt in km from m

            cond = cond and (geo[i,j,32] < max_limb_alt) and (geo[i,j,17] > min_lat) and (geo[i,j,17] < max_lat)

            if(cond):
                el = 0
                pix[0][el] = nome
                el+=1
                pix[0][el] = i
                el+=1
                pix[0][el] = j
                el+=1
                for l in range(13):
                    pix[0][el] = lbls[i][l][0]
                    el+=1
                for l in range(38):
                    pix[0][el] = geo[i,j,l]
                    el+=1

                pix['wl'][0] = wls
                pix['spe'][0] = spe[i,j,:]

                fondo = jirfu.fondojir(wls,spe[i,j,:])
                ind_h3p = jirfu.ind_h3p(wls,spe[i,j,:],fondo)

                pix['fondo'] = fondo
                pix['ind_h3p'] = ind_h3p

                fu = jirfu.checkqual(wls,spe[i,j,:],fondo)

                if(ind_h3p < thres_h3p): continue
                if(fu == 0): continue

                found = True
                pixcu = np.append(pixcu, pix)
    pixcu = pixcu[1:]
    pixtot = np.append(pixtot,pixcu)
    print(len(pixtot))

pixtot = pixtot[1:]
pickle.dump(pixtot,open(cart_out+nomeout,'w'))
