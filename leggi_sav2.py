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

cart = '/home/fede/Scrivania/Jiram/DATA/TEST_marisa/'
cart2 = '/home/fede/Scrivania/Jiram/DATA/TEST_marisa/Spe_prova_S_bis/'

nomeout = 'pix_nadir_S_bis.pic'
thres = 0.01 # soglia per l'indice di H3+

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

with open(cart+'lista_S','r') as lista:
    for line in lista:
        nome = line.rstrip()
        nomefi = line.rstrip()+'.sav'
        print(nome)
        #nome = '160826_205434_JUNO_SPE_00.sav'

        cubo = io.readsav(cart+nomefi)
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

                #cond = (geo[i,j,17] > 60.0 or geo[i,j,17] < -60.0) and (geo[i,j,28] == 1)# or geo[i,j,32] < 3000.0)
                #cond = geo[i,j,17] > 60.0 and (geo[i,j,28] == 1)# or geo[i,j,32] < 3000.0)
                cond = geo[i,j,17] < -60.0 and (geo[i,j,28] == 1)# or geo[i,j,32] < 3000.0)
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

                    if(ind_h3p < thres): continue
                    if(fu == 0): continue

                    found = True
                    pixcu = np.append(pixcu, pix)
        pixcu = pixcu[1:]
        if(found and len(pixcu)>10): pixtot = np.append(pixtot,pixcu)
        print(len(pixtot))

pixtot = pixtot[1:]
pickle.dump(pixtot,open(cart+nomeout,'w'))
