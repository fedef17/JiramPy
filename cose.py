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
#
# cart2 = '/home/fede/Scrivania/Dotto/AbstrArt/CH4_HCN_climatology/DATA/'
# cubo_ch4 = io.readsav(cart+'PIXs_HCN-CH4-C2H2_season.sav')

pixs,pixs500,col,col_c,err_c,temp,err_t,chi,off,shi,ch4 = jirfu.leggi_console('N')

edges=[]
for pi in zip(pixs500.pc_lon_1, pixs500.pc_lon_2, pixs500.pc_lon_3, pixs500.pc_lon_4, pixs500.pc_lat_1, pixs500.pc_lat_2, pixs500.pc_lat_3, pixs500.pc_lat_4):
    vec = [[pi[i],pi[i+4]] for i in range(4)]
    edges.append(vec)

cart = '/home/fede/Scrivania/Jiram/ANALISI/'
jirfu.stereoplot(pixs500.pc_lon,pixs500.pc_lat,col_c,polo='N',minu=0.2e12,maxu=3e12,title='CIAOAOAOA',live=True,aur_model='both',edges=edges,nomefi=cart+'Aur_N_pixels.pdf')

sys.exit()

def read_sim_jir(filename):
    infile = open(filename, 'r')
    data = [line.split() for line in infile]
    data_arr = np.array(data)
    wl = np.array([float(r) for r in data_arr[:, 0]])
    data = np.array([float(r) for r in data_arr[:, 1]])
    sim = np.array([float(r) for r in data_arr[:, 2]])
    infile.close()
    return wl, data, sim

wls = []
datas = []
sims = []

cart = '/home/federico/JIRAM/JIRAM_MAP_S/Results_aur_S_nadir_70_ch4/'
ntot = 32338
for i in range(ntot):
    fi = 'OUT_{0:05d}/sim_fin.dat'.format(i)
    wl,data,sim = read_sim_jir(cart+fi)
    wls.append(wl)
    datas.append(data)
    sims.append(sim)

pickle.dump([wls,datas,sims],open(cart+'sims.pic','w'))

# pixs,pixs500,col,col_c,err_c,temp,err_t,chi,off,shi,ch4 = jirfu.leggi_console('N')
#
# cart = '/home/fede/Scrivania/Jiram/ANALISI/'
#
# integr = np.zeros(len(pixs))
#
# for ww,sp,of,i in zip(pixs500.wl,pixs500.spe,off,range(len(pixs500))):
#     if(ww is None): continue
#     print(i)
#     mask = jirfu.findspi(ww,sp)
#     integr[i] = jirfu.integr_h3p(ww,sp*mask,of,w1=3200,w2=3700)
#
# fi = open(cart+'integr_N.dat','w')
# for i in range(len(pixs)):
#     if np.isnan(integr[i]) or integr[i]==0:
#         continue
#     fi.write(('{:10.4e}'+5*'{:8.2f}'+'\n').format(integr[i],pixs500.emiss_angle[i],pixs500.incid_angle[i],pixs500.pc_lat[i],pixs500.pc_lon[i],pixs500.solar_time[i]))
#
# fi.close()