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

print('ciao')

cart = '/home/fede/Scrivania/Jiram/ANALISI/CH4_corr/FIT_CH4_TEMP/'
Ts = ['200','350','500','650','800']

sims_S = []
for t in Ts:
    wl_S, obs_S, sim = jirfu.read_sim_jir(cart+'OUT_ch4_SUD_T'+t+'/'+'sim_fin.dat')
    sims_S.append(sim)

sims_N = []
for t in Ts:
    wl_N, obs_N, sim = jirfu.read_sim_jir(cart+'OUT_ch4_NORD_T'+t+'/'+'sim_fin.dat')
    sims_N.append(sim)

sims_S_nlte = []
for t in Ts:
    wl_S, obs_S, sim = jirfu.read_sim_jir(cart+'OUT_ch4_SUDNLTE_T'+t+'/'+'sim_fin.dat')
    sims_S_nlte.append(sim)

names = ['Sim at '+t+' K' for t in Ts]
names_nlte = ['Sim at '+t+' K, r = 12' for t in Ts]
xlim = [3150,3550]
sbm.plot_spect_sim(cart+'Spect_SUD.pdf', wl_S, obs_S, sims_S, names, xscale = xlim, err=2.5e-7, title = r'CH$_4$ at South Pole')
sbm.plot_spect_sim(cart+'Spect_NORD.pdf', wl_N, obs_N, sims_N, names, xscale = xlim, err=2.5e-7, title = r'CH$_4$ at North Pole')
sbm.plot_spect_sim(cart+'Spect_SUD_nlte.pdf', wl_S, obs_S, sims_S_nlte, names, xscale = xlim, err=2.5e-7, title = r'CH$_4$ at South Pole, nlte factor r = 12')

xlim2 = [3280,3360]
sbm.plot_spect_sim(cart+'Spect_SUD_zoom.pdf', wl_S, obs_S, sims_S, names, xscale = xlim2, err=2.5e-7, title = r'CH$_4$ at South Pole')
sbm.plot_spect_sim(cart+'Spect_NORD_zoom.pdf', wl_N, obs_N, sims_N, names, xscale = xlim2, err=2.5e-7, title = r'CH$_4$ at North Pole')
sbm.plot_spect_sim(cart+'Spect_SUD_nlte_zoom.pdf', wl_S, obs_S, sims_S_nlte, names, xscale = xlim2, err=2.5e-7, title = r'CH$_4$ at South Pole, nlte factor r = 12')

sbm.plot_spect_sim(cart+'Spect_SUD_less.pdf', wl_S, obs_S, sims_S[::2], names[::2], xscale = xlim, err=2.5e-7, title = r'CH$_4$ at South Pole')
sbm.plot_spect_sim(cart+'Spect_NORD_less.pdf', wl_N, obs_N, sims_N[::2], names[::2], xscale = xlim, err=2.5e-7, title = r'CH$_4$ at North Pole')
sbm.plot_spect_sim(cart+'Spect_SUD_nlte_less.pdf', wl_S, obs_S, sims_S_nlte[::2], names[::2], xscale = xlim, err=2.5e-7, title = r'CH$_4$ at South Pole, nlte factor r = 12')


cart = '/home/fede/Scrivania/Jiram/ANALISI/SOLAR_OFFSET/'

#wl, obs_HI, obs_LO = pickle.load(open(cart + 'Spet_HILO.pic','r'))
wl, obs_HI, obs_LO = pickle.load(open(cart + 'Spet_HILO2.pic','r'))

esp = -4
fig = pl.figure(figsize=(8, 6), dpi=150)
pl.ylabel(r'Radiance ($\times 10^{{{}}}$ $W {{m}}^{{-2}} {{nm}}^{{-1}} {{sr}}^{{-1}}$)'.format(esp))
pl.xlabel('Wavelength (nm)')
pl.grid()
pl.plot(wl,obs_HI/10**esp,label='Day')
pl.plot(wl,obs_LO/10**esp,label='Night')
pl.legend(loc=1,fancybox=True)
fig.savefig(cart+'Spettro_HILO_scatt.pdf', format='pdf', dpi=150)
pl.close()

esp = -3
#oks = (wl > 3000.0) & (wl < 4500.0)
fig = pl.figure(figsize=(8, 6), dpi=150)
pl.ylabel(r'Radiance ($\times 10^{{{}}}\ W\ {{m}}^{{-2}}\ {{\mu m}}^{{-1}}\ {{sr}}^{{-1}}$)'.format(esp))
pl.xlabel(r'Wavelength ($\mu$m)')
#pl.grid()
pl.xlim(3.,4.5)
pl.ylim(-0.5,3.5)
# pl.plot(wl[oks],obs_HI[oks]/10**esp,label='Day')
# pl.plot(wl[oks],obs_LO[oks]/10**esp,label='Night')
pl.plot(wl*1e-3,1e3*obs_HI/10**esp,label='Day',color='red')
pl.plot(wl*1e-3,1e3*obs_LO/10**esp,label='Night',color='blue')
pl.legend(loc=1,fancybox=True)
fig.savefig(cart+'Spettro_HILO_scatt_zoom.pdf', format='pdf', dpi=150)

esp = -3
#oks = (wl > 3000.0) & (wl < 4500.0)
fig = pl.figure(figsize=(8, 6), dpi=150)
pl.ylabel(r'Radiance ($\times 10^{{{}}}$ $W {{m}}^{{-2}} {{nm}}^{{-1}} {{sr}}^{{-1}}$)'.format(esp))
pl.xlabel(r'Wavelength ($\mu$m)')
#pl.grid()
pl.xlim(3.,4.5)
pl.ylim(-0.5,3.5)
# pl.plot(wl[oks],obs_HI[oks]/10**esp,label='Day')
# pl.plot(wl[oks],obs_LO[oks]/10**esp,label='Night')
pl.plot(wl*1e-3,1e3*obs_HI/10**esp,label='Day')
pl.plot(wl*1e-3,1e3*obs_LO/10**esp,label='Night')
pl.legend(loc=1,fancybox=True)
fig.savefig(cart+'Spettro_HILO_scatt_zoom_bello.pdf', format='pdf', dpi=150)
