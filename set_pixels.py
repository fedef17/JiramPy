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

pixs_S,col_S,col_c_S,err_c_S,temp_S,err_t_S,chi_S,off_S,shi_S,ch4_S,err_ch4_S,nonan_S = jirfu.leggi_console('S')
pixs_N,col_N,col_c_N,err_c_N,temp_N,err_t_N,chi_N,off_N,shi_N,ch4_N,err_ch4_N,nonan_N = jirfu.leggi_console('N')

cartout = '/home/fede/Scrivania/Jiram/ANALISI/CH4_corr/'
try:
    os.stat(cartout)
except:
    os.mkdir(cartout)

pixarr_S = np.array([])
for pix in pixs_S:
    pix1 = jirfu.JirPix(keys = pix.dtype.names, things = pix)
    pixarr_S = np.append(pixarr_S,pix1)

pixarr_N = np.array([])
for pix in pixs_N:
    pix1 = jirfu.JirPix(keys = pix.dtype.names, things = pix)
    pixarr_N = np.append(pixarr_N,pix1)

# set_new = jirfu.JirSet(pixarr, descr = 'Aurora S, nuovo run con ch4 corretto')
# set_new.read_res(cartres)
#
# set_old = jirfu.JirSet(pixarr, descr = 'Aurora S, vecchio run con il ch4 sbagliato')
# set_old.read_res(cartpix)

cartlam_S = '/home/fede/Scrivania/Jiram/DATA/JM0003_all/aur_S_new_lampo/'
cartlam_N = '/home/fede/Scrivania/Jiram/DATA/JM0003_all/aur_N_new_lampo/'
cartlam_N80 = '/home/fede/Scrivania/Jiram/DATA/JM0003_all/aur_N_new_lampo/'

set_S = jirfu.JirSet(pixarr_S, descr = 'Aurora S, nuovo run lampo')
set_S.set_res([col_S,err_c_S,temp_S,err_t_S,chi_S,off_S,ch4_S,err_ch4_S])

set_N = jirfu.JirSet(pixarr_N, descr = 'Aurora N, nuovo run lampo')
set_N.set_res([col_N,err_c_N,temp_N,err_t_N,chi_N,off_N,ch4_N,err_ch4_N])

# set_new_sims = pickle.load(open(cartres+'sims.pic','r'))
# set_old_sims = pickle.load(open(cartpix+'sims.pic','r'))
set_lampo_sims_S = pickle.load(open(cartlam_S+'sims.pic','r'))
set_lampo_sims_N = pickle.load(open(cartlam_N+'sims.pic','r'))
set_lampo_sims_N80 = pickle.load(open(cartlam_N80+'sims.pic','r'))

wlsim_S = np.array(set_lampo_sims_S[0])[nonan_S]
simobs_S = np.array(set_lampo_sims_S[1])[nonan_S]
sim_S = np.array(set_lampo_sims_S[2])[nonan_S]
wlsim_N = np.array(set_lampo_sims_N[0]+set_lampo_sims_N80[0])[nonan_N]
simobs_N = np.array(set_lampo_sims_N[1]+set_lampo_sims_N80[1])[nonan_N]
sim_N = np.array(set_lampo_sims_N[2]+set_lampo_sims_N80[2])[nonan_N]

print(len(sim_N),len(pixs_N))
print(len(sim_S),len(pixs_S))
print('ciao')

# for i,wls,data,sim in zip(range(len(set_S)),set_lampo_sims_S[0],set_lampo_sims_S[1],set_lampo_sims_S[2]):
#     setattr(set_S[i],'wlsim',wls)
#     setattr(set_S[i],'obs',data)
#     setattr(set_S[i],'sim',sim)
#     setattr(set_S[i],'resid',data-sim)
#
# for i,wls,data,sim in zip(range(len(set_N)),set_lampo_sims_N[0]+set_lampo_sims_N80[0],set_lampo_sims_N[1]+set_lampo_sims_N80[1],set_lampo_sims_N[2]+set_lampo_sims_N80[2]):
#     setattr(set_N[i],'wlsim',wls)
#     setattr(set_N[i],'obs',data)
#     setattr(set_N[i],'sim',sim)
#     setattr(set_N[i],'resid',data-sim)

print('ciaociao')

cond_S = (set_S.ch4_col > 2e13) & (set_S.h3p_col < 5e11)
cond_N = (set_N.ch4_col > 0.5e13) & (set_N.h3p_col < 1e12)
print(len(set_S),len(set_S[cond_S]))
print(len(set_N),len(set_N[cond_N]))
#cond = (set_lampo.ch4_col > 3e13) & (set_lampo.h3p_col < 0.5e12)

# wlsim_S = set_S.wlsim[1]
# wlsim_N = set_N.wlsim[1]
#
# sims_S = np.nanmean(set_S.sim[cond_S],axis=0)
# sims_N = np.nanmean(set_N.sim[cond_N],axis=0)
#
# resid_S = np.nanmean(set_S.resid[cond_S],axis=0)
# resid_N = np.nanmean(set_N.resid[cond_N],axis=0)
#
# _dataS = np.nanmean(set_S.spe[cond_S],axis=0)
# _dataN = np.nanmean(set_N.spe[cond_N],axis=0)

wlsim_S = wlsim_S[1]
wlsim_N = wlsim_N[1]

print(np.shape(sim_S),np.shape(sim_N))
sims_ch4_S = np.nanmean(sim_S[cond_S],axis=0)
sims_ch4_N = np.nanmean(sim_N[cond_N],axis=0)
print(sims_ch4_N)

# resid_S = np.nanmean(set_S.resid[cond_S],axis=0)
# resid_N = np.nanmean(set_N.resid[cond_N],axis=0)

_dataS = np.nanmean(set_S.spe[cond_S],axis=0)
_dataN = np.nanmean(set_N.spe[cond_N],axis=0)

jirfu.write_obs_JIR(set_S.wl[1],_dataS,np.ones(len(_dataS),dtype=int),cartout+'CH4_spet_S.dat',comment='Average Ch4 spectrum for South pole (set_lampo.ch4_col > 2e13) & (set_lampo.h3p_col < 5e11)')

jirfu.write_obs_JIR(set_N.wl[1],_dataN,np.ones(len(_dataN),dtype=int),cartout+'CH4_spet_N.dat',comment='Average Ch4 spectrum for North pole (set_lampo.ch4_col > 1e13) & (set_lampo.h3p_col < 5e11)')

pl.figure(1)
pl.plot(set_S.wl[0],_dataS,label='Obs_S',linewidth=2)
pl.plot(wlsim_S, sims_ch4_S, label = 'Sims_S',linewidth=2)
pl.legend()
#pl.plot(wlsim_S, resid_S, label = 'Resid_S',linewidth=1)
pl.figure(2)
pl.plot(set_N.wl[0],_dataN,label='Obs_N',linewidth=2)
pl.plot(wlsim_N, sims_ch4_N, label = 'Sims_N',linewidth=2)
#pl.plot(wlsim_N, resid_N, label = 'Resid_N',linewidth=1)
pl.legend()
pl.show()


cond_S_day = (set_S.offset > 1.5*np.nanmean(set_S.offset))
cond_N_day = (set_N.offset > 1.5*np.nanmean(set_N.offset))
cond_S_night = (set_S.offset < 0.5*np.nanmean(set_S.offset))
cond_N_night = (set_N.offset < 0.5*np.nanmean(set_N.offset))
print(len(set_S),len(set_S[cond_S_day]),len(set_N[cond_N_night]))
print(len(set_N),len(set_N[cond_N_day]),len(set_N[cond_N_night]))
#cond = (set_lampo.ch4_col > 3e13) & (set_lampo.h3p_col < 0.5e12)

# wlsim_S = set_S.wlsim[1]
# wlsim_N = set_N.wlsim[1]
#
# sims_S = np.nanmean(set_S.sim[cond_S],axis=0)
# sims_N = np.nanmean(set_N.sim[cond_N],axis=0)
#
# resid_S = np.nanmean(set_S.resid[cond_S],axis=0)
# resid_N = np.nanmean(set_N.resid[cond_N],axis=0)
#
# _dataS = np.nanmean(set_S.spe[cond_S],axis=0)
# _dataN = np.nanmean(set_N.spe[cond_N],axis=0)

wlsim_S = wlsim_S[1]
wlsim_N = wlsim_N[1]

print(np.shape(sim_S),np.shape(sim_N))
sims_S_day = np.nanmean(sim_S[cond_S_day],axis=0)
sims_N_day = np.nanmean(sim_N[cond_N_day],axis=0)
sims_S_night = np.nanmean(sim_S[cond_S_night],axis=0)
sims_N_night = np.nanmean(sim_N[cond_N_night],axis=0)
print(sims_ch4_N)

# resid_S = np.nanmean(set_S.resid[cond_S],axis=0)
# resid_N = np.nanmean(set_N.resid[cond_N],axis=0)

_dataS_day = np.nanmean(set_S.spe[cond_S_day],axis=0)
_dataN_day = np.nanmean(set_N.spe[cond_N_day],axis=0)
_dataS_night = np.nanmean(set_S.spe[cond_S_night],axis=0)
_dataN_night = np.nanmean(set_N.spe[cond_N_night],axis=0)

cartout2 = '/home/fede/Scrivania/Jiram/ANALISI/SOLAR_OFFSET/'

jirfu.write_obs_JIR(set_S.wl[1],_dataS_day,np.ones(len(_dataS_day),dtype=int),cartout2+'HO_spet_S.dat',comment='Average spectrum with large offset at South pole (set_lampo.offset > 1.5*mean)')
jirfu.write_obs_JIR(set_S.wl[1],_dataS_night,np.ones(len(_dataS_night),dtype=int),cartout2+'LO_spet_S.dat',comment='Average spectrum with small offset at South pole (set_lampo.offset < 0.5*mean)')
jirfu.write_obs_JIR(set_N.wl[1],_dataN_day,np.ones(len(_dataN_day),dtype=int),cartout2+'HO_spet_N.dat',comment='Average spectrum with large offset at Nouth pole (set_lampo.offset > 1.5*mean)')
jirfu.write_obs_JIR(set_N.wl[1],_dataN_night,np.ones(len(_dataN_night),dtype=int),cartout2+'LO_spet_N.dat',comment='Average spectrum with small offset at Nouth pole (set_lampo.offset < 0.5*mean)')

pl.figure(1)
pl.plot(set_S.wl[0],_dataS_day,label='Obs_S_HO',linewidth=2)
pl.plot(set_S.wl[0],_dataS_night,label='Obs_S_LO',linewidth=2)
pl.legend()
#pl.plot(wlsim_S, resid_S, label = 'Resid_S',linewidth=1)
pl.figure(2)
pl.plot(set_N.wl[0],_dataN_day,label='Obs_N_HO',linewidth=2)
pl.plot(set_N.wl[0],_dataN_night,label='Obs_N_LO',linewidth=2)
#pl.plot(wlsim_N, resid_N, label = 'Resid_N',linewidth=1)
pl.legend()
pl.show()


cartout = '/home/fede/Scrivania/Jiram/ANALISI/SOLAR_OFFSET/'

set_N.integr_tot(key='emiss',range=[4700,5000])
set_N.integr_tot(key='sscatt2',range=[2700,2750])

okHI = (set_N.sscatt2 > 0.0005) & (set_N.emiss < 0.05)
okLO = (set_N.sscatt2 < 0.0001) & (set_N.emiss < 0.05)

# _data_HI2 = np.nanmean(set_N.spe[okHI],axis=0)
# _data_LO2 = np.nanmean(set_N.spe[okLO],axis=0)
pickle.dump([set_N[0].wl,set_N.spe[okHI][10],set_N.spe[okLO][10]],file=open(cartout+'Spet_HILO2.pic','w'))

# off=1e-6
# pl.plot(wlsim,-np.ones(len(wlsim))*off,linestyle=':',color='black')
# pl.plot(wlsim,_new,label='Sim_new')
# pl.plot(wlsim,_old,label='Sim_old')
# pl.plot(wlsim,_lampo,label='Sim_lampo')
# pl.plot(wlsim,mea_new-off,label='Resid_new')
# pl.plot(wlsim,mea_old-off,label='Resid_old')
# pl.plot(wlsim,mea_lampo-off,label='Resid_lampo')
# pl.grid()
# pl.legend()
# pl.show()

# new_ = np.mean([sim for sim in set_new_sims[2]], axis=0)
# old_ = np.mean([sim for sim in set_old_sims[2]], axis=0)
# lampo_ = np.mean([sim for sim in set_lampo_sims[2]], axis=0)


# pickle.dump(set_old,open(cartpix+'set_pix.pic','w'))
# pickle.dump(set_new,open(cartres+'set_pix.pic','w'))
# pickle.dump(set_lampo,open(cartlam+'set_pix.pic','w'))

# print('ciao')
#
# pdf = PdfPages(cartout+'Differenze_ch4_corr.pdf')
#
# titolo = 'New - Old run (after ch4 correction)'
# lab = 'new_'
#
#
# fig = pl.figure(figsize=(8, 6), dpi=150)
# pl.hist((set_new.h3p_col-set_old.h3p_col),bins=40,range= [-4e12,4e12])
# pl.title(titolo)
# pl.xlabel('Diff in h3p columns [cm-2]')
# pl.grid()
# fig.savefig(cartout+lab+'diff_S_h3p_col.eps', format='eps', dpi=150)
# pdf.savefig(fig)
# pl.close()
#
# fig = pl.figure(figsize=(8, 6), dpi=150)
# pl.hist((set_new.ch4_col/set_old.ch4_col),bins=np.logspace(0,3,50))
# pl.title(titolo)
# pl.xlabel('Ratio new/old ch4 columns (new temp = 210 K, old = 250 K)')
# pl.xscale('log')
# pl.grid()
# fig.savefig(cartout+lab+'diff_S_ch4_col.eps', format='eps', dpi=150)
# pdf.savefig(fig)
# pl.close()
#
#
# fig = pl.figure(figsize=(8, 6), dpi=150)
# pl.hist((set_new.h3p_temp-set_old.h3p_temp),bins=40,range= [-300,300])
# pl.title(titolo)
# pl.xlabel('Diff in h3p temp [K]')
# pl.grid()
# fig.savefig(cartout+lab+'diff_S_h3p_temp.eps', format='eps', dpi=150)
# pdf.savefig(fig)
# pl.close()
#
#
# fig = pl.figure(figsize=(8, 6), dpi=150)
# pl.grid(color='white')
# pl.hist2d((set_new.h3p_temp-set_old.h3p_temp),set_new.h3p_col,bins=[40,40],range= [[-300,300],[0,1e13]])
# pl.colorbar()
# pl.title(titolo)
# pl.xlabel('Diff in h3p temp [K]')
# pl.ylabel('h3p col [cm-2]')
# pl.grid()
# fig.savefig(cartout+lab+'diff_2d_S_h3p_temp.eps', format='eps', dpi=150)
# pdf.savefig(fig)
# pl.close()
#
#
# set_new_2 = set_new
# set_new = set_lampo
# titolo = 'Lampo - Old run (after ch4 correction)'
# lab = 'lampo_'
#
# fig = pl.figure(figsize=(8, 6), dpi=150)
# pl.hist((set_new.h3p_col-set_old.h3p_col),bins=40,range= [-4e12,4e12])
# pl.title(titolo)
# pl.xlabel('Diff in h3p columns [cm-2]')
# pl.grid()
# fig.savefig(cartout+lab+'diff_S_h3p_col.eps', format='eps', dpi=150)
# pdf.savefig(fig)
# pl.close()
#
# fig = pl.figure(figsize=(8, 6), dpi=150)
# pl.hist((set_new.ch4_col/set_old.ch4_col),bins=np.logspace(-5,-1,50))
# pl.title(titolo)
# pl.xlabel('Ratio new/old ch4 columns (new temp = 500 K, old = 250 K)')
# pl.xscale('log')
# pl.grid()
# fig.savefig(cartout+lab+'diff_S_ch4_col.eps', format='eps', dpi=150)
# pdf.savefig(fig)
# pl.close()
#
# fig = pl.figure(figsize=(8, 6), dpi=150)
# pl.hist((set_new.h3p_temp-set_old.h3p_temp),bins=40,range= [-300,300])
# pl.title(titolo)
# pl.xlabel('Diff in h3p temp [K]')
# pl.grid()
# fig.savefig(cartout+lab+'diff_S_h3p_temp.eps', format='eps', dpi=150)
# pdf.savefig(fig)
# pl.close()
#
#
# fig = pl.figure(figsize=(8, 6), dpi=150)
# pl.grid(color='white')
# pl.hist2d((set_new.h3p_temp-set_old.h3p_temp),set_new.h3p_col,bins=[40,40],range= [[-300,300],[0,1e13]])
# pl.colorbar()
# pl.title(titolo)
# pl.xlabel('Diff in h3p temp [K]')
# pl.ylabel('h3p col [cm-2]')
# pl.grid()
# fig.savefig(cartout+lab+'diff_2d_S_h3p_temp.eps', format='eps', dpi=150)
# pdf.savefig(fig)
# pl.close()
#
#
# set_old = set_new_2
# set_new = set_lampo
# titolo = 'Lampo - New run (after ch4 correction)'
# lab = 'lampo_'
#
# fig = pl.figure(figsize=(8, 6), dpi=150)
# pl.hist((set_new.h3p_col-set_old.h3p_col),bins=40,range= [-4e12,4e12])
# pl.title(titolo)
# pl.xlabel('Diff in h3p columns [cm-2]')
# pl.grid()
# fig.savefig(cartout+lab+'diff_S_h3p_col.eps', format='eps', dpi=150)
# pdf.savefig(fig)
# pl.close()
#
# fig = pl.figure(figsize=(8, 6), dpi=150)
# pl.hist((set_new.ch4_col/set_old.ch4_col),bins=np.logspace(-7,-3,50))
# pl.title(titolo)
# pl.xlabel('Ratio new/old ch4 columns (new temp = 500 K, old = 210 K)')
# pl.xscale('log')
# pl.grid()
# fig.savefig(cartout+lab+'diff_S_ch4_col.eps', format='eps', dpi=150)
# pdf.savefig(fig)
# pl.close()
#
# fig = pl.figure(figsize=(8, 6), dpi=150)
# pl.hist((set_new.h3p_temp-set_old.h3p_temp),bins=40,range= [-300,300])
# pl.title(titolo)
# pl.xlabel('Diff in h3p temp [K]')
# pl.grid()
# fig.savefig(cartout+lab+'diff_S_h3p_temp.eps', format='eps', dpi=150)
# pdf.savefig(fig)
# pl.close()
#
#
# fig = pl.figure(figsize=(8, 6), dpi=150)
# pl.grid(color='white')
# pl.hist2d((set_new.h3p_temp-set_old.h3p_temp),set_new.h3p_col,bins=[40,40],range= [[-300,300],[0,1e13]])
# pl.colorbar()
# pl.title(titolo)
# pl.xlabel('Diff in h3p temp [K]')
# pl.ylabel('h3p col [cm-2]')
# pl.grid()
# fig.savefig(cartout+lab+'diff_2d_S_h3p_temp.eps', format='eps', dpi=150)
# pdf.savefig(fig)
# pl.close()
#
# pdf.close()
