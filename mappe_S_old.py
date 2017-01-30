#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
import matplotlib.pyplot as pl
import math as mt
import spect_base_module as sbm
import jirfu
from subprocess import call
import pickle
import scipy.io as io
from scipy.interpolate import PchipInterpolator as spline
from matplotlib.backends.backend_pdf import PdfPages


####################################################################################
#############                     SETTINGS                            ##############
####################################################################################

# Da dove leggo i risultati:
# cart = '/home/federico/JIRAM/JIRAM_MAP_S/Results_new_S_lampo/'

cart = '/home/fede/Scrivania/Jiram/DATA/JM0003_all/aur_S_nadir_70_CH4/'
cartres = cart

# Da dove leggo il file .pic originale e il .pic per i 500 km:
# picfile = '/home/federico/JIRAM/MAPPE/all_pixs_S_0km.pic'
# picfile500 = '/home/federico/JIRAM/MAPPE/all_pixs_S_500km.pic'

picfile500 = '/home/fede/Scrivania/Jiram/DATA/JM0003_all/aur_S_nadir_70_CH4/all_pixs_500.pic'
picfile = '/home/fede/Scrivania/Jiram/DATA/JM0003_all/aur_S_nadir_70_CH4/all_pixs.pic'
#picfile80 = '/home/fede/Scrivania/Jiram/DATA/JM0003_all/aur_N_nadir_80_CH4/all_pixs.pic'

# Dove metto le mappe:
# cartout = '/home/federico/JIRAM/MAPPE/MAPPE_S_lampo/'
cartout = '/home/fede/Scrivania/Jiram/ANALISI/MAPPE/MAPPE_S_lampo/'
if not os.path.exists(cartout): os.makedirs(cartout)

polo = 'S'
lch4 = True  # ho fittato il ch4?
l500 = True  # così fa tutto con le geometrie a 500 km
lshi = False # ho fittato la wl shift?

# settings base

maxchi = 20 # massimo chi quadro considerato
mcol = 4e12  # massimo della scala per la colonna di H3+
mincol = 0.
npix = 50  # risoluzione della mappa (npix x npix pixels)
aurm = 'both'  # così stampa entrambi gli anelli. Se 'stat' solo quello statistico, se 'VIP4' solo quello vip4.
tmax = 1100
cstep = 0.25e12
tstep = 25

cbarlabcol = r'$H_3^+$ column ($\times 10^{{{}}}$ cm$^{{-2}}$)'
cbarlabtemp = r'$H_3^+$ temperature (K)'
cbarlabch4 = r'CH$_4$ column ($\times 10^{{{}}}$ cm$^{{-2}}$)'

cbarlabint = r'Integrated intensity ($\times 10^{{{}}}$ $W {{m}}^{{-2}} {{sr}}^{{-1}}$)'
titleintegr = 'Integrated intensity between 3.35 and 3.75 nm'

titlecol = r'$H_3^+$ column : FILT {}'
titletemp = r'$H_3^+$ temperature : FILT {}'
titlech4 = r'CH$_4$ column : FILT {}'

####################################################################################
####################                                   #############################
#############                   LEGGE                         #######################
###################                                    #############################
####################################################################################

pixs = pickle.load(open(picfile,'r'))

if l500:
    # si va a prendere le geometrie per 500 km
    p500 = pickle.load(open(picfile500,'r'))
    pixs500 = p500.view(np.recarray)
    print(len(pixs),len(pixs500))
    #print(pixs[1139])
    #print(pixs500[1139])

ii,col1, err_col = jirfu.read_res_jir_3(cartres+'CD-H3p.dat')
iit,temp1, err_temp = jirfu.read_res_jir_3(cartres+'VT-H3p.dat')
iic,chi1 = jirfu.read_res_jir_4(cartres+'chisq.dat')
iio,off1 = jirfu.read_res_jir_4(cartres+'offset_ok.dat')
if lshi: iis,shi1 = jirfu.read_res_jir_4(cartres+'shift_ok.dat')
if lch4: ii4,ch41, err_ch4 = jirfu.read_res_jir_3(cartres+'CD-CH4.dat')

col = np.zeros(len(pixs))
err_c = np.zeros(len(pixs))
temp = np.zeros(len(pixs))
err_t = np.zeros(len(pixs))
chi = np.zeros(len(pixs))
off = np.zeros(len(pixs))
shi = np.zeros(len(pixs))
ch4 = np.zeros(len(pixs))

col[ii] = col1
err_c[ii] = err_col
temp[iit] = temp1
err_t[iit] = err_temp
chi[iic] = chi1
chi[chi == 0.0] = np.nan
off[iio] = off1
if lshi: shi[iis] = shi1
if lch4: ch4[ii4] = ch41

#print(type(pixs.cubo))
cubs = np.unique(pixs500.cubo)
i=0
for cu in cubs:
    if not np.all(np.isnan(pixs500[pixs500.cubo == cu].pc_lon)):
        em = np.mean(pixs500[pixs500.cubo == cu].emiss_angle)
        so = np.mean(pixs500[pixs500.cubo == cu].solar_time)
        n = len(pixs500[pixs500.cubo == cu])
        print('{}  {}  {:5.1f}  {:5.1f} {}'.format(i,cu,em,so,n))

######## fine letture ################################################################################

# Faccio le mappe: #######################################

cond2 = ((chi > maxchi) | (col <= 0.0) | (off < 0.0))

col[cond2] = float('NaN')
print(len(col[cond2]))
temp[cond2] = float('NaN')
chi[cond2] = float('NaN')
off[cond2] = float('NaN')
err_t[cond2] = float('NaN')
err_c[cond2] = float('NaN')

if lch4: ch4[cond2] = float('NaN')
if lch4: ch4[(ch4<0.0)] = 0.0

# fi = open(cart+'lista_ok_ch4.dat','w')
# i=0
# coch4 = ch4 > 3e15
# iii = np.arange(len(pixs))
# for pi,i in zip(pixs[~cond2 & coch4],iii[~cond2 & coch4]):
#     fi.write('{0:05d}\n'.format(i))
# fi.close()
# sys.exit()

if(l500): pixs = pixs500

indx = np.arange(len(pixs))

nonan = (~np.isnan(col)) & (~np.isnan(pixs.pc_lon))
col = col[nonan]
temp = temp[nonan]
chi = chi[nonan]
off = off[nonan]
err_t = err_t[nonan]
err_c = err_c[nonan]
if lch4: ch4 = ch4[nonan]
pixs = pixs[nonan]
indx = indx[nonan]

col_c = col*np.cos(np.pi*pixs.emiss_angle/180.0)
err_cc = err_c*np.cos(np.pi*pixs.emiss_angle/180.0)

lon = pixs.pc_lon
lat = pixs.pc_lat
lab = '_ok1'

edges=[]
for pi in zip(pixs500.pc_lon_1, pixs500.pc_lon_2, pixs500.pc_lon_3, pixs500.pc_lon_4, pixs500.pc_lat_1, pixs500.pc_lat_2, pixs500.pc_lat_3, pixs500.pc_lat_4):
    vec = [[pi[i],pi[i+4]] for i in range(4)]
    edges.append(vec)

##################################################################

# Condizione:

pdf_col = PdfPages(cartout+'MAPPE_h3p_col.pdf')
pdf_temp = PdfPages(cartout+'MAPPE_h3p_temp.pdf')
pdf_scatt_1 = PdfPages(cartout+'MAPPE_scatt_col_temp.pdf')
pdf_scatt_2 = PdfPages(cartout+'MAPPE_scatt_errc_errt.pdf')
pdf_scatt_3 = PdfPages(cartout+'MAPPE_scatt_errcRel_errt.pdf')

#for filt in [0,1,2,3,4]:
for filt in [0,1]:
    print('TOT: {}'.format(len(pixs)))

    if filt == 0:
        coto = (~np.isnan(col))
        lab = '_nofil'
        print('OK'+lab+': {}'.format(len(pixs[coto])))
        cartout2 = cartout + lab + '/'
        if not os.path.exists(cartout2): os.makedirs(cartout2)
    elif filt == 1:
        coto = (~np.isnan(col)) & (err_t < 100)
        lab = '_fil1'
        print('OK'+lab+': {}'.format(len(pixs[coto])))
        cartout2 = cartout + lab + '/'
        if not os.path.exists(cartout2): os.makedirs(cartout2)
    elif filt == 2:
        coto = (~np.isnan(col)) & (err_c/col < 0.3) & (err_t < 100)
        lab = '_fil2'
        print('OK'+lab+': {}'.format(len(pixs[coto])))
        cartout2 = cartout + lab + '/'
        if not os.path.exists(cartout2): os.makedirs(cartout2)
    elif filt == 3:
        coto = (~np.isnan(col)) & (err_c < 0.5e12) & (err_t < 100)
        lab = '_fil3'
        print('OK'+lab+': {}'.format(len(pixs[coto])))
        cartout2 = cartout + lab + '/'
        if not os.path.exists(cartout2): os.makedirs(cartout2)
    elif filt == 4:
        coto = (~np.isnan(col)) & (col > 1.5e12) & (err_t < 100)
        lab = '_fil4'
        print('OK'+lab+': {}'.format(len(pixs[coto])))
        cartout2 = cartout + lab + '/'
        if not os.path.exists(cartout2): os.makedirs(cartout2)
    ##################################################################

    jirfu.stereomap2(lon[coto],lat[coto],col_c[coto],nomefi=cartout2+'MAP_h3p_col'+lab+'.eps',title=titlecol.format(filt),polo=polo,cbarlabel=cbarlabcol,cbarmult=12,cbarform='%.1f',step=cstep,min=mincol,max=mcol, show=False, xres=npix, image=False, aur_model=aurm,pdf=pdf_col)
    jirfu.stereomap2(lon[coto],lat[coto],err_cc[coto]/col_c[coto],nomefi=cartout2+'ERR_h3p_col'+lab+'.eps',title='Relative error H3+ column [mol/cm2] : FILT {}'.format(filt),polo=polo,min=0.,max=1., show=False,xres=npix)
    jirfu.stereomap2(lon[coto],lat[coto],col_c[coto],nomefi=cartout2+'Npoints_ok'+lab+'.eps',title='Number of measurements : FILT {}'.format(filt),polo=polo,show=False,xres=50,npoints=True,aur_model=aurm)#,addpoints=True)

    co = pixs.emiss_angle < 60
    jirfu.stereomap2(lon[co & coto],lat[co & coto],col_c[co & coto],nomefi=cartout2+'MAP_h3p_col'+lab+'_emissm60.eps',title='H3+ column [mol/cm2] : emiss < 60',polo=polo,cbarlabel=cbarlabcol,cbarmult=12,cbarform='%.1f',step=cstep,min=mincol,max=mcol, show=False, xres=npix,aur_model=aurm,addpoints=True)
    #jirfu.stereoplot(lon[co & coto],lat[co & coto],col_c[co & coto],nomefi=cartout2+'PLOT_h3p_col'+lab+'_emissm60.eps',title='H3+ column [mol/cm2] : emiss < 60 ',polo=polo,min=mincol,max=mcol,show=False,aur_model=aurm)
    co = pixs.emiss_angle > 60
    jirfu.stereomap2(lon[co & coto],lat[co & coto],col_c[co & coto],nomefi=cartout2+'MAP_h3p_col'+lab+'_emissp60.eps',title='H3+ column [mol/cm2] : emiss > 60',polo=polo,cbarlabel=cbarlabcol,cbarmult=12,cbarform='%.1f',step=cstep,min=mincol,max=mcol, show=False, xres=npix,aur_model=aurm,addpoints=True)
    #jirfu.stereoplot(lon[co & coto],lat[co & coto],col_c[co & coto],nomefi=cartout2+'PLOT_h3p_col'+lab+'_emissp60.eps',title='H3+ column [mol/cm2] : emiss > 60',polo=polo,min=mincol,max=mcol,show=False,aur_model=aurm)

    cond4 = (pixs.alt_surf < 4e8)
    nome = cartout2+'MAP_h3p_col_close'+lab+'.eps'
    jirfu.stereomap2(lon[cond4 & coto],lat[cond4 & coto],col_c[cond4 & coto],nomefi=nome,title='H3+ column [mol/cm2] -> alt_surf < 4e5 km',polo=polo,cbarlabel=cbarlabcol,cbarmult=12,cbarform='%.1f',step=cstep,min=mincol,max=mcol,show=False,xres=npix,aur_model=aurm)
    cond4 = (pixs.alt_surf > 4e8)
    nome = cartout2+'MAP_h3p_col_far'+lab+'.eps'
    jirfu.stereomap2(lon[cond4 & coto],lat[cond4 & coto],col_c[cond4 & coto],nomefi=nome,title='H3+ column [mol/cm2] -> alt_surf > 4e5 km',polo=polo,cbarlabel=cbarlabcol,cbarmult=12,cbarform='%.1f',step=cstep,min=mincol,max=mcol,show=False,xres=npix,aur_model=aurm)
    cond4 = (pixs.solar_time < 12) & (pixs.solar_time > 6)
    nome = cartout2+'MAP_h3p_col_am'+lab+'.eps'
    jirfu.stereomap2(lon[cond4 & coto],lat[cond4 & coto],col_c[cond4 & coto],nomefi=nome,title='H3+ column [mol/cm2] -> 6h < solar time < 12h',polo=polo,cbarlabel=cbarlabcol,cbarmult=12,cbarform='%.1f',step=cstep,min=mincol,max=mcol,show=False,xres=npix,aur_model=aurm)
    cond4 = (pixs.solar_time > 12) & (pixs.solar_time < 18)
    nome = cartout2+'MAP_h3p_col_pm'+lab+'.eps'
    jirfu.stereomap2(lon[cond4 & coto],lat[cond4 & coto],col_c[cond4 & coto],nomefi=nome,title='H3+ column [mol/cm2] -> 12h < solar time < 18h',polo=polo,cbarlabel=cbarlabcol,cbarmult=12,cbarform='%.1f',step=cstep,min=mincol,max=mcol,show=False,xres=npix,aur_model=aurm)
    cond4 = (pixs.solar_time < 24) & (pixs.solar_time > 18)
    nome = cartout2+'MAP_h3p_col_ni1'+lab+'.eps'
    jirfu.stereomap2(lon[cond4 & coto],lat[cond4 & coto],col_c[cond4 & coto],nomefi=nome,title='H3+ column [mol/cm2] -> 18h < solar time < 24h',polo=polo,cbarlabel=cbarlabcol,cbarmult=12,cbarform='%.1f',step=cstep,min=mincol,max=mcol,show=False,xres=npix,aur_model=aurm)
    cond4 = (pixs.solar_time < 6) & (pixs.solar_time > 0)
    nome = cartout2+'MAP_h3p_col_ni2'+lab+'.eps'
    jirfu.stereomap2(lon[cond4 & coto],lat[cond4 & coto],col_c[cond4 & coto],nomefi=nome,title='H3+ column [mol/cm2] -> 0h < solar time < 6h',polo=polo,cbarlabel=cbarlabcol,cbarmult=12,cbarform='%.1f',step=cstep,min=mincol,max=mcol,show=False,xres=npix,aur_model=aurm)

    # giu = err_t < 100
    lab2 = ''
    jirfu.stereomap2(lon[coto],lat[coto],temp[coto],nomefi=cartout2+'MAP_h3p_temp'+lab2+lab+'.eps',title=titletemp.format(filt), polo=polo,cbarlabel=cbarlabtemp,cbarform='%.0f',step=tstep,min=700.,max=tmax, show=False,xres=npix,aur_model=aurm,pdf=pdf_temp)
    jirfu.stereomap2(lon[coto],lat[coto],err_t[coto],nomefi=cartout2+'ERR_h3p_temp'+lab2+lab+'.eps',title='Error H3+ temperature [K] : FILT {}'.format(filt), polo=polo,min=0.,max=200., show=False,xres=npix,aur_model=aurm)
    # giu = err_t < 80
    # lab2 = '_err80K'
    # jirfu.stereomap2(lon[giu],lat[giu],temp[giu],nomefi=cartout2+'MAP_h3p_temp'+lab2+lab+'.eps',title='H3+ temperature [K]', polo=polo,min=700.,max=tmax, show=False,xres=npix,aur_model=aurm)
    # jirfu.stereomap2(lon[giu],lat[giu],err_t[giu],nomefi=cartout2+'ERR_h3p_temp'+lab2+lab+'.eps',title='Error H3+ temperature [K]', polo=polo, min=0.,max=200., show=False,xres=npix,aur_model=aurm)
    # giu = err_t < 50
    # lab2 = '_err50K'
    # jirfu.stereomap2(lon[giu],lat[giu],temp[giu],nomefi=cartout2+'MAP_h3p_temp'+lab2+lab+'.eps',title='H3+ temperature [K]', polo=polo,min=700.,max=tmax, show=False,xres=npix,aur_model=aurm)
    # jirfu.stereomap2(lon[giu],lat[giu],err_t[giu],nomefi=cartout2+'ERR_h3p_temp'+lab2+lab+'.eps',title='Error H3+ temperature [K]', polo=polo, min=0.,max=200., show=False,xres=npix,aur_model=aurm)
    # condt = (col > 0e12)
    # lab2 = '_nofilt'
    # jirfu.stereomap2(lon[condt],lat[condt],temp[condt],nomefi=cartout2+'MAP_h3p_temp'+lab2+lab+'.eps',title='H3+ temperature [K]', polo=polo,min=700.,max=tmax, show=False,xres=npix,aur_model=aurm)
    # jirfu.stereomap2(lon[condt],lat[condt],err_t[condt],nomefi=cartout2+'ERR_h3p_temp'+lab2+lab+'.eps',title='Error H3+ temperature [K]', polo=polo, min=0.,max=200., show=False,xres=npix,aur_model=aurm)

    jirfu.stereomap2(lon[coto],lat[coto],ch4[coto],nomefi=cartout2+'MAP_ch4_col'+lab+'.eps',title=titlech4,cbarlabel=cbarlabch4,cbarform='%.1f',cbarmult=13,min=0,max=1.5e13,step=1.5e12,polo=polo, show=False,xres=npix,aur_model=aurm)


    ## l'emissione integrata:
    condpi = ~np.isnan(pixs.pc_lon)
    zu = pixs[condpi & coto]
    integr = np.zeros(len(zu))
    for ww,spp,of,i in zip(zu.wl,zu.spe,off[condpi & coto],range(len(zu))):
	nano = np.isnan(spp)
	spp[nano] = 0.0
        mask = jirfu.findspi(ww,spp)
        integr[i] = jirfu.integr_h3p(ww,spp*mask,of)

    integr=integr*np.cos(np.pi*zu.emiss_angle/180.0)
    jirfu.stereomap2(lon,lat,integr,nomefi=cartout2+'MAP_integr_h3p_conslit'+lab+'.eps',title=titleintegr,step=1.5e-5,
                    polo=polo,cbarlabel=cbarlabint,cbarmult=-4,cbarform='%.1f',min=0.0,max=1.5e-4,addpoints=True,show=False,xres=npix,aur_model=aurm)
    jirfu.stereomap2(lon,lat,integr,nomefi=cartout2+'MAP_integr_h3p'+lab+'.eps',title=titleintegr,step=1.5e-5,
                    polo=polo,cbarlabel=cbarlabint,cbarmult=-4,cbarform='%.1f',min=0.0,max=1.5e-4,addpoints=False,show=False,xres=npix,aur_model=aurm)


    # i plot punto per punto
    jirfu.stereoplot(lon[coto],lat[coto],col_c[coto],nomefi=cartout2+'h3p_col'+lab+'.eps',title='H3+ column [mol/cm2]',polo=polo,cbarlabel=cbarlabcol,cbarmult=12,cbarform='%.1f',step=cstep,min=mincol,max=mcol,aur_model=aurm)
    jirfu.stereoplot(lon[coto],lat[coto],col[coto],nomefi=cartout2+'h3p_col_nocorr'+lab+'.eps',title='H3+ column [mol/cm2] (non corrected)',polo=polo,cbarlabel=cbarlabcol,cbarmult=12,cbarform='%.1f',step=cstep,min=mincol,max=8e12,aur_model=aurm)
    jirfu.stereoplot(lon[coto],lat[coto],temp[coto],nomefi=cartout2+'h3p_temp'+lab+'.eps',title='H3+ temperature [K]',polo=polo,cbarlabel=cbarlabtemp,cbarform='%.0f',step=tstep,min=700.,max=tmax,aur_model=aurm)
    if lch4: jirfu.stereoplot(lon[coto],lat[coto],ch4[coto],nomefi=cartout2+'ch4_col'+lab+'.eps',title='Effective CH4 column [mol/cm2]',polo=polo,aur_model=aurm)
    jirfu.stereoplot(lon[coto],lat[coto],off[coto],nomefi=cartout2+'offset'+lab+'.eps',polo=polo,title='Offset [W/m2/nm/sr]',aur_model=aurm)


    jirfu.stereoplot(lon[coto],lat[coto],pixs.emiss_angle[coto],nomefi=cartout2+'emiss_angle'+lab+'.eps',polo=polo,title='Emission angle [deg]',aur_model=aurm)
    jirfu.stereoplot(lon[coto],lat[coto],pixs.slant_dist[coto],nomefi=cartout2+'slant_dist'+lab+'.eps',polo=polo,title='Slant distance (km)',aur_model=aurm)
    jirfu.stereoplot(lon[coto],lat[coto],pixs.solar_time[coto],nomefi=cartout2+'solar_time'+lab+'.eps',polo=polo,title='Solar time [hr]',aur_model=aurm)
    jirfu.stereoplot(lon[coto],lat[coto],pixs.incid_angle[coto],nomefi=cartout2+'incid_angle'+lab+'.eps',polo=polo,title='Incidence angle [deg]',aur_model=aurm)
    jirfu.stereoplot(lon[coto],lat[coto],pixs.phase_angle[coto],nomefi=cartout2+'phase_angle'+lab+'.eps',polo=polo,title='Phase angle [deg]',aur_model=aurm)
    jirfu.stereoplot(lon[coto],lat[coto],np.arange(len(pixs))[coto],nomefi=cartout2+'index'+lab+'.eps',polo=polo,title='Index',aur_model=aurm)


    # Alcuni scatter plots
    yy = lambda x: 100/0.4*x
    jirfu.scatter(col,temp,err_t,cond = coto,lat=pixs.pc_lat,lon=pixs.pc_lon,ylim=(700,1200),xlim=(0,1e13),pdf=pdf_scatt_1,polo=polo,
                  nome=cartout2+'SP_col_temp_errt'+lab+'.eps',xlabel='H3+ col (cm-2)',ylabel='Temp (K)',title='color = Temp. err (K) : FILT {}'.format(filt))
    jirfu.scatter(err_c,err_t,temp,cond=coto,lat=pixs.pc_lat,lon=pixs.pc_lon,ylim=(0,200),xlim=(0,1.5e12),clim=(700,tmax),pdf=pdf_scatt_2,polo=polo,
                  nome=cartout2+'SP_errc_errt_temp'+lab+'.eps',xlabel='Error H3+ col (cm-2)',ylabel='Error Temp (K)',title='color = Temp (K) : FILT {}'.format(filt))
    jirfu.scatter(err_c/col,err_t,temp,cond = coto,lat=pixs.pc_lat,lon=pixs.pc_lon,ylim=(0,200),xlim=(0,1),clim=(700,tmax),pdf=pdf_scatt_3,polo=polo,
                  nome=cartout2+'SP_errcRel_errt_temp'+lab+'.eps',xlabel='Rel. err. H3+ col',ylabel='Error Temp (K)',title='color = Temp (K) : FILT {}'.format(filt))
    jirfu.scatter(err_c/col,err_t,col,cond = coto,lat=pixs.pc_lat,lon=pixs.pc_lon,ylim=(0,200),xlim=(0,1),clim=(0,8e12),pdf=pdf_scatt_3,polo=polo,
                  nome=cartout2+'SP_errcRel_errt_temp'+lab+'.eps',xlabel='Rel. err. H3+ col',ylabel='Error Temp (K)',title='color = H3+ col. (cm-2) : FILT {}'.format(filt))

    if filt == 1:
        gi = err_t < yy(err_c/col)
        jirfu.scatter(err_c/col,err_t,temp,cond = coto & gi,lat=pixs.pc_lat,lon=pixs.pc_lon,ylim=(0,200),xlim=(0,1),clim=(700,tmax),polo=polo,pdf=pdf_scatt_3,
                      nome=cartout2+'SP_errcRel_errt_temp_giu'+lab+'.eps',xlabel='Rel. err. H3+ col (cm-2)',ylabel='Error temp (K)',title='color = Temp (K) : FILT {}'.format(filt))
        gi = err_t > yy(err_c/col)
        jirfu.scatter(err_c/col,err_t,temp,cond = coto & gi,lat=pixs.pc_lat,lon=pixs.pc_lon,ylim=(0,200),xlim=(0,1),clim=(700,tmax),polo=polo,pdf=pdf_scatt_3,
                      nome=cartout2+'SP_errcRel_errt_temp_su'+lab+'.eps',xlabel='Rel. err. H3+ col (cm-2)',ylabel='Error temp (K)',title='color = Temp (K) : FILT {}'.format(filt))

    # stampo i risultati
    fi = open(cart+'tutti'+lab+'.dat','w')
    i=0
    fu='{:05d}'+'{:>28s}'+2*'{:4d}'+5*'{:7.2f}'+'{:10.2e}'+3*'{:11.3e}'+2*'{:9.2f}'+'{:7.2f}'+'{:9.2e}'+'{:11.3e}'+'{:7.2f}'+'\n'
    fun='{:>5s}'+'{:>28s}'+2*'{:>4s}'+5*'{:>7s}'+'{:>10s}'+3*'{:>11s}'+2*'{:>9s}'+'{:>7s}'+'{:>9s}'+'{:>11s}'+'{:>7s}'+'\n'
    fi.write('Results for cart {:40s}\n'.format(cart))
    fi.write('Legend: num -> id number, cubo -> obs. cube, i -> line, j -> sample, lon -> planetocentric longitude, '
             'lat -> planetocentric latitude, emi -> emission angle, sza -> sza, stime -> local solar time, '
             'dist -> slant distance [m], col -> retrieved h3p column [cm-2], err_c -> retrieval error, '
             'col_c -> column corrected for emission angle, temp -> retrieved h3p temp [K], err_t -> retr. error, '
             'chi -> chi squared, off -> offset [W/m^2/sr/nm], ch4 -> retrieved ch4 column [cm-2], wls -> wl shift (nm) \n')
    fi.write('{:1s}\n'.format(' '))
    fi.write(fun.format('num','cubo','i','j','lon','lat','emi','sza','stime','dist','col','err_c','col_c','temp',
                        'err_t','chi','off','ch4','wls'))
    fi.write('{:1s}\n'.format('#'))
    for pi, ind in zip(pixs, indx):
        fi.write(fu.format(ind,pi['cubo'],pi['line'],pi['sample'],pi['pc_lon'],pi['pc_lat'],pi['emiss_angle'],
                           pi['incid_angle'],pi['solar_time'],pi['slant_dist'],
                           col[i],err_c[i],col_c[i],temp[i],err_t[i],chi[i],off[i],ch4[i],shi[i])) # num,lon,lat,emi,sza,loc,dist,))
    fi.close()

pdf_col.close()
pdf_temp.close()
pdf_scatt_1.close()
pdf_scatt_2.close()
pdf_scatt_3.close()
