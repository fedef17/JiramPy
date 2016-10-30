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

cart = '/home/fede/Scrivania/Jiram/DATA/JM0003_all/aur_N_nadir_70_CH4/'
cartres = cart
cartout = cart + 'MAP_ok2/'
cart80 = '/home/fede/Scrivania/Jiram/DATA/JM0003_all/aur_N_nadir_80_CH4/'
polo = 'N'
lch4 = True
l500 = True  # così fa tutto con le geometrie a 500 km
lN80 = True #True # così aggiungo le nuove misure fino a 80 gradi di emission angle

emissMAX = 75.0

###### LEGGE ########################################################################

pixs = pickle.load(open(cart+'all_pixs.pic','r'))

if l500:
    # si va a prendere le geometrie per 500 km
    p500 = pickle.load(open(cart+'all_pixs_500.pic','r'))
    pixs500 = p500.view(np.recarray)
    print(len(pixs),len(pixs500))
    print(pixs[1139])
    print(pixs500[1139])

ii,col1, err_col = jirfu.read_res_jir_3(cartres+'CD-H3p.dat')
iit,temp1, err_temp = jirfu.read_res_jir_3(cartres+'VT-H3p.dat')
iic,chi1 = jirfu.read_res_jir_4(cartres+'chisq.dat')
iio,off1 = jirfu.read_res_jir_4(cartres+'offset_ok.dat')
iis,shi1 = jirfu.read_res_jir_4(cartres+'shift_ok.dat')
if lch4: ii4,ch41, err_ch4 = jirfu.read_res_jir_3(cartres+'CD-CH4.dat')

col = np.zeros(len(pixs))
err_c = np.zeros(len(pixs))
temp = np.zeros(len(pixs))
err_t = np.zeros(len(pixs))
chi = np.zeros(len(pixs))
off = np.zeros(len(pixs))
shi = np.zeros(len(pixs))
if lch4: ch4 = np.zeros(len(pixs))

col[ii] = col1
err_c[ii] = err_col
temp[iit] = temp1
err_t[iit] = err_temp
chi[iic] = chi1
chi[chi == 0.0] = np.nan
off[iio] = off1
shi[iis] = shi1
if lch4: ch4[ii4] = ch41

if lN80:
    p80 = pickle.load(open(cart80+'all_pixs.pic','r'))
    ii,col1, err_col = jirfu.read_res_jir_3(cart80+'CD-H3p.dat')
    iit,temp1, err_temp = jirfu.read_res_jir_3(cart80+'VT-H3p.dat')
    iic,chi1 = jirfu.read_res_jir_4(cart80+'chisq.dat')
    iio,off1 = jirfu.read_res_jir_4(cart80+'offset_ok.dat')
    iis,shi1 = jirfu.read_res_jir_4(cart80+'shift_ok.dat')
    if lch4: ii4,ch41, err_ch4 = jirfu.read_res_jir_3(cart80+'CD-CH4.dat')

    col80 = np.zeros(len(p80))
    err_c80 = np.zeros(len(p80))
    temp80 = np.zeros(len(p80))
    err_t80 = np.zeros(len(p80))
    chi80 = np.zeros(len(p80))
    off80 = np.zeros(len(p80))
    shi80 = np.zeros(len(p80))
    if lch4: ch480 = np.zeros(len(p80))

    col80[ii] = col1
    err_c80[ii] = err_col
    temp80[iit] = temp1
    err_t80[iit] = err_temp
    chi80[iic] = chi1
    chi80[chi80 == 0.0] = np.nan
    off80[iio] = off1
    shi80[iis] = shi1
    if lch4: ch480[ii4] = ch41

    print(type(pixs),len(pixs))
    p80 = p80.view(np.recarray)
    co80 = (p80.emiss_angle < emissMAX)
    p80 = p80[co80]
    pixs = np.append(pixs,p80)
    pixs500 = np.append(pixs500,p80)
    print(type(pixs),len(pixs))
    pixs = pixs.view(np.recarray)
    pixs500 = pixs500.view(np.recarray)
    print(type(pixs),len(pixs))
    col = np.append(col,col80[co80])
    err_c = np.append(err_c,err_c80[co80])
    temp = np.append(temp,temp80[co80])
    err_t = np.append(err_t,err_t80[co80])
    chi = np.append(chi,chi80[co80])
    off = np.append(off,off80[co80])
    shi = np.append(shi,shi80[co80])
    if lch4: ch4 = np.append(ch4,ch480[co80])

print(type(pixs.cubo))
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

cond2 = ((chi > 4) | (col <= 0.0) | (off < 0.0))

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
# iii = np.a
# range(len(pixs))
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

mcol = 3e12
npix = 40
aurm = 'both'


##################################################################

# Condizione:

pdf_col = PdfPages(cartout+'MAPPE_h3p_col.pdf')
pdf_temp = PdfPages(cartout+'MAPPE_h3p_temp.pdf')
pdf_scatt_1 = PdfPages(cartout+'MAPPE_scatt_col_temp.pdf')
pdf_scatt_2 = PdfPages(cartout+'MAPPE_scatt_errc_errt.pdf')
pdf_scatt_3 = PdfPages(cartout+'MAPPE_scatt_errcRel_errt.pdf')

for filt in [1,2,3,4]:
    print('TOT: {}'.format(len(pixs)))

    if filt == 1:
        coto = (~np.isnan(col))
        lab = '_fil1'
        print('OK'+lab+': {}'.format(len(pixs[coto])))
    elif filt == 2:
        coto = (~np.isnan(col)) & (err_c/col < 0.3) & (err_t < 100)
        lab = '_fil2'
        print('OK'+lab+': {}'.format(len(pixs[coto])))
    elif filt == 3:
        coto = (~np.isnan(col)) & (err_c < 0.5e12) & (err_t < 100)
        lab = '_fil3'
        print('OK'+lab+': {}'.format(len(pixs[coto])))
    elif filt == 4:
        coto = (~np.isnan(col)) & (col > 1.5e12) & (err_t < 100)
        lab = '_fil4'
        print('OK'+lab+': {}'.format(len(pixs[coto])))
    ##################################################################

    jirfu.stereomap2(lon[coto],lat[coto],col_c[coto],nomefi=cartout+'MAP_h3p_col'+lab+'.eps',title='H3+ column [mol/cm2] : FILT {}'.format(filt),polo=polo,min=0.2e12,max=mcol, show=False, xres=npix, image=False, aur_model=aurm,pdf=pdf_col)
    jirfu.stereomap2(lon[coto],lat[coto],err_cc[coto]/col_c[coto],nomefi=cartout+'ERR_h3p_col'+lab+'.eps',title='Relative error H3+ column [mol/cm2] : FILT {}'.format(filt),polo=polo,min=0.,max=1., show=False,xres=npix)
    jirfu.stereomap2(lon[coto],lat[coto],col_c[coto],nomefi=cartout+'Npoints_ok'+lab+'.eps',title='Number of measurements : FILT {}'.format(filt),polo=polo,show=False,xres=50,npoints=True,aur_model=aurm)#,addpoints=True)

    co = pixs.emiss_angle < 60
    jirfu.stereomap2(lon[co & coto],lat[co & coto],col_c[co & coto],nomefi=cartout+'MAP_h3p_col'+lab+'_emissm60.eps',title='H3+ column [mol/cm2] : emiss < 60',polo=polo,min=0.2e12,max=mcol, show=False, xres=npix,aur_model=aurm,addpoints=True)
    #jirfu.stereoplot(lon[co & coto],lat[co & coto],col_c[co & coto],nomefi=cartout+'PLOT_h3p_col'+lab+'_emissm60.eps',title='H3+ column [mol/cm2] : emiss < 60 ',polo=polo,min=0.2e12,max=mcol,show=False,aur_model=aurm)
    co = pixs.emiss_angle > 60
    jirfu.stereomap2(lon[co & coto],lat[co & coto],col_c[co & coto],nomefi=cartout+'MAP_h3p_col'+lab+'_emissp60.eps',title='H3+ column [mol/cm2] : emiss > 60',polo=polo,min=0.2e12,max=mcol, show=False, xres=npix,aur_model=aurm,addpoints=True)
    #jirfu.stereoplot(lon[co & coto],lat[co & coto],col_c[co & coto],nomefi=cartout+'PLOT_h3p_col'+lab+'_emissp60.eps',title='H3+ column [mol/cm2] : emiss > 60',polo=polo,min=0.2e12,max=mcol,show=False,aur_model=aurm)

    cond4 = (pixs.alt_surf < 4e8)
    nome = cartout+'MAP_h3p_col_close'+lab+'.eps'
    jirfu.stereomap2(lon[cond4 & coto],lat[cond4 & coto],col_c[cond4 & coto],nomefi=nome,title='H3+ column [mol/cm2] -> alt_surf < 4e5 km',polo=polo,min=0.2e12,max=mcol,show=False,xres=npix,aur_model=aurm)
    cond4 = (pixs.alt_surf > 4e8)
    nome = cartout+'MAP_h3p_col_far'+lab+'.eps'
    jirfu.stereomap2(lon[cond4 & coto],lat[cond4 & coto],col_c[cond4 & coto],nomefi=nome,title='H3+ column [mol/cm2] -> alt_surf > 4e5 km',polo=polo,min=0.2e12,max=mcol,show=False,xres=npix,aur_model=aurm)
    cond4 = (pixs.solar_time < 12) & (pixs.solar_time > 6)
    nome = cartout+'MAP_h3p_col_am'+lab+'.eps'
    jirfu.stereomap2(lon[cond4 & coto],lat[cond4 & coto],col_c[cond4 & coto],nomefi=nome,title='H3+ column [mol/cm2] -> 6h < solar time < 12h',polo=polo,min=0.2e12,max=mcol,show=False,xres=npix,aur_model=aurm)
    cond4 = (pixs.solar_time > 12) & (pixs.solar_time < 18)
    nome = cartout+'MAP_h3p_col_pm'+lab+'.eps'
    jirfu.stereomap2(lon[cond4 & coto],lat[cond4 & coto],col_c[cond4 & coto],nomefi=nome,title='H3+ column [mol/cm2] -> 12h < solar time < 18h',polo=polo,min=0.2e12,max=mcol,show=False,xres=npix,aur_model=aurm)
    cond4 = (pixs.solar_time < 24) & (pixs.solar_time > 18)
    nome = cartout+'MAP_h3p_col_ni1'+lab+'.eps'
    jirfu.stereomap2(lon[cond4 & coto],lat[cond4 & coto],col_c[cond4 & coto],nomefi=nome,title='H3+ column [mol/cm2] -> 18h < solar time < 24h',polo=polo,min=0.2e12,max=mcol,show=False,xres=npix,aur_model=aurm)
    cond4 = (pixs.solar_time < 6) & (pixs.solar_time > 0)
    nome = cartout+'MAP_h3p_col_ni2'+lab+'.eps'
    jirfu.stereomap2(lon[cond4 & coto],lat[cond4 & coto],col_c[cond4 & coto],nomefi=nome,title='H3+ column [mol/cm2] -> 0h < solar time < 6h',polo=polo,min=0.2e12,max=mcol,show=False,xres=npix,aur_model=aurm)

    # giu = err_t < 100
    lab2 = ''
    jirfu.stereomap2(lon[coto],lat[coto],temp[coto],nomefi=cartout+'MAP_h3p_temp'+lab2+lab+'.eps',title='H3+ temperature [K] : FILT {}'.format(filt), polo=polo,min=700.,max=1100, show=False,xres=npix,aur_model=aurm,pdf=pdf_temp)
    jirfu.stereomap2(lon[coto],lat[coto],err_t[coto],nomefi=cartout+'ERR_h3p_temp'+lab2+lab+'.eps',title='Error H3+ temperature [K] : FILT {}'.format(filt), polo=polo, min=0.,max=200., show=False,xres=npix,aur_model=aurm)
    # giu = err_t < 80
    # lab2 = '_err80K'
    # jirfu.stereomap2(lon[giu],lat[giu],temp[giu],nomefi=cartout+'MAP_h3p_temp'+lab2+lab+'.eps',title='H3+ temperature [K]', polo=polo,min=700.,max=1100, show=False,xres=npix,aur_model=aurm)
    # jirfu.stereomap2(lon[giu],lat[giu],err_t[giu],nomefi=cartout+'ERR_h3p_temp'+lab2+lab+'.eps',title='Error H3+ temperature [K]', polo=polo, min=0.,max=200., show=False,xres=npix,aur_model=aurm)
    # giu = err_t < 50
    # lab2 = '_err50K'
    # jirfu.stereomap2(lon[giu],lat[giu],temp[giu],nomefi=cartout+'MAP_h3p_temp'+lab2+lab+'.eps',title='H3+ temperature [K]', polo=polo,min=700.,max=1100, show=False,xres=npix,aur_model=aurm)
    # jirfu.stereomap2(lon[giu],lat[giu],err_t[giu],nomefi=cartout+'ERR_h3p_temp'+lab2+lab+'.eps',title='Error H3+ temperature [K]', polo=polo, min=0.,max=200., show=False,xres=npix,aur_model=aurm)
    # condt = (col > 0e12)
    # lab2 = '_nofilt'
    # jirfu.stereomap2(lon[condt],lat[condt],temp[condt],nomefi=cartout+'MAP_h3p_temp'+lab2+lab+'.eps',title='H3+ temperature [K]', polo=polo,min=700.,max=1100, show=False,xres=npix,aur_model=aurm)
    # jirfu.stereomap2(lon[condt],lat[condt],err_t[condt],nomefi=cartout+'ERR_h3p_temp'+lab2+lab+'.eps',title='Error H3+ temperature [K]', polo=polo, min=0.,max=200., show=False,xres=npix,aur_model=aurm)

    jirfu.stereomap2(lon[coto],lat[coto],ch4[coto],nomefi=cartout+'MAP_ch4_col'+lab+'.eps',title='Effective CH4 column [mol/cm2]',polo=polo, show=False,xres=npix,aur_model=aurm)


    ## l'emissione integrata:
    condpi = ~np.isnan(pixs.pc_lon)
    zu = pixs[condpi & coto]
    integr = np.zeros(len(zu))
    for ww,spp,of,i in zip(zu.wl,zu.spe,off[condpi & coto],range(len(zu))):
        mask = jirfu.findspi(ww,spp)
        integr[i] = jirfu.integr_h3p(ww,spp*mask,of)

    integr=integr*np.cos(np.pi*zu.emiss_angle/180.0)
    jirfu.stereomap2(lon,lat,integr,nomefi=cartout+'MAP_integr_h3p_conslit'+lab+'.eps',title='H3+ integrated intensity [W/m2/sr]',
                    polo=polo,form=True,addpoints=True,show=False,xres=npix,aur_model=aurm)
    jirfu.stereomap2(lon,lat,integr,nomefi=cartout+'MAP_integr_h3p'+lab+'.eps',title='H3+ integrated intensity [W/m2/sr]',
                    polo=polo,form=True,addpoints=False,show=False,xres=npix,aur_model=aurm)

    # pl.close()
    # pl.scatter(temp[condpi],integr,c=err_t[condpi],s=2,)
    # pl.colorbar()
    # pl.show()
    #sys.exit()

    # i plot punto per punto
    jirfu.stereoplot(lon[coto],lat[coto],col_c[coto],nomefi=cartout+'h3p_col'+lab+'.eps',title='H3+ column [mol/cm2]',polo=polo,min=0.2e12,max=mcol,aur_model=aurm)
    jirfu.stereoplot(lon[coto],lat[coto],col[coto],nomefi=cartout+'h3p_col_nocorr'+lab+'.eps',title='H3+ column [mol/cm2] (non corrected)',polo=polo,min=0.2e12,max=8e12,aur_model=aurm)
    jirfu.stereoplot(lon[coto],lat[coto],temp[coto],nomefi=cartout+'h3p_temp'+lab+'.eps',title='H3+ temperature [K]',polo=polo,min=700.,max=1100,aur_model=aurm)
    if lch4: jirfu.stereoplot(lon[coto],lat[coto],ch4[coto],nomefi=cartout+'ch4_col'+lab+'.eps',title='Effective CH4 column [mol/cm2]',polo=polo,aur_model=aurm)
    jirfu.stereoplot(lon[coto],lat[coto],off[coto],nomefi=cartout+'offset'+lab+'.eps',polo=polo,title='Offset [W/m2/nm/sr]',aur_model=aurm)


    jirfu.stereoplot(lon[coto],lat[coto],pixs.emiss_angle[coto],nomefi=cartout+'emiss_angle'+lab+'.eps',polo=polo,title='Emission angle [deg]',aur_model=aurm)
    jirfu.stereoplot(lon[coto],lat[coto],pixs.slant_dist[coto],nomefi=cartout+'slant_dist'+lab+'.eps',polo=polo,title='Slant distance (km)',aur_model=aurm)
    jirfu.stereoplot(lon[coto],lat[coto],pixs.solar_time[coto],nomefi=cartout+'solar_time'+lab+'.eps',polo=polo,title='Solar time [hr]',aur_model=aurm)
    jirfu.stereoplot(lon[coto],lat[coto],pixs.incid_angle[coto],nomefi=cartout+'incid_angle'+lab+'.eps',polo=polo,title='Incidence angle [deg]',aur_model=aurm)
    jirfu.stereoplot(lon[coto],lat[coto],pixs.phase_angle[coto],nomefi=cartout+'phase_angle'+lab+'.eps',polo=polo,title='Phase angle [deg]',aur_model=aurm)
    jirfu.stereoplot(lon[coto],lat[coto],np.arange(len(pixs))[coto],nomefi=cartout+'index'+lab+'.eps',polo=polo,title='Index',aur_model=aurm)


    # Alcuni scatter plots
    yy = lambda x: 100/0.4*x
    jirfu.scatter(col,temp,err_t,cond = coto,lat=pixs.pc_lat,lon=pixs.pc_lon,ylim=(700,1200),xlim=(0,1e13),pdf=pdf_scatt_1,polo=polo,
                  nome=cartout+'SP_col_temp_errt'+lab+'.eps',xlabel='H3+ col (cm-2)',ylabel='Temp (K)',title='color = Temp. err (K) : FILT {}'.format(filt))
    jirfu.scatter(err_c,err_t,temp,cond=coto,lat=pixs.pc_lat,lon=pixs.pc_lon,ylim=(0,200),xlim=(0,1.5e12),clim=(700,1100),pdf=pdf_scatt_2,polo=polo,
                  nome=cartout+'SP_errc_errt_temp'+lab+'.eps',xlabel='Error H3+ col (cm-2)',ylabel='Error Temp (K)',title='color = Temp (K) : FILT {}'.format(filt))
    jirfu.scatter(err_c/col,err_t,temp,cond = coto,lat=pixs.pc_lat,lon=pixs.pc_lon,ylim=(0,200),xlim=(0,1),clim=(700,1100),pdf=pdf_scatt_3,polo=polo,
                  nome=cartout+'SP_errcRel_errt_temp'+lab+'.eps',xlabel='Rel. err. H3+ col',ylabel='Error Temp (K)',title='color = Temp (K) : FILT {}'.format(filt))
    jirfu.scatter(err_c/col,err_t,col,cond = coto,lat=pixs.pc_lat,lon=pixs.pc_lon,ylim=(0,200),xlim=(0,1),clim=(0,8e12),pdf=pdf_scatt_3,polo=polo,
                  nome=cartout+'SP_errcRel_errt_temp'+lab+'.eps',xlabel='Rel. err. H3+ col',ylabel='Error Temp (K)',title='color = H3+ col. (cm-2) : FILT {}'.format(filt))

    if filt == 1:
        gi = err_t < yy(err_c/col)
        jirfu.scatter(err_c/col,err_t,temp,cond = coto & gi,lat=pixs.pc_lat,lon=pixs.pc_lon,ylim=(0,200),xlim=(0,1),clim=(700,1100),polo=polo,pdf=pdf_scatt_3,
                      nome=cartout+'SP_errcRel_errt_temp_giu'+lab+'.eps',xlabel='Rel. err. H3+ col (cm-2)',ylabel='Error temp (K)',title='color = Temp (K) : FILT {}'.format(filt))
        gi = err_t > yy(err_c/col)
        jirfu.scatter(err_c/col,err_t,temp,cond = coto & gi,lat=pixs.pc_lat,lon=pixs.pc_lon,ylim=(0,200),xlim=(0,1),clim=(700,1100),polo=polo,pdf=pdf_scatt_3,
                      nome=cartout+'SP_errcRel_errt_temp_su'+lab+'.eps',xlabel='Rel. err. H3+ col (cm-2)',ylabel='Error temp (K)',title='color = Temp (K) : FILT {}'.format(filt))

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