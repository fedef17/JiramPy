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
import matplotlib.gridspec as gridspec

################# SETUP ##################

# File from which the input variables are read:
input_file = '/home/fedefab/Scrivania/Research/Dotto/Git/JiramPy/input_mappe.in'

# Tipi di filtering da applicare ai risultati. Lista di valori da 0 a 4. 0 = nessun filtro, 1 = filtro articoli 2017. Per farli tutti mettere [0,1,2,3,4]
list_filts = [0,1]

##################################################

keys = 'cartres cartout picfile polo lch4 lshi maxchi mcol mincol tmax tmin npix aurm emissMAX cstep tstep errmax_T maxint maxch4 maxch4_c stepint stepch4 stepch4_c lim_ratio_max lim_ratio_min step_ratio npix_ratio'
keys = keys.split()
types = [str,str,str,str,bool,bool,float,float,float,float,float,int,str,float,float,float,float,float,float,float,float,float,float,float,float,float,int]

values = sbm.read_inputs(input_file,keys,itype=types,verbose=True)

cartres = values['cartres']
picfile = values['picfile']
cartout = values['cartout']
if not os.path.exists(cartout): os.makedirs(cartout)

polo = values['polo']
lch4 = values['lch4']
lshi = values['lshi']

maxchi = values['maxchi']
mcol = values['mcol']
mincol = values['mincol']
tmax = values['tmax']
tmin = values['tmin']
npix = values['npix']
aurm = values['aurm']
emissMAX = values['emissMAX']
cstep = values['cstep']
tstep = values['tstep']
errmax_T = values['errmax_T']
maxint = values['maxint']
maxch4 = values['maxch4']
maxch4_c = values['maxch4_c']
stepint = values['stepint']
stepch4 = values['stepch4']
stepch4_c = values['stepch4_c']
lim_ratio_max = values['lim_ratio_max']
lim_ratio_min = values['lim_ratio_min']
step_ratio = values['step_ratio']
npix_ratio = values['npix_ratio']

cbarlabcol = r'$H_3^+$ column ($\times 10^{{{}}}$ cm$^{{-2}}$)'
cbarlaberrcol = 'Relative error'
cbarlabtemp = r'$H_3^+$ temperature (K)'
cbarlaberrtemp = r'Error on $H_3^+$ temperature (K)'
cbarlabch4 = r'CH$_4$ column ($\times 10^{{{}}}$ cm$^{{-2}}$)'
cbarlabint = r'Integrated intensity ($\times 10^{{{}}}$ $W {{m}}^{{-2}} {{sr}}^{{-1}}$)'

titleintegr = 'Integrated intensity'
titlecol = r'$H_3^+$ column'
titleerrcol = r'Relative error on $H_3^+$ column'
titletemp = r'$H_3^+$ temperature'
titleerrtemp = r'Error on $H_3^+$ temperature'
titlech4 = r'CH$_4$ column'

titleintegr = ''
titlecol = ''
titleerrcol = ''
titletemp = ''
titleerrtemp = ''
titlech4 = ''

###################################

pixs = pickle.load(open(picfile,'r'))

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

cubs = np.unique(pixs.cubo)
i=0
for cu in cubs:
    if not np.all(np.isnan(pixs[pixs.cubo == cu].pc_lon)):
        em = np.mean(pixs[pixs.cubo == cu].emiss_angle)
        so = np.mean(pixs[pixs.cubo == cu].solar_time)
        n = len(pixs[pixs.cubo == cu])
        print('{}  {}  {:5.1f}  {:5.1f} {}'.format(i,cu,em,so,n))

######## fine letture ###################
# Faccio le mappe:

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

indx = np.arange(len(pixs))

nonan = (~np.isnan(col)) & (~np.isnan(pixs.pc_lon))
col = col[nonan]
temp = temp[nonan]
chi = chi[nonan]
off = off[nonan]
err_t = err_t[nonan]
err_c = err_c[nonan]
ch4 = ch4[nonan]
shi = shi[nonan]
pixs = pixs[nonan]
indx = indx[nonan]

print('CONTROLLOOOOO', np.sum(pixs.incid_angle == 0))

col_c = col*np.cos(np.pi*pixs.emiss_angle/180.0)
err_cc = err_c*np.cos(np.pi*pixs.emiss_angle/180.0)
ch4_c = ch4*np.cos(np.pi*pixs.emiss_angle/180.0)

lon = pixs.pc_lon
lat = pixs.pc_lat
lab = '_ok1'

edges=[]
for pi in zip(pixs.pc_lon_1, pixs.pc_lon_2, pixs.pc_lon_3, pixs.pc_lon_4, pixs.pc_lat_1, pixs.pc_lat_2, pixs.pc_lat_3, pixs.pc_lat_4):
    vec = [[pi[i],pi[i+4]] for i in range(4)]
    edges.append(vec)

edges = np.array(edges)
print('EDGES')
print(len(edges),len(pixs))

integr = np.zeros(len(pixs))
ratio_32_35 = np.zeros(len(pixs))
ratio_32_37 = np.zeros(len(pixs))
ratio_35_37 = np.zeros(len(pixs))
erratio_32_35 = np.zeros(len(pixs))
erratio_32_37 = np.zeros(len(pixs))
erratio_35_37 = np.zeros(len(pixs))
for ww,spp,of,i in zip(pixs.wl,pixs.spe,off,range(len(pixs))):
    nano = np.isnan(spp)
    spp[nano] = 0.0
    mask = jirfu.findspi(ww,spp)
    integr[i] = jirfu.integr_h3p(ww,spp*mask,of)
    ratio_32_35[i], erratio_32_35[i] = jirfu.ratio32_35(ww, spp*mask, of)
    ratio_32_37[i], erratio_32_37[i] = jirfu.ratio32_37(ww, spp*mask, of)
    ratio_35_37[i], erratio_35_37[i] = jirfu.ratio35_37(ww, spp*mask, of)

integr=integr*np.cos(np.pi*pixs.emiss_angle/180.0)

##################################################################

# Condizione:

pdf_col = PdfPages(cartout+'MAPPE_h3p_col.pdf')
pdf_temp = PdfPages(cartout+'MAPPE_h3p_temp.pdf')
pdf_scatt_1 = PdfPages(cartout+'MAPPE_scatt_col_temp.pdf')
pdf_scatt_2 = PdfPages(cartout+'MAPPE_scatt_errc_errt.pdf')
pdf_scatt_3 = PdfPages(cartout+'MAPPE_scatt_errcRel_errt.pdf')

#for filt in [0,1,2,3,4]:
for filt in list_filts:
    print('TOT: {}'.format(len(pixs)))

    if filt == 0:
        coto = (~np.isnan(col))
        lab = '_nofil'
        print('OK'+lab+': {}'.format(len(pixs[coto])))
        cartout2 = cartout + lab + '/'
        if not os.path.exists(cartout2): os.makedirs(cartout2)
    elif filt == 1:
        coto = (~np.isnan(col)) & (err_t < errmax_T)
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

# I RATIOS

    if polo == 'S':
        blats = [-60.0]*4
    else:
        blats = [55.0]*4

    npix2 = npix_ratio
    cmappa = 'jet'
    # titleratio = r'Band intensity ratio: I({}$\,\mu$m)/I({}$\,\mu$m)'
    # jirfu.stereomap2(lon[coto],lat[coto],ratio_32_35[coto],nomefi=cartout2+'RATIO_MAP_32_35_'+lab+'.pdf',title=titleratio.format(3.2,3.5),step=0.05,polo=polo,cbarlabel='',cbarform='%.1f',minu=0.5,maxu=1.0,addpoints=False,show=False,xres=npix2,aur_model=aurm, cmap = cmappa, blats = blats)
    # jirfu.stereomap2(lon[coto],lat[coto],ratio_32_37[coto],nomefi=cartout2+'RATIO_MAP_32_37_'+lab+'.pdf',title=titleratio.format(3.2,3.7),step=0.05,polo=polo,cbarlabel='',cbarform='%.1f',minu=0.5,maxu=1.0,addpoints=False,show=False,xres=npix2,aur_model=aurm, cmap = cmappa, blats = blats)
    # jirfu.stereomap2(lon[coto],lat[coto],ratio_35_37[coto],nomefi=cartout2+'RATIO_MAP_35_37_'+lab+'_VIP4.pdf',title=titleratio.format(3.5,3.7),step=step_ratio,polo=polo,cbarlabel='',cbarform='%.2f',minu=lim_ratio_min,maxu=lim_ratio_max,addpoints=False,show=False,xres=npix2,aur_model='VIP4', cmap = cmappa, blats = blats, print_values = True)
    # jirfu.stereomap2(lon[coto],lat[coto],ratio_35_37[coto],nomefi=cartout2+'RATIO_MAP_35_37_'+lab+'_stat.pdf',title=titleratio.format(3.5,3.7),step=step_ratio,polo=polo,cbarlabel='',cbarform='%.2f',minu=lim_ratio_min,maxu=lim_ratio_max,addpoints=False,show=False,xres=npix2,aur_model='stat', cmap = cmappa, blats = blats, print_values = True)
    #
    # titleratio = r'Error on band intensity ratio: I({}$\,\mu$m)/I({}$\,\mu$m)'
    # jirfu.stereomap2(lon[coto],lat[coto],erratio_32_35[coto],nomefi=cartout2+'ERRATIO_MAP_32_35_'+lab+'.pdf',title=titleratio.format(3.2,3.5),step=0.05,polo=polo,cbarlabel='',cbarform='%.1f',minu=0.0,maxu=0.5,addpoints=False,show=False,xres=npix2,aur_model=aurm, cmap = cmappa, blats = blats)
    # jirfu.stereomap2(lon[coto],lat[coto],erratio_32_37[coto],nomefi=cartout2+'ERRATIO_MAP_32_37_'+lab+'.pdf',title=titleratio.format(3.2,3.7),step=0.05,polo=polo,cbarlabel='',cbarform='%.1f',minu=0.0,maxu=0.5,addpoints=False,show=False,xres=npix2,aur_model=aurm, cmap = cmappa, blats = blats)
    # jirfu.stereomap2(lon[coto],lat[coto],erratio_35_37[coto],nomefi=cartout2+'ERRATIO_MAP_35_37_'+lab+'_VIP4.pdf',title=titleratio.format(3.5,3.7),step=0.05,polo=polo,cbarlabel='',cbarform='%.1f',minu=0.0,maxu=1.0,addpoints=False,show=False,xres=npix2,aur_model='VIP4', cmap = cmappa, blats = blats, print_values = True)
    # jirfu.stereomap2(lon[coto],lat[coto],erratio_35_37[coto],nomefi=cartout2+'ERRATIO_MAP_35_37_'+lab+'_stat.pdf',title=titleratio.format(3.5,3.7),step=0.05,polo=polo,cbarlabel='',cbarform='%.1f',minu=0.0,maxu=1.0,addpoints=False,show=False,xres=npix2,aur_model='stat', cmap = cmappa, blats = blats, print_values = True)


    fig = pl.figure(figsize=(16, 6), dpi=150)
    gs = gridspec.GridSpec(1, 2)
    gs.update(hspace=0.25,wspace=0.25)
    ax = pl.subplot(gs[0])
    jirfu.stereomap2(lon[coto],lat[coto],col_c[coto],nomefi=cartout2+'MAP_h3p_col'+lab+'.pdf',title=titlecol,polo=polo,cbarlabel=cbarlabcol,cbarmult=12,cbarform='%.1f',step=cstep,minu=mincol,maxu=mcol, show=False, xres=npix, image=False, aur_model=aurm,pdf=pdf_col,axu=ax)
    ax = pl.subplot(gs[1])
    jirfu.stereomap2(lon[coto],lat[coto],err_cc[coto]/col_c[coto],nomefi=cartout2+'ERR_h3p_col'+lab+'.pdf',title=titleerrcol,cbarlabel=cbarlaberrcol,polo=polo,show=False,xres=npix,minu=0.,maxu=0.6,step=0.05,axu=ax, aur_model=aurm)
    fig.savefig(cartout2+'MAPERR_col'+lab+'.pdf', format='pdf', dpi=150)
    pl.close(fig)

    # jirfu.stereomap2(lon[coto],lat[coto],col_c[coto],nomefi=cartout2+'Npoints_ok'+lab+'.eps',title='Number of measurements',polo=polo,show=False,xres=npix,npoints=True,aur_model=aurm)#,addpoints=True)
    #
    # cond4 = (pixs.solar_time < 12) & (pixs.solar_time > 4)
    # nome = cartout2+'MAP_h3p_col_am'+lab+'.eps'
    # jirfu.stereomap2(lon[cond4 & coto],lat[cond4 & coto],col_c[cond4 & coto],nomefi=nome,title=titlecol+' 4h <--> 12h',polo=polo,cbarlabel=cbarlabcol,cbarmult=12,cbarform='%.1f',step=cstep,minu=mincol,maxu=mcol,show=False,xres=npix,aur_model=aurm)
    # cond4 = (pixs.solar_time > 12) & (pixs.solar_time < 20)
    # nome = cartout2+'MAP_h3p_col_pm'+lab+'.eps'
    # jirfu.stereomap2(lon[cond4 & coto],lat[cond4 & coto],col_c[cond4 & coto],nomefi=nome,title=titlecol+' 12h <--> 20h',polo=polo,cbarlabel=cbarlabcol,cbarmult=12,cbarform='%.1f',step=cstep,minu=mincol,maxu=mcol,show=False,xres=npix,aur_model=aurm)
    # cond4 = (pixs.solar_time > 20) | (pixs.solar_time < 4)
    # nome = cartout2+'MAP_h3p_col_night'+lab+'.eps'
    # jirfu.stereomap2(lon[cond4 & coto],lat[cond4 & coto],col_c[cond4 & coto],nomefi=nome,title=titlecol+' 20h <--> 4h',polo=polo,cbarlabel=cbarlabcol,cbarmult=12,cbarform='%.1f',step=cstep,minu=mincol,maxu=mcol,show=False,xres=npix,aur_model=aurm)

    jirfu.stereomap2(lon[coto],lat[coto],temp[coto],nomefi=cartout2+'MAP_h3p_temp'+lab+'.pdf',title=titletemp, polo=polo,cbarlabel=cbarlabtemp,cbarform='%.0f',step=tstep,minu=tmin,maxu=tmax, show=False,xres=npix,aur_model=aurm)
    jirfu.stereomap2(lon[coto],lat[coto],col_c[coto],nomefi=cartout2+'MAP_h3p_col'+lab+'.pdf',title=titlecol,polo=polo,cbarlabel=cbarlabcol,cbarmult=12,cbarform='%.1f',step=cstep,minu=mincol,maxu=mcol, show=False, xres=npix, image=False, aur_model=aurm)

    lab2 = ''
    fig = pl.figure(figsize=(16, 6), dpi=150)
    gs = gridspec.GridSpec(1, 2)
    gs.update(hspace=0.25,wspace=0.25)
    ax = pl.subplot(gs[0])
    jirfu.stereomap2(lon[coto],lat[coto],temp[coto],nomefi=cartout2+'MAP_h3p_temp'+lab2+lab+'.pdf',title=titletemp, polo=polo,cbarlabel=cbarlabtemp,cbarform='%.0f',step=tstep,minu=tmin,maxu=tmax, show=False,xres=npix,aur_model=aurm,pdf=pdf_temp,axu=ax)
    ax = pl.subplot(gs[1])
    jirfu.stereomap2(lon[coto],lat[coto],err_t[coto],nomefi=cartout2+'ERR_h3p_temp'+lab2+lab+'.pdf',title=titleerrtemp,cbarlabel=cbarlaberrtemp,cbarform='%.0f', polo=polo, show=False,xres=npix,aur_model=aurm,minu=0.,maxu=80,step=5,axu=ax)
    fig.savefig(cartout2+'MAPERR_temp'+lab+'.pdf', format='pdf', dpi=150)
    pl.close(fig)

    jirfu.stereomap2(lon[coto],lat[coto],ch4[coto],nomefi=cartout2+'MAP_ch4_col'+lab+'.pdf',title=titlech4,cbarlabel=cbarlabch4,cbarform='%.1f',cbarmult=13,polo=polo, show=False,xres=npix,aur_model=aurm,minu=0.0,maxu=maxch4,step=stepch4)
    jirfu.stereomap2(lon[coto],lat[coto],ch4_c[coto],nomefi=cartout2+'MAP_ch4_col_c'+lab+'.pdf',title=titlech4,cbarlabel=cbarlabch4,cbarform='%.1f',cbarmult=13,polo=polo, show=False,xres=npix,aur_model=aurm,minu=0.0,maxu=maxch4_c,step=stepch4_c,salta=3)
    jirfu.stereomap2(lon[coto],lat[coto],ch4_c[coto],nomefi=cartout2+'MAP_ch4_col_norm'+lab+'.pdf',title=r'Relative CH$_4$ abundance',cbarlabel=r'Relative CH$_4$ abundance',cbarform='%.1f',polo=polo, show=False,xres=npix,aur_model=aurm,minu=0.0,maxu=maxch4_c,step=stepch4_c,salta=3, normalize = True, print_values = True)

    coto2 = (~np.isnan(col)) & (err_t < 100) & (pixs.pc_lon < 30) & (pixs.pc_lon > 0)
    pl.ion()
    jirfu.scatter(pixs.solar_time,col_c,temp,cond = coto2,lat=pixs.pc_lat,lon=pixs.pc_lon,ylim=(0.5,3.e12),xlim=(0,24),pdf=pdf_scatt_1,polo=polo, nome=cartout2+'SP_colc_lon_soltim'+lab+'.pdf',ylabel=r'H$_3^+$ col (cm^{-2})',xlabel='Solar time (hr)',title='color = Temp (K)', live = True)

    # # # i plot punto per punto
    # print('CONTROLLO',len(lon),len(lat),len(col),len(coto))
    # jirfu.stereoplot(lon[coto],lat[coto],integr[coto],nomefi=cartout2+'integr_h3p_'+lab+'.pdf',title=titleintegr,step=stepint,polo=polo,cbarlabel=cbarlabint,cbarmult=-4,cbarform='%.1f',minu=0.0,maxu=maxint,show=False,aur_model=aurm,edges=edges[coto])
    jirfu.stereoplot(lon[coto],lat[coto],col_c[coto],nomefi=cartout2+'h3p_col'+lab+'.pdf',title=titlecol,polo=polo,cbarlabel=cbarlabcol,cbarmult=12,cbarform='%.1f',step=cstep,minu=mincol,maxu=mcol,aur_model=aurm,edges=edges[coto], live = True)
    sys.exit()
    # jirfu.stereoplot(lon[coto],lat[coto],col[coto],nomefi=cartout2+'h3p_col_nocorr'+lab+'.pdf',title=titlecol+' (not corrected)',polo=polo,cbarlabel=cbarlabcol,cbarmult=12,cbarform='%.1f',step=cstep,minu=mincol,maxu=8e12,aur_model=aurm,edges=edges[coto])
    # jirfu.stereoplot(lon[coto],lat[coto],temp[coto],nomefi=cartout2+'h3p_temp'+lab+'.pdf',title=titletemp,polo=polo,cbarlabel=cbarlabtemp,cbarform='%.0f',step=tstep,minu=700.,maxu=tmax,aur_model=aurm,edges=edges[coto])
    # if lch4: jirfu.stereoplot(lon[coto],lat[coto],ch4[coto],nomefi=cartout2+'ch4_col'+lab+'.pdf',title=titlech4,polo=polo,aur_model=aurm,edges=edges[coto],cbarlabel=cbarlabch4,cbarform='%.1f',cbarmult=13,minu=0.0,maxu=maxch4)
    # jirfu.stereoplot(lon[coto],lat[coto],off[coto],nomefi=cartout2+'offset'+lab+'.pdf',polo=polo,title='Offset [W/m2/nm/sr]',aur_model=aurm,edges=edges[coto])
    #
    #
    # jirfu.stereoplot(lon[coto],lat[coto],pixs.emiss_angle[coto],nomefi=cartout2+'emiss_angle'+lab+'.pdf',polo=polo,title='Emission angle [deg]',aur_model=aurm,edges=edges[coto])
    # jirfu.stereoplot(lon[coto],lat[coto],pixs.slant_dist[coto],nomefi=cartout2+'slant_dist'+lab+'.pdf',polo=polo,title='Slant distance (km)',aur_model=aurm,edges=edges[coto])
    # jirfu.stereoplot(lon[coto],lat[coto],pixs.solar_time[coto],nomefi=cartout2+'solar_time'+lab+'.pdf',polo=polo,title='Solar time [hr]',aur_model=aurm,edges=edges[coto])
    # jirfu.stereoplot(lon[coto],lat[coto],pixs.incid_angle[coto],nomefi=cartout2+'incid_angle'+lab+'.pdf',polo=polo,title='Incidence angle [deg]',aur_model=aurm,edges=edges[coto],invert_cmap=True)
    # jirfu.stereoplot(lon[coto],lat[coto],pixs.phase_angle[coto],nomefi=cartout2+'phase_angle'+lab+'.pdf',polo=polo,title='Phase angle [deg]',aur_model=aurm,edges=edges[coto])
    # #jirfu.stereoplot(lon[coto],lat[coto],np.arange(len(pixs))[coto],nomefi=cartout2+'index'+lab+'.pdf',polo=polo,title='Index',aur_model=aurm,edges=edges[coto])
    #

    jirfu.stereomap2(lon[coto],lat[coto],integr[coto],nomefi=cartout2+'MAP_integr_h3p'+lab+'.pdf',title=titleintegr,step=stepint,polo=polo,cbarlabel=cbarlabint,cbarmult=-4,cbarform='%.1f',minu=0.0,maxu=maxint,addpoints=False,show=False,xres=npix,aur_model=aurm)
    #jirfu.stereomap2(lon[coto],lat[coto],ch4[coto],nomefi=cartout2+'MAP_ch4_col'+lab+'.pdf',title=titlech4,cbarlabel=cbarlabch4,cbarform='%.1f',cbarmult=13,polo=polo, show=False,xres=npix,aur_model=aurm,minu=0.0,maxu=maxch4,step=stepch4)

    # Faccio il panel
    fig = pl.figure(figsize=(8, 6), dpi=150)
    gs = gridspec.GridSpec(2, 2)
    gs.update(hspace=0.25,wspace=0.25)
    ax = pl.subplot(gs[0])
    jirfu.stereoplot(lon[coto],lat[coto],pixs.emiss_angle[coto],nomefi=cartout2+'emiss_angle'+lab+'.pdf',polo=polo,cbarlabel='Emission angle (deg)',aur_model=aurm,edges=edges[coto],axu=ax,style='medium')
    ax = pl.subplot(gs[1])
    jirfu.stereoplot(lon[coto],lat[coto],pixs.incid_angle[coto],nomefi=cartout2+'incid_angle'+lab+'.pdf',polo=polo,cbarlabel='Incidence angle (deg)',aur_model=aurm,edges=edges[coto],invert_cmap=True,axu=ax,style='medium')
    ax = pl.subplot(gs[2])
    jirfu.stereoplot(lon[coto],lat[coto],integr[coto],nomefi=cartout2+'integr_h3p_'+lab+'.eps',title=titleintegr,step=stepint,polo=polo,cbarlabel=cbarlabint,cbarmult=-4,cbarform='%.1f',minu=0.0,maxu=maxint,show=False,aur_model=aurm,edges=edges[coto],axu=ax,style='medium')
    ax = pl.subplot(gs[3])
    jirfu.stereomap2(lon[coto],lat[coto],integr[coto],nomefi=cartout2+'MAP_integr_h3p'+lab+'.eps',title=titleintegr,step=stepint,polo=polo,cbarlabel=cbarlabint,cbarmult=-4,cbarform='%.1f',minu=0.0,maxu=maxint,addpoints=False,show=False,xres=npix,aur_model=aurm,axu=ax,style='medium')
    #pl.tight_layout()
    fig.savefig(cartout2+'panel_'+lab+'.pdf', format='pdf', dpi=150)

    jirfu.stereoplot(lon[coto],lat[coto],pixs.emiss_angle[coto],nomefi=cartout2+'emiss_angle'+lab+'.pdf',polo=polo,cbarlabel='Emission angle (deg)',aur_model=aurm,edges=edges[coto])
    jirfu.stereoplot(lon[coto],lat[coto],pixs.incid_angle[coto],nomefi=cartout2+'incid_angle'+lab+'.pdf',polo=polo,cbarlabel='Incidence angle (deg)',aur_model=aurm,edges=edges[coto],invert_cmap=True)
    jirfu.stereoplot(lon[coto],lat[coto],integr[coto],nomefi=cartout2+'integr_h3p_'+lab+'.pdf',title=titleintegr,step=stepint,polo=polo,cbarlabel=cbarlabint,cbarmult=-4,cbarform='%.1f',minu=0.0,maxu=maxint,show=False,aur_model=aurm,edges=edges[coto])
    jirfu.stereoplot(lon[coto],lat[coto],pixs.solar_time[coto]/2.4,nomefi=cartout2+'solar_time'+lab+'.pdf',polo=polo,cbarlabel='Solar time (hr)',aur_model=aurm,edges=edges[coto],invert_cmap=True)


    #
    # # Alcuni scatter plots
    # yy = lambda x: 100/0.4*x
    # jirfu.scatter(col,temp,err_t,cond = coto,lat=pixs.pc_lat,lon=pixs.pc_lon,ylim=(700,1200),xlim=(0,1e13),pdf=pdf_scatt_1,polo=polo, nome=cartout2+'SP_col_temp_errt'+lab+'.pdf',xlabel=r'H$_3^+$ col (cm^{-2})',ylabel='Temp (K)',title='color = Temp. err (K) : FILT {}')
    # jirfu.scatter(err_c,err_t,temp,cond=coto,lat=pixs.pc_lat,lon=pixs.pc_lon,ylim=(0,200),xlim=(0,1.5e12),clim=(700,tmax),pdf=pdf_scatt_2,polo=polo,
    #               nome=cartout2+'SP_errc_errt_temp'+lab+'.eps',xlabel='Error H3+ col (cm-2)',ylabel='Error Temp (K)',title='color = Temp (K) : FILT {}')
    # jirfu.scatter(err_c/col,err_t,temp,cond = coto,lat=pixs.pc_lat,lon=pixs.pc_lon,ylim=(0,200),xlim=(0,1),clim=(700,tmax),pdf=pdf_scatt_3,polo=polo,
    #               nome=cartout2+'SP_errcRel_errt_temp'+lab+'.eps',xlabel='Rel. err. H3+ col',ylabel='Error Temp (K)',title='color = Temp (K) : FILT {}')
    # jirfu.scatter(err_c/col,err_t,col,cond = coto,lat=pixs.pc_lat,lon=pixs.pc_lon,ylim=(0,200),xlim=(0,1),clim=(0,8e12),pdf=pdf_scatt_3,polo=polo,
    #               nome=cartout2+'SP_errcRel_errt_temp'+lab+'.eps',xlabel='Rel. err. H3+ col',ylabel='Error Temp (K)',title='color = H3+ col. (cm-2) : FILT {}')
    #
    # if filt == 1:
    #     gi = err_t < yy(err_c/col)
    #     jirfu.scatter(err_c/col,err_t,temp,cond = coto & gi,lat=pixs.pc_lat,lon=pixs.pc_lon,ylim=(0,200),xlim=(0,1),clim=(700,tmax),polo=polo,pdf=pdf_scatt_3,
    #                   nome=cartout2+'SP_errcRel_errt_temp_giu'+lab+'.eps',xlabel='Rel. err. H3+ col (cm-2)',ylabel='Error temp (K)',title='color = Temp (K) : FILT {}')
    #     gi = err_t > yy(err_c/col)
    #     jirfu.scatter(err_c/col,err_t,temp,cond = coto & gi,lat=pixs.pc_lat,lon=pixs.pc_lon,ylim=(0,200),xlim=(0,1),clim=(700,tmax),polo=polo,pdf=pdf_scatt_3,
    #                   nome=cartout2+'SP_errcRel_errt_temp_su'+lab+'.eps',xlabel='Rel. err. H3+ col (cm-2)',ylabel='Error temp (K)',title='color = Temp (K) : FILT {}')

    # stampo i risultati
    fi = open(cartout2+'tutti'+lab+'.dat','w')
    i = 0
    fu='{:05d}'+'{:>28s}'+2*'{:4d}'+5*'{:7.2f}'+'{:10.2e}'+3*'{:11.3e}'+2*'{:9.2f}'+'{:7.2f}'+'{:9.2e}'+'{:11.3e}'+'{:7.2f}'+'{:11.3e}'+'\n'
    fun='{:>5s}'+'{:>28s}'+2*'{:>4s}'+5*'{:>7s}'+'{:>10s}'+3*'{:>11s}'+2*'{:>9s}'+'{:>7s}'+'{:>9s}'+'{:>11s}'+'{:>7s}'+'{:>11s}'+'\n'
    fi.write('Results for cart {:40s}\n'.format(cartres))
    fi.write('Legend: num -> id number, cubo -> obs. cube, i -> line, j -> sample, lon -> planetocentric longitude, '
             'lat -> planetocentric latitude, emi -> emission angle, sza -> sza, stime -> local solar time, '
             'dist -> slant distance [m], col -> retrieved h3p column [cm-2], err_c -> retrieval error, '
             'col_c -> column corrected for emission angle, temp -> retrieved h3p temp [K], err_t -> retr. error, '
             'chi -> chi squared, off -> offset [W/m^2/sr/nm], ch4 -> retrieved ch4 column [cm-2], wls -> wl shift (nm), integr -> integrated intensity btw 3.35 and 3.75 nm [W/m^2/sr] (corrected for emiss_angle)\n')
    fi.write('{:1s}\n'.format(' '))
    fi.write(fun.format('num','cubo','i','j','lon','lat','emi','sza','stime','dist','col','err_c','col_c','temp',
                        'err_t','chi','off','ch4','wls','integr'))
    fi.write('{:1s}\n'.format('#'))
    for pi, ind in zip(pixs[coto], indx[coto]):
        fi.write(fu.format(ind,pi['cubo'],pi['line'],pi['sample'],pi['pc_lon'],pi['pc_lat'],pi['emiss_angle'],
                           pi['incid_angle'],pi['solar_time'],pi['slant_dist'],
                           col[coto][i],err_c[coto][i],col_c[coto][i],temp[coto][i],err_t[coto][i],chi[coto][i],off[coto][i],ch4[coto][i],shi[coto][i],integr[coto][i])) # num,lon,lat,emi,sza,loc,dist,))
        i+=1
    fi.close()

pdf_col.close()
pdf_temp.close()
pdf_scatt_1.close()
pdf_scatt_2.close()
pdf_scatt_3.close()
