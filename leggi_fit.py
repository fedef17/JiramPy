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
cart = '/home/fede/Scrivania/Jiram/DATA/JM0003_all/aur_N_nadir_70_CH4/'
cartres = cart
polo = 'N'
lch4 = True

aur_lon_0,aur_lat,aur_lon,aur_theta = jirfu.leggi_map_aur(polo)

######################################################################################à

pixs = pickle.load(open(cart+'all_pixs.pic','r'))

print(type(pixs.cubo))
cubs = np.unique(pixs.cubo)
print(cubs)

ii,col1, err_col = jirfu.read_res_jir_3(cartres+'CD-H3p.dat')
iit,temp1, err_temp = jirfu.read_res_jir_3(cartres+'VT-H3p.dat')
iic,chi1 = jirfu.read_res_jir_4(cartres+'chisq.dat')
iio,off1 = jirfu.read_res_jir_4(cartres+'offset_ok.dat')
if lch4: ii4,ch41, err_ch4 = jirfu.read_res_jir_3(cartres+'CD-CH4.dat')

col = np.zeros(len(pixs))
err_c = np.zeros(len(pixs))
temp = np.zeros(len(pixs))
err_t = np.zeros(len(pixs))
chi = np.zeros(len(pixs))
off = np.zeros(len(pixs))
if lch4: ch4 = np.zeros(len(pixs))

col[ii] = col1
err_c[ii] = err_col
temp[iit] = temp1
err_t[iit] = err_temp
chi[iic] = chi1
off[iio] = off1
if lch4: ch4[ii4] = ch41


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

# cond3 = (ch4 > 1e16)
# spech4 = np.mean(pixs[cond3].spe)
# # pl.plot(pixs.wl[0],spech4)
# # pl.show()
# zu = pixs[cond3]
# ofzu = off[cond3]
# zu = zu.view(np.recarray)
# for wl1,spe1,ofz in zip(zu.wl,zu.spe,ofzu):
#     print(ofz)
#     mask = jirfu.findspi(wl1,spe1)
#     pl.plot(wl1,spe1*mask-1e3*ofz)
# pl.show()
# sys.exit()


#######################

# nn = (pixs.incid_angle > 80) & (pixs.ind_h3p > 0.02)
# night = np.mean(pixs[nn].spe)
# pickle.dump(night,open('night.pic','w'))
#
# dd = (pixs.incid_angle < 70) & (pixs.ind_h3p > 0.02)
# day = np.mean(pixs[dd].spe)
# print(len(day))
# # cond3 = ((chi < 3) | (col >= 0.0) | (off > 0.0))
# # i=0
# # zu = pixs[cond3]
# # zu = zu.view(np.recarray)
# # for wl1,spe1 in zip(zu.wl[0:10000:100],zu.spe[0:10000:100]):
# #     mask = jirfu.findspi(wl1,spe1)
# #     pl.plot(wl1,spe1*mask,label=i)
# #     i+=1
# pl.plot(pixs.wl[0],day,label='day')
# pl.plot(pixs.wl[0],night,label='night')
# pl.grid()
# pl.xlabel('Wavelength (nm)')
# pl.ylabel('Radiance (W/m2/um/sr)')
# #pl.legend()
# pl.show()
# sys.exit()
# ###

### PROVA DIMENSIONE PIXEL
# ii = 1500
# vec1 = np.array([pixs.pc_lon_1[ii],pixs.pc_lon_2[ii],pixs.pc_lon_3[ii],pixs.pc_lon_4[ii]])
# vec2 = np.array([pixs.pc_lat_1[ii],pixs.pc_lat_2[ii],pixs.pc_lat_3[ii],pixs.pc_lat_4[ii]])
# print(vec1,vec2)
# #jirfu.stereopos(vec1,vec2,cartres+'prova_pixel_1.eps',title='Prova pixel',polo=polo)
# print(pixs.emiss_angle[ii])
# ii = 6500
# vec1 = np.array([pixs.pc_lon_1[ii],pixs.pc_lon_2[ii],pixs.pc_lon_3[ii],pixs.pc_lon_4[ii]])
# vec2 = np.array([pixs.pc_lat_1[ii],pixs.pc_lat_2[ii],pixs.pc_lat_3[ii],pixs.pc_lat_4[ii]])
# print(vec1,vec2)
# print(pixs.emiss_angle[ii])
# #jirfu.stereopos(vec1,vec2,cartres+'prova_pixel_2.eps',title='Prova pixel',polo=polo)
# sys.exit()

col_c = col*np.cos(np.pi*pixs.emiss_angle/180.0)
err_cc = err_c*np.cos(np.pi*pixs.emiss_angle/180.0)
print(pixs.start_time[0])
# pl.plot(pixs.alt_surf)
# pl.show()
#sbm.plotcorr(pixs.emiss_angle,col,cartres+'corr_emiss_col.eps',xlabel='Emission angle (deg)',ylabel='Column density (cm-2)')
#sbm.plotcorr(pixs.emiss_angle,col_c,cartres+'corr_emiss_col_c.eps',xlabel='Emission angle (deg)',ylabel='Column density (cm-2)')

# System III longitude increases towards W.
# note: the planetographic longitude goes towards W. So 30°=30°W=330°. So lon = 360-lon (il plot prende la lon che cresce verso Est)
# note: the planetocentric longitude goeas towards E. So ok for the plot.

# lon = 360-pixs.pg_lon
# lat = pixs.pg_lat
# lab = '_PG'
#
#
# print(ch4)
# jirfu.stereoplot(lon,lat,pixs.ind_h3p*np.cos(np.pi*pixs.emiss_angle/180.0),cartres+'ind_h3p'+lab+'.eps',min=0.01,max=0.025,title='H3+ integrated intensity [W/m2/sr]',polo=polo)
# jirfu.stereoplot(lon,lat,col_c,cartres+'h3p_col_'+lab+'.eps',title='H3+ column [mol/cm2] (corrected)',polo=polo,min=0.2e12,max=2.9e12)
# jirfu.stereoplot(lon,lat,col,cartres+'h3p_col_nocorr'+lab+'.eps',title='H3+ column [mol/cm2] (original)',polo=polo,min=0.2e12,max=8e12)
# jirfu.stereoplot(lon,lat,temp,cartres+'h3p_temp'+lab+'.eps',title='H3+ temperature [K]',polo=polo,min=700.,max=1100)
# if lch4: jirfu.stereoplot(lon,lat,ch4,cartres+'ch4_col'+lab+'.eps',title='Effective CH4 column [mol/cm2]',polo=polo)
# jirfu.stereoplot(lon,lat,pixs.emiss_angle,cartres+'emiss_angle'+lab+'.eps',polo=polo,title='Emission angle [deg]')
# jirfu.stereoplot(lon,lat,pixs.solar_time,cartres+'solar_time'+lab+'.eps',polo=polo,title='Solar time [hr]')
# jirfu.stereoplot(lon,lat,off,cartres+'offset'+lab+'.eps',polo=polo,title='Offset [W/m2/nm/sr]')
# jirfu.stereoplot(lon,lat,pixs.incid_angle,cartres+'incid_angle'+lab+'.eps',polo=polo,title='Incidence angle [deg]')
# jirfu.stereoplot(lon,lat,pixs.phase_angle,cartres+'phase_angle'+lab+'.eps',polo=polo,title='Phase angle [deg]')
# jirfu.stereoplot(lon,lat,np.arange(len(pixs)),cartres+'index'+lab+'.eps',polo=polo,title='Index')

# sys.exit()

lon = pixs.pc_lon
lat = pixs.pc_lat
lab = '_PC'

if polo == 'N':
    mcol = 3e12
else:
    mcol = 4e12

cartres = cartres + 'MAP_' + polo +'/'

#print(pixs[100])
#sys.exit()
cond4 = (pixs.solar_time < 12)
print(len(pixs[cond4]))
# nome = cartres+'MAP_h3p_col_am.eps'
# jirfu.stereomap(lon[cond4],lat[cond4],col_c[cond4],nome,title='H3+ column [mol/cm2] -> local time < 12',polo=polo,min=0.2e12,max=4.0e12,show=True)
# cond4 = (pixs.solar_time > 12)
# nome = cartres+'MAP_h3p_col_pm.eps'
# jirfu.stereomap(lon[cond4],lat[cond4],col_c[cond4],nome,title='H3+ column [mol/cm2] -> local time > 12',polo=polo,min=0.2e12,max=4.0e12,show=True)
cond4 = (pixs.alt_surf < 4e8)
nome = cartres+'MAP_h3p_col_close.eps'
jirfu.stereomap(lon[cond4],lat[cond4],col_c[cond4],nome,title='H3+ column [mol/cm2] -> dist < 4e5 km',polo=polo,min=0.2e12,max=mcol,show=False)
cond4 = (pixs.alt_surf > 4e8)
nome = cartres+'MAP_h3p_col_far.eps'
jirfu.stereomap(lon[cond4],lat[cond4],col_c[cond4],nome,title='H3+ column [mol/cm2] -> dist > 4e5 km',polo=polo,min=0.2e12,max=mcol,show=False)


# Mappe rapporto tra 3.2 e 3.54:
# zu = pixs
# ratio = np.zeros(len(zu))
# err = np.zeros(len(zu))
# for ww,spp,of,i in zip(zu.wl,zu.spe,off,range(len(zu))):
#     ratio[i], err[i] = jirfu.ratio32_35(ww,spp,of)
# jirfu.stereomap(lon,lat,ratio,cartres+'MAP_ratio_32_35.eps',title='H3+ ratio 3.2um/3.54um',
#                 polo=polo,form=True,addpoints=False,show=True,min=0.2,max=1.0,ncont=10)
# jirfu.stereomap(lon,lat,err/ratio,cartres+'ERR_ratio_32_35.eps',title='Relative error on H3+ ratio 3.2m/3.54um',
#                 polo=polo,form=True,addpoints=False,show=True,min=0.0,max=1.0,ncont=10,divide=True)
# print(len(ratio[~np.isnan(ratio)]))
# #jirfu.stereoplot(lon,lat,ratio,cartres+'Plot_ratio_32_35.eps',title='H3+ ratio 3.2um/3.54um',
# #                polo=polo,show=True,min=0.2,max=1.0)
# for ww,spp,of,i in zip(zu.wl,zu.spe,off,range(len(zu))):
#     ratio[i], err[i] = jirfu.ratio35_37(ww,spp,of)
# jirfu.stereomap(lon,lat,ratio,cartres+'MAP_ratio_35_37.eps',title='H3+ ratio 3.54um/3.67um',
#                 polo=polo,form=True,addpoints=False,show=True,min=1.0,max=2,ncont=10)
# jirfu.stereomap(lon,lat,err/ratio,cartres+'ERR_ratio_35_37.eps',title='Relative error on H3+ ratio 3.54um/3.67um',
#                 polo=polo,form=True,addpoints=False,show=True,min=0.0,max=1,ncont=10,divide=True)
# #jirfu.stereoplot(lon,lat,ratio,cartres+'Plot_ratio_35_37.eps',title='H3+ ratio 3.54um/3.67um',
# #                polo=polo,show=True,min=0.0,max=1.0)
# print(len(ratio[~np.isnan(ratio)]))
# for ww,spp,of,i in zip(zu.wl,zu.spe,off,range(len(zu))):
#     ratio[i], err[i] = jirfu.ratio32_37(ww,spp,of)
# jirfu.stereomap(lon,lat,ratio,cartres+'MAP_ratio_32_37.eps',title='H3+ ratio 3.2um/3.67um',
#                 polo=polo,form=True,addpoints=False,show=True,min=0.5,max=1.5,ncont=10)
# jirfu.stereomap(lon,lat,err/ratio,cartres+'ERR_ratio_32_37.eps',title='Relative error on H3+ ratio 3.2um/3.67um',
#                 polo=polo,form=True,addpoints=False,show=True,min=0.0,max=1.0,ncont=10,divide=True)
# #jirfu.stereoplot(lon,lat,ratio,cartres+'Plot_ratio_32_37.eps',title='H3+ ratio 3.2um/3.67um',
# #                polo=polo,show=True,min=0.5,max=1.5)
# print(len(ratio[~np.isnan(ratio)]))
# #sys.exit()
# # jirfu.stereomap(lon,lat,ratio,cartres+'MAP_ratio_32_35_nopoints.eps',title='H3+ ratio 3.2um/3.54um',
# #                 polo=polo,form=True,show=True)#,min=5e-1,max=1.0,ncont=10)
# # cond = pixs.alt_surf < 4e8
# # jirfu.stereomap(lon[cond],lat[cond],ratio[cond],cartres+'MAP_ratio_32_35_close.eps',title='H3+ ratio 3.2um/3.54um -> dist < 4e5 km',
# #                 polo=polo,form=True,addpoints=True,show=True,min=5e-1,max=1.0,ncont=10)
# # cond = pixs.alt_surf > 4e8
# # jirfu.stereomap(lon[cond],lat[cond],ratio[cond],cartres+'MAP_ratio_32_35_far.eps',title='H3+ ratio 3.2um/3.54um -> dist > 4e5 km',
# #                 polo=polo,form=True,show=True,min=5e-1,max=1.0,ncont=10)
#
# cond1 = ratio > 1.5
# spehot = np.mean(pixs[cond1].spe)
# cond1 = ratio < 1.5
# specold = np.mean(pixs[cond1].spe)
# pl.plot(pixs.wl[0],spehot)
# pl.plot(pixs.wl[0],specold)
# pl.show()
# #sys.exit()
# # sys.exit()

jirfu.stereomap(lon,lat,col_c,cartres+'MAP_h3p_col.eps',title='H3+ column [mol/cm2]',polo=polo,min=0.2e12,max=mcol, show=False)
jirfu.stereomap(lon,lat,err_cc/col_c,cartres+'ERR_h3p_col.eps',title='Relative error H3+ column [mol/cm2]',polo=polo,min=0.,max=1., show=False)
condt = (col > 0e12)
#condt = np.ones(len(col),dtype=bool)
condt = (col > 0e12)
lab = 'nofilt'
jirfu.stereomap(lon[condt],lat[condt],temp[condt],cartres+'MAP_h3p_temp'+lab+'.eps',title='H3+ temperature [K]', polo=polo,min=700.,max=1100, show=False)
jirfu.stereomap(lon[condt],lat[condt],err_t[condt],cartres+'ERR_h3p_temp'+lab+'.eps',title='Error H3+ temperature [K]', polo=polo, min=0.,max=200., show=False)
condt = (col > 1e12)
lab = '_1e12'
jirfu.stereomap(lon[condt],lat[condt],temp[condt],cartres+'MAP_h3p_temp'+lab+'.eps',title='H3+ temperature [K]', polo=polo,min=700.,max=1100, show=False)
jirfu.stereomap(lon[condt],lat[condt],err_t[condt],cartres+'ERR_h3p_temp'+lab+'.eps',title='Error H3+ temperature [K]', polo=polo, min=0.,max=200., show=False)
condt = (col > 2e12)
lab = '_2e12'
jirfu.stereomap(lon[condt],lat[condt],temp[condt],cartres+'MAP_h3p_temp'+lab+'.eps',title='H3+ temperature [K]', polo=polo,min=700.,max=1100, show=False)
jirfu.stereomap(lon[condt],lat[condt],err_t[condt],cartres+'ERR_h3p_temp'+lab+'.eps',title='Error H3+ temperature [K]', polo=polo, min=0.,max=200., show=False)

#sys.exit()
if lch4: jirfu.stereomap(lon,lat,ch4,cartres+'MAP_ch4_col.eps',title='Effective CH4 column [mol/cm2]',polo=polo, show=False)
print(np.min(pixs.ind_h3p))
#cond3 = ~np.isnan(off)
#zu = pixs[cond3]
zu = pixs
integr = np.zeros(len(zu))
for ww,spp,of,i in zip(zu.wl,zu.spe,off,range(len(zu))):
    mask = jirfu.findspi(ww,spp)
    integr[i] = jirfu.integr_h3p(ww,spp*mask,of)

integr=integr*np.cos(np.pi*zu.emiss_angle/180.0)
jirfu.stereomap(lon,lat,integr,cartres+'MAP_integr_h3p_conslit.eps',title='H3+ integrated intensity [W/m2/sr]',
                polo=polo,form=True,addpoints=True,show=False)
jirfu.stereomap(lon,lat,integr,cartres+'MAP_integr_h3p.eps',title='H3+ integrated intensity [W/m2/sr]',
                polo=polo,form=True,addpoints=False,show=False)

# cond3 = ~np.isnan(off)
# jirfu.stereomap(lon[cond3],lat[cond3],integr[cond3],cartres+'MAP_integr_h3p_solook.eps',title='H3+ integrated intensity [W/m2/sr]',
#                 polo=polo,form=True,addpoints=False,show=True)
#jirfu.stereomap(lon[cond3],lat[cond3],integr*np.cos(np.pi*zu.emiss_angle/180.0),cartres+'MAP_integr_h3p_nocorr.eps',title='H3+ integrated intensity [W/m2/sr]',polo=polo,form=True)

# sys.exit()

lab=''

#jirfu.stereoplot(lon,lat,pixs.ind_h3p*np.cos(np.pi*pixs.emiss_angle/180.0),cartres+'ind_h3p'+lab+'.eps',min=0.01,max=0.025,title='H3+ integrated intensity [W/m2/sr]',polo=polo)
jirfu.stereoplot(lon,lat,col_c,cartres+'h3p_col'+lab+'.eps',title='H3+ column [mol/cm2] (corrected)',polo=polo,min=0.2e12,max=mcol)
#jirfu.stereoplot(lon,lat,col,cartres+'h3p_col_nocorr'+lab+'.eps',title='H3+ column [mol/cm2] (original)',polo=polo,min=0.2e12,max=8e12)
jirfu.stereoplot(lon,lat,temp,cartres+'h3p_temp'+lab+'.eps',title='H3+ temperature [K]',polo=polo,min=700.,max=1100)
if lch4: jirfu.stereoplot(lon,lat,ch4,cartres+'ch4_col'+lab+'.eps',title='Effective CH4 column [mol/cm2]',polo=polo)
jirfu.stereoplot(lon,lat,pixs.emiss_angle,cartres+'emiss_angle'+lab+'.eps',polo=polo,title='Emission angle [deg]')
jirfu.stereoplot(lon,lat,pixs.solar_time,cartres+'solar_time'+lab+'.eps',polo=polo,title='Solar time [hr]')
jirfu.stereoplot(lon,lat,off,cartres+'offset'+lab+'.eps',polo=polo,title='Offset [W/m2/nm/sr]')
jirfu.stereoplot(lon,lat,pixs.incid_angle,cartres+'incid_angle'+lab+'.eps',polo=polo,title='Incidence angle [deg]')
jirfu.stereoplot(lon,lat,np.arange(len(pixs)),cartres+'index'+lab+'.eps',polo=polo,title='Index')

sys.exit()

emin = np.min(pixs.emiss_angle)
emax = np.max(pixs.emiss_angle)
omin = np.min(off)
omax = np.max(off)
for cub in cubs:
    cond = (pixs.cubo == cub)
    gimp = pixs[cond]
    if(len(gimp) < 20): continue

    lab = '_PC'
    lon = gimp.pc_lon
    lat = gimp.pc_lat

    nome = cartres+cub+'_h3p_col'+lab+'.eps'
    jirfu.stereoplot(lon,lat,col_c[cond],nome,title='H3+ column [mol/cm2] (corr) ->'+cub,polo=polo,min=0.2e12,max=2.9e12)
    nome = cartres+cub+'_h3p_temp'+lab+'.eps'
    jirfu.stereoplot(lon,lat,temp[cond],nome,title='H3+ temperature [K] ->'+cub,polo=polo,min=700.,max=1100)
    nome = cartres+cub+'_emiss_angle'+lab+'.eps'
    jirfu.stereoplot(lon,lat,gimp.emiss_angle,nome,polo=polo,title='Emiss. ang [deg] ->'+cub,min=emin,max=emax)
    nome = cartres+cub+'_offset'+lab+'.eps'
    jirfu.stereoplot(lon,lat,gimp.emiss_angle,nome,polo=polo,title='Emiss. ang [deg] ->'+cub,min=omin,max=omax)

# sys.exit()
#### CONTROLLO BIAS TEMPERATURA:
#
# cond1 = (temp > 1200)
# cond2 = (temp < 1050)
# uno = np.mean(pixs[cond1].spe)
# due = np.mean(pixs[cond2].spe)
# pl.plot(pixs.wl[0],uno, label='T > 1200 K')
# pl.plot(pixs.wl[0],due, label='T < 1050 K')
# pl.plot(pixs.wl[0],2*uno-due-0.001, label='2*uno-due')
# pl.grid()
# pl.legend(fontsize='small')
# pl.show()

# fig = pl.figure(figsize=(8, 6), dpi=150)
# pl.title('H3+ column [cm-2]')
# pl.scatter(pixs.pc_lon,pixs.pc_lat,c=col_c)
# pl.show()
# pl.colorbar(orientation="horizontal")
# fig.savefig(cartres+'col_h3p.eps', format='eps', dpi=150)
# pl.close()