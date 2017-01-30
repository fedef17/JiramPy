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

# from mayavi import mlab
# from mpl_toolkits.basemap import Basemap

# import spiceypy.spiceypy as spice
# spice.furnsh('/home/fede/Scrivania/Jiram/DATA/KERNELS_JIRAM/Kernels_jm0003/jm0003.mk')

# cart = '/home/fede/Scrivania/Jiram/DATA/TEST_marisa/Spe_prova_N_bis/'
# cartres = cart + 'Res_CH4/'
# cart = '/home/fede/Scrivania/Jiram/DATA/JM0003_all/aur_S_nadir_70/'
cart = '/home/fede/Scrivania/Jiram/DATA/JM0003_all/aur_N_nadir_70_CH4/'
cartres = cart
cartout = cart + 'MAP_orto/'
cart80 = '/home/fede/Scrivania/Jiram/DATA/JM0003_all/aur_N_nadir_80_CH4/'
polo = 'N'
lch4 = True
l500 = True  # così fa tutto con le geometrie a 500 km
lN80 = False  # True # così aggiungo le nuove misure fino a 80 gradi di emission angle
emissMAX = 80.0
if polo == 'S': lN80 = False
if l500 == False: lN80 = False

aur_lon_0, aur_lat, aur_lon, aur_theta = jirfu.leggi_map_aur(polo)

######################################################################################à

pixs = pickle.load(open(cart + 'all_pixs.pic', 'r'))

if l500:
    # si va a prendere le geometrie per 500 km
    p500 = pickle.load(open(cart + 'all_pixs_500.pic', 'r'))
    pixs500 = p500.view(np.recarray)
    print(len(pixs), len(pixs500))
    print(pixs[1139])
    print(pixs500[1139])

# sys.exit()

ii, col1, err_col = jirfu.read_res_jir_3(cartres + 'CD-H3p.dat')
iit, temp1, err_temp = jirfu.read_res_jir_3(cartres + 'VT-H3p.dat')
iic, chi1 = jirfu.read_res_jir_4(cartres + 'chisq.dat')
iio, off1 = jirfu.read_res_jir_4(cartres + 'offset_ok.dat')
iis, shi1 = jirfu.read_res_jir_4(cartres + 'shift_ok.dat')
if lch4: ii4, ch41, err_ch4 = jirfu.read_res_jir_3(cartres + 'CD-CH4.dat')

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
    p80 = pickle.load(open(cart80 + 'all_pixs.pic', 'r'))
    ii, col1, err_col = jirfu.read_res_jir_3(cart80 + 'CD-H3p.dat')
    iit, temp1, err_temp = jirfu.read_res_jir_3(cart80 + 'VT-H3p.dat')
    iic, chi1 = jirfu.read_res_jir_4(cart80 + 'chisq.dat')
    iio, off1 = jirfu.read_res_jir_4(cart80 + 'offset_ok.dat')
    iis, shi1 = jirfu.read_res_jir_4(cart80 + 'shift_ok.dat')
    if lch4: ii4, ch41, err_ch4 = jirfu.read_res_jir_3(cart80 + 'CD-CH4.dat')

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

    print(type(pixs), len(pixs))
    p80 = p80.view(np.recarray)
    co80 = (p80.emiss_angle < emissMAX)
    p80 = p80[co80]
    pixs = np.append(pixs, p80)
    pixs500 = np.append(pixs500, p80)
    print(type(pixs), len(pixs))
    pixs = pixs.view(np.recarray)
    pixs500 = pixs500.view(np.recarray)
    print(type(pixs), len(pixs))
    col = np.append(col, col80[co80])
    err_c = np.append(err_c, err_c80[co80])
    temp = np.append(temp, temp80[co80])
    err_t = np.append(err_t, err_t80[co80])
    chi = np.append(chi, chi80[co80])
    off = np.append(off, off80[co80])
    shi = np.append(shi, shi80[co80])
    if lch4: ch4 = np.append(ch4, ch480[co80])

print(type(pixs.cubo))
cubs = np.unique(pixs500.cubo)
i = 0
for cu in cubs:
    if not np.all(np.isnan(pixs500[pixs500.cubo == cu].pc_lon)):
        em = np.mean(pixs500[pixs500.cubo == cu].emiss_angle)
        so = np.mean(pixs500[pixs500.cubo == cu].solar_time)
        n = len(pixs500[pixs500.cubo == cu])
        print('{}  {}  {:5.1f}  {:5.1f} {}'.format(i, cu, em, so, n))
        # print( pixs500[pixs500.cubo == cu].pc_lon)
        i += 1
# sys.exit()

col_c = col * np.cos(np.pi * pixs.emiss_angle / 180.0)
err_cc = err_c * np.cos(np.pi * pixs.emiss_angle / 180.0)

fi = open(cart + 'tutti.dat', 'w')
i = 0
fu = '{:05d}' + '{:>28s}' + 2 * '{:4d}' + 5 * '{:7.2f}' + '{:10.2e}' + 3 * '{:11.3e}' + 2 * '{:9.2f}' + '{:7.2f}' + '{:9.2e}' + '{:11.3e}' + '{:7.2f}' + '\n'
fun = '{:>5s}' + '{:>28s}' + 2 * '{:>4s}' + 5 * '{:>7s}' + '{:>10s}' + 3 * '{:>11s}' + 2 * '{:>9s}' + '{:>7s}' + '{:>9s}' + '{:>11s}' + '{:>7s}' + '\n'
fi.write('Results for cart {:40s}\n'.format(cart))
fi.write('Legend: num -> id number, cubo -> obs. cube, i -> line, j -> sample, lon -> planetocentric longitude, '
         'lat -> planetocentric latitude, emi -> emission angle, sza -> sza, stime -> local solar time, '
         'dist -> slant distance [m], col -> retrieved h3p column [cm-2], err_c -> retrieval error, '
         'col_c -> column corrected for emission angle, temp -> retrieved h3p temp [K], err_t -> retr. error, '
         'chi -> chi squared, off -> offset [W/m^2/sr/nm], ch4 -> retrieved ch4 column [cm-2], wls -> wl shift (nm) \n')
fi.write('{:1s}\n'.format(' '))
fi.write(
    fun.format('num', 'cubo', 'i', 'j', 'lon', 'lat', 'emi', 'sza', 'stime', 'dist', 'col', 'err_c', 'col_c', 'temp',
               'err_t', 'chi', 'off', 'ch4', 'wls'))
fi.write('{:1s}\n'.format('#'))
for pi in pixs500:
    # if (chi[i] > 4 or col[i] < 0.0 or off[i] < 0.0):
    #     i+=1
    #     continue
    fi.write(fu.format(i, pi['cubo'], pi['line'], pi['sample'], pi['pc_lon'], pi['pc_lat'], pi['emiss_angle'],
                       pi['incid_angle'], pi['solar_time'], pi['slant_dist'],
                       col[i], err_c[i], col_c[i], temp[i], err_t[i], chi[i], off[i], ch4[i],
                       shi[i]))  # num,lon,lat,emi,sza,loc,dist,))
    i += 1
fi.close()
sys.exit()

cond2 = ((chi > 4) | (col <= 0.0) | (off < 0.0))
sbm.plotcorr(pixs[~cond2].sample, shi[~cond2], cart + 'WL_shift.eps', xlabel='sample', ylabel='wl shift (nm)',
             xlim=[0, 260], ylim=[-3, +3])
condshi = (chi < 4) & (col > 2e12)
sbm.plotcorr(pixs[condshi].sample, shi[condshi], cart + 'WL_shift_selected.eps', xlabel='sample',
             ylabel='wl shift (nm)', xlim=[0, 260], ylim=[-3, +3])

col[cond2] = float('NaN')
print(len(col[cond2]))
temp[cond2] = float('NaN')
chi[cond2] = float('NaN')
off[cond2] = float('NaN')
err_t[cond2] = float('NaN')
err_c[cond2] = float('NaN')

if lch4: ch4[cond2] = float('NaN')
if lch4: ch4[(ch4 < 0.0)] = 0.0

# lab=''
# mcol=4e12
# npix=100
# jirfu.stereomap2(pixs.pc_lon,pixs.pc_lat,col_c,cartout+'MAP_h3p_col'+lab+'_image.eps',title='H3+ column [mol/cm2]',polo=polo,min=0.2e12,max=mcol, show=True, xres=npix, image=False, aur_model='stat')
# sys.exit()

# giu = (err_t < 100)
# pl.scatter(col_c[giu],temp[giu],c=err_t[giu],s=2,edgecolors='none',vmin=0,vmax=200)
# pl.xlabel('H3p col')
# pl.ylabel('H3p temp')
# #pl.scatter(temp,err_t,s=1)
# pl.colorbar()
# pl.grid()

# his,bin,pat = pl.hist(err_t,bins=20,range=(0,300))
# pl.show()


# pickle.dump([pixs,col,col_c,err_c,temp,err_t,off,chi,ch4],open('res_'+polo+'.pic','w'))

# sys.exit()


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
# #jirfu.stereopos(vec1,vec2,cartout+'prova_pixel_1.eps',title='Prova pixel',polo=polo)
# print(pixs.emiss_angle[ii])
# ii = 6500
# vec1 = np.array([pixs.pc_lon_1[ii],pixs.pc_lon_2[ii],pixs.pc_lon_3[ii],pixs.pc_lon_4[ii]])
# vec2 = np.array([pixs.pc_lat_1[ii],pixs.pc_lat_2[ii],pixs.pc_lat_3[ii],pixs.pc_lat_4[ii]])
# print(vec1,vec2)
# print(pixs.emiss_angle[ii])
# #jirfu.stereopos(vec1,vec2,cartout+'prova_pixel_2.eps',title='Prova pixel',polo=polo)
# sys.exit()

col_c = col * np.cos(np.pi * pixs.emiss_angle / 180.0)
err_cc = err_c * np.cos(np.pi * pixs.emiss_angle / 180.0)
if l500:
    col_c = col * np.cos(np.pi * pixs500.emiss_angle / 180.0)
    err_cc = err_c * np.cos(np.pi * pixs500.emiss_angle / 180.0)
print(pixs.start_time[0])
# pl.plot(pixs.alt_surf)
# pl.show()
# sbm.plotcorr(pixs.emiss_angle,col,cartout+'corr_emiss_col.eps',xlabel='Emission angle (deg)',ylabel='Column density (cm-2)')
# sbm.plotcorr(pixs.emiss_angle,col_c,cartout+'corr_emiss_col_c.eps',xlabel='Emission angle (deg)',ylabel='Column density (cm-2)')

# System III longitude increases towards W.
# note: the planetographic longitude goes towards W. So 30°=30°W=330°. So lon = 360-lon (il plot prende la lon che cresce verso Est)
# note: the planetocentric longitude goeas towards E. So ok for the plot.

# lon = 360-pixs.pg_lon
# lat = pixs.pg_lat
# lab = '_PG'
#
#
# print(ch4)
# jirfu.stereoplot(lon,lat,pixs.ind_h3p*np.cos(np.pi*pixs.emiss_angle/180.0),cartout+'ind_h3p'+lab+'.eps',min=0.01,max=0.025,title='H3+ integrated intensity [W/m2/sr]',polo=polo)
# jirfu.stereoplot(lon,lat,col_c,cartout+'h3p_col_'+lab+'.eps',title='H3+ column [mol/cm2] (corrected)',polo=polo,min=0.2e12,max=2.9e12)
# jirfu.stereoplot(lon,lat,col,cartout+'h3p_col_nocorr'+lab+'.eps',title='H3+ column [mol/cm2] (original)',polo=polo,min=0.2e12,max=8e12)
# jirfu.stereoplot(lon,lat,temp,cartout+'h3p_temp'+lab+'.eps',title='H3+ temperature [K]',polo=polo,min=700.,max=1100)
# if lch4: jirfu.stereoplot(lon,lat,ch4,cartout+'ch4_col'+lab+'.eps',title='Effective CH4 column [mol/cm2]',polo=polo)
# jirfu.stereoplot(lon,lat,pixs.emiss_angle,cartout+'emiss_angle'+lab+'.eps',polo=polo,title='Emission angle [deg]')
# jirfu.stereoplot(lon,lat,pixs.solar_time,cartout+'solar_time'+lab+'.eps',polo=polo,title='Solar time [hr]')
# jirfu.stereoplot(lon,lat,off,cartout+'offset'+lab+'.eps',polo=polo,title='Offset [W/m2/nm/sr]')
# jirfu.stereoplot(lon,lat,pixs.incid_angle,cartout+'incid_angle'+lab+'.eps',polo=polo,title='Incidence angle [deg]')
# jirfu.stereoplot(lon,lat,pixs.phase_angle,cartout+'phase_angle'+lab+'.eps',polo=polo,title='Phase angle [deg]')
# jirfu.stereoplot(lon,lat,np.arange(len(pixs)),cartout+'index'+lab+'.eps',polo=polo,title='Index')

# sys.exit()

lon = pixs.pc_lon
lat = pixs.pc_lat
lab = ''
if l500:
    lon = pixs500.pc_lon
    lat = pixs500.pc_lat
    lab = '_500'

if lN80:
    lab = '_lN80'

if polo == 'N':
    mcol = 3e12
    npix = 40
else:
    mcol = 4e12
    npix = 50

# cartout = cartout + 'MAP_' + polo +'/'

# print(pixs[100])
# sys.exit()
cond4 = (pixs.solar_time < 12)
print(len(pixs[cond4]))
# nome = cartout+'MAP_h3p_col_am.eps'
# jirfu.stereomap(lon[cond4],lat[cond4],col_c[cond4],nome,title='H3+ column [mol/cm2] -> local time < 12',polo=polo,min=0.2e12,max=4.0e12,show=True)
# cond4 = (pixs.solar_time > 12)
# nome = cartout+'MAP_h3p_col_pm.eps'
# jirfu.stereomap(lon[cond4],lat[cond4],col_c[cond4],nome,title='H3+ column [mol/cm2] -> local time > 12',polo=polo,min=0.2e12,max=4.0e12,show=True)

# Mappe rapporto tra 3.2 e 3.54:
# zu = pixs
# ratio = np.zeros(len(zu))
# err = np.zeros(len(zu))
# for ww,spp,of,i in zip(zu.wl,zu.spe,off,range(len(zu))):
#     ratio[i], err[i] = jirfu.ratio32_35(ww,spp,of)
# jirfu.stereomap(lon,lat,ratio,cartout+'MAP_ratio_32_35.eps',title='H3+ ratio 3.2um/3.54um',
#                 polo=polo,form=True,addpoints=False,show=True,min=0.2,max=1.0,ncont=10)
# jirfu.stereomap(lon,lat,err/ratio,cartout+'ERR_ratio_32_35.eps',title='Relative error on H3+ ratio 3.2m/3.54um',
#                 polo=polo,form=True,addpoints=False,show=True,min=0.0,max=1.0,ncont=10,divide=True)
# print(len(ratio[~np.isnan(ratio)]))
# #jirfu.stereoplot(lon,lat,ratio,cartout+'Plot_ratio_32_35.eps',title='H3+ ratio 3.2um/3.54um',
# #                polo=polo,show=True,min=0.2,max=1.0)
# for ww,spp,of,i in zip(zu.wl,zu.spe,off,range(len(zu))):
#     ratio[i], err[i] = jirfu.ratio35_37(ww,spp,of)
# jirfu.stereomap(lon,lat,ratio,cartout+'MAP_ratio_35_37.eps',title='H3+ ratio 3.54um/3.67um',
#                 polo=polo,form=True,addpoints=False,show=True,min=1.0,max=2,ncont=10)
# jirfu.stereomap(lon,lat,err/ratio,cartout+'ERR_ratio_35_37.eps',title='Relative error on H3+ ratio 3.54um/3.67um',
#                 polo=polo,form=True,addpoints=False,show=True,min=0.0,max=1,ncont=10,divide=True)
# #jirfu.stereoplot(lon,lat,ratio,cartout+'Plot_ratio_35_37.eps',title='H3+ ratio 3.54um/3.67um',
# #                polo=polo,show=True,min=0.0,max=1.0)
# print(len(ratio[~np.isnan(ratio)]))
# for ww,spp,of,i in zip(zu.wl,zu.spe,off,range(len(zu))):
#     ratio[i], err[i] = jirfu.ratio32_37(ww,spp,of)
# jirfu.stereomap(lon,lat,ratio,cartout+'MAP_ratio_32_37.eps',title='H3+ ratio 3.2um/3.67um',
#                 polo=polo,form=True,addpoints=False,show=True,min=0.5,max=1.5,ncont=10)
# jirfu.stereomap(lon,lat,err/ratio,cartout+'ERR_ratio_32_37.eps',title='Relative error on H3+ ratio 3.2um/3.67um',
#                 polo=polo,form=True,addpoints=False,show=True,min=0.0,max=1.0,ncont=10,divide=True)
# #jirfu.stereoplot(lon,lat,ratio,cartout+'Plot_ratio_32_37.eps',title='H3+ ratio 3.2um/3.67um',
# #                polo=polo,show=True,min=0.5,max=1.5)
# print(len(ratio[~np.isnan(ratio)]))
# #sys.exit()
# # jirfu.stereomap(lon,lat,ratio,cartout+'MAP_ratio_32_35_nopoints.eps',title='H3+ ratio 3.2um/3.54um',
# #                 polo=polo,form=True,show=True)#,min=5e-1,max=1.0,ncont=10)
# # cond = pixs.alt_surf < 4e8
# # jirfu.stereomap(lon[cond],lat[cond],ratio[cond],cartout+'MAP_ratio_32_35_close.eps',title='H3+ ratio 3.2um/3.54um -> dist < 4e5 km',
# #                 polo=polo,form=True,addpoints=True,show=True,min=5e-1,max=1.0,ncont=10)
# # cond = pixs.alt_surf > 4e8
# # jirfu.stereomap(lon[cond],lat[cond],ratio[cond],cartout+'MAP_ratio_32_35_far.eps',title='H3+ ratio 3.2um/3.54um -> dist > 4e5 km',
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

# pl.scatter(col_c,temp,c=err_t,s=2,edgecolors='none',vmin=0,vmax=200)
# pl.colorbar()
# pl.grid()
# pl.show()
# pl.close()

# MAPPE:

jirfu.stereomap2(lon, lat, col_c, nomefi=cartout + 'MAP_h3p_col' + lab + '_image.eps', title='H3+ column [mol/cm2]',
                 polo=polo, min=0.2e12, max=mcol, show=True, xres=npix, image=False, aur_model='both')
sys.exit()
co = err_t < 60
jirfu.stereomap2(lon, lat, col_c, nomefi=cartout + 'MAP_h3p_col' + lab + '_errTm60.eps', title='H3+ column [mol/cm2]',
                 polo=polo, min=0.2e12, max=mcol, show=True, xres=npix, image=False, addpoints=True, condpo=co)
co = err_t > 80
jirfu.stereomap2(lon, lat, col_c, nomefi=cartout + 'MAP_h3p_col' + lab + '_errTp60.eps', title='H3+ column [mol/cm2]',
                 polo=polo, min=0.2e12, max=mcol, show=True, xres=npix, image=False, addpoints=True, condpo=co)
jirfu.stereomap2(lon, lat, err_cc / col_c, nomefi=cartout + 'ERR_h3p_col' + lab + '.eps',
                 title='Relative error H3+ column [mol/cm2]', polo=polo, min=0., max=1., show=False, xres=npix)

jirfu.stereomap2(lon, lat, col_c, nomefi=cartout + 'Npoints_ok' + lab + '.eps', title='Number of measurements',
                 polo=polo, show=True, xres=50, npoints=True)  # ,addpoints=True)

co = pixs500.emiss_angle < 60
jirfu.stereomap2(lon[co], lat[co], col_c[co], nomefi=cartout + 'MAP_h3p_col' + lab + '_emissm60.eps',
                 title='H3+ column [mol/cm2]', polo=polo, min=0.2e12, max=mcol, show=True, xres=npix)
jirfu.stereoplot(lon[co], lat[co], col_c[co], nomefi=cartout + 'PLOT_h3p_col' + lab + '_emissm60.eps',
                 title='H3+ column [mol/cm2]', polo=polo, min=0.2e12, max=mcol, show=True)
co = pixs500.emiss_angle > 60
jirfu.stereomap2(lon[co], lat[co], col_c[co], nomefi=cartout + 'MAP_h3p_col' + lab + '_emissp60.eps',
                 title='H3+ column [mol/cm2]', polo=polo, min=0.2e12, max=mcol, show=True, xres=npix)
jirfu.stereoplot(lon[co], lat[co], col_c[co], nomefi=cartout + 'PLOT_h3p_col' + lab + '_emissp60.eps',
                 title='H3+ column [mol/cm2]', polo=polo, min=0.2e12, max=mcol, show=True)

sys.exit()

cond4 = (pixs.alt_surf < 4e8)
nome = cartout + 'MAP_h3p_col_close' + lab + '.eps'
jirfu.stereomap2(lon[cond4], lat[cond4], col_c[cond4], nomefi=nome, title='H3+ column [mol/cm2] -> dist < 4e5 km',
                 polo=polo, min=0.2e12, max=mcol, show=False, xres=npix)
cond4 = (pixs.alt_surf > 4e8)
nome = cartout + 'MAP_h3p_col_far' + lab + '.eps'
jirfu.stereomap2(lon[cond4], lat[cond4], col_c[cond4], nomefi=nome, title='H3+ column [mol/cm2] -> dist > 4e5 km',
                 polo=polo, min=0.2e12, max=mcol, show=False, xres=npix)

condt = (col > 0e12)
# condt = np.ones(len(col),dtype=bool)
giu = err_t < 100
lab2 = '_err100K'
jirfu.stereomap2(lon[giu], lat[giu], temp[giu], nomefi=cartout + 'MAP_h3p_temp' + lab2 + lab + '.eps',
                 title='H3+ temperature [K]', polo=polo, min=700., max=1100, show=False, xres=npix)
jirfu.stereomap2(lon[giu], lat[giu], err_t[giu], nomefi=cartout + 'ERR_h3p_temp' + lab2 + lab + '.eps',
                 title='Error H3+ temperature [K]', polo=polo, min=0., max=200., show=False, xres=npix)
giu = err_t < 80
lab2 = '_err80K'
jirfu.stereomap2(lon[giu], lat[giu], temp[giu], nomefi=cartout + 'MAP_h3p_temp' + lab2 + lab + '.eps',
                 title='H3+ temperature [K]', polo=polo, min=700., max=1100, show=False, xres=npix)
jirfu.stereomap2(lon[giu], lat[giu], err_t[giu], nomefi=cartout + 'ERR_h3p_temp' + lab2 + lab + '.eps',
                 title='Error H3+ temperature [K]', polo=polo, min=0., max=200., show=False, xres=npix)

condt = (col > 0e12)
lab2 = '_nofilt'
jirfu.stereomap2(lon[condt], lat[condt], temp[condt], nomefi=cartout + 'MAP_h3p_temp' + lab2 + lab + '.eps',
                 title='H3+ temperature [K]', polo=polo, min=700., max=1100, show=False, xres=npix)
jirfu.stereomap2(lon[condt], lat[condt], err_t[condt], nomefi=cartout + 'ERR_h3p_temp' + lab2 + lab + '.eps',
                 title='Error H3+ temperature [K]', polo=polo, min=0., max=200., show=False, xres=npix)
# condt = (col > 1e12)
# lab2 = '_1e12'
# jirfu.stereomap2(lon[condt],lat[condt],temp[condt],nomefi=cartout+'MAP_h3p_temp'+lab2+lab+'.eps',title='H3+ temperature [K]', polo=polo,min=700.,max=1100, show=False,xres=npix)
# jirfu.stereomap2(lon[condt],lat[condt],err_t[condt],nomefi=cartout+'ERR_h3p_temp'+lab2+lab+'.eps',title='Error H3+ temperature [K]', polo=polo, min=0.,max=200., show=False,xres=npix)
# condt = (col > 2e12)
# lab2 = '_2e12'
# jirfu.stereomap2(lon[condt],lat[condt],temp[condt],nomefi=cartout+'MAP_h3p_temp'+lab2+lab+'.eps',title='H3+ temperature [K]', polo=polo,min=700.,max=1100, show=False,xres=npix)
# jirfu.stereomap2(lon[condt],lat[condt],err_t[condt],nomefi=cartout+'ERR_h3p_temp'+lab2+lab+'.eps',title='Error H3+ temperature [K]', polo=polo, min=0.,max=200., show=False,xres=npix)

# sys.exit()

if lch4: jirfu.stereomap2(lon, lat, ch4, nomefi=cartout + 'MAP_ch4_col' + lab + '.eps',
                          title='Effective CH4 column [mol/cm2]', polo=polo, show=False, xres=npix)
print(np.min(pixs.ind_h3p))
# cond3 = ~np.isnan(off)
# zu = pixs[cond3]
zu = pixs
integr = np.zeros(len(zu))
for ww, spp, of, i in zip(zu.wl, zu.spe, off, range(len(zu))):
    mask = jirfu.findspi(ww, spp)
    integr[i] = jirfu.integr_h3p(ww, spp * mask, of)

integr = integr * np.cos(np.pi * zu.emiss_angle / 180.0)
jirfu.stereomap2(lon, lat, integr, nomefi=cartout + 'MAP_integr_h3p_conslit' + lab + '.eps',
                 title='H3+ integrated intensity [W/m2/sr]',
                 polo=polo, form=True, addpoints=True, show=False, xres=npix)
jirfu.stereomap2(lon, lat, integr, nomefi=cartout + 'MAP_integr_h3p' + lab + '.eps',
                 title='H3+ integrated intensity [W/m2/sr]',
                 polo=polo, form=True, addpoints=False, show=False, xres=npix)

# cond3 = ~np.isnan(off)
# jirfu.stereomap2(lon[cond3],lat[cond3],integr[cond3],nomefi=cartout+'MAP_integr_h3p_solook.eps',title='H3+ integrated intensity [W/m2/sr]',
#                 polo=polo,form=True,addpoints=False,show=True)
# jirfu.stereomap2(lon[cond3],lat[cond3],integr*np.cos(np.pi*zu.emiss_angle/180.0),nomefi=cartout+'MAP_integr_h3p_nocorr.eps',title='H3+ integrated intensity [W/m2/sr]',polo=polo,form=True)

sys.exit()

# lab=''

# jirfu.stereoplot(lon,lat,pixs.ind_h3p*np.cos(np.pi*pixs.emiss_angle/180.0),nomefi=cartout+'ind_h3p'+lab+'.eps',min=0.01,max=0.025,title='H3+ integrated intensity [W/m2/sr]',polo=polo)
jirfu.stereoplot(lon, lat, col_c, nomefi=cartout + 'h3p_col' + lab + '.eps', title='H3+ column [mol/cm2] (corrected)',
                 polo=polo, min=0.2e12, max=mcol)
# jirfu.stereoplot(lon,lat,col,nomefi=cartout+'h3p_col_nocorr'+lab+'.eps',title='H3+ column [mol/cm2] (original)',polo=polo,min=0.2e12,max=8e12)
jirfu.stereoplot(lon, lat, temp, nomefi=cartout + 'h3p_temp' + lab + '.eps', title='H3+ temperature [K]', polo=polo,
                 min=700., max=1100)
if lch4: jirfu.stereoplot(lon, lat, ch4, nomefi=cartout + 'ch4_col' + lab + '.eps',
                          title='Effective CH4 column [mol/cm2]', polo=polo)
jirfu.stereoplot(lon, lat, pixs.emiss_angle, nomefi=cartout + 'emiss_angle' + lab + '.eps', polo=polo,
                 title='Emission angle [deg]')
jirfu.stereoplot(lon, lat, pixs.solar_time, nomefi=cartout + 'solar_time' + lab + '.eps', polo=polo,
                 title='Solar time [hr]')
jirfu.stereoplot(lon, lat, off, nomefi=cartout + 'offset' + lab + '.eps', polo=polo, title='Offset [W/m2/nm/sr]')
jirfu.stereoplot(lon, lat, pixs.incid_angle, nomefi=cartout + 'incid_angle' + lab + '.eps', polo=polo,
                 title='Incidence angle [deg]')
jirfu.stereoplot(lon, lat, pixs.phase_angle, nomefi=cartout + 'phase_angle' + lab + '.eps', polo=polo,
                 title='Phase angle [deg]')
jirfu.stereoplot(lon, lat, np.arange(len(pixs)), nomefi=cartout + 'index' + lab + '.eps', polo=polo, title='Index')

sys.exit()

emin = np.min(pixs.emiss_angle)
emax = np.max(pixs.emiss_angle)
omin = np.min(off)
omax = np.max(off)
for cub in cubs:
    cond = (pixs.cubo == cub)
    gimp = pixs[cond]
    if (len(gimp) < 20): continue

    lab = '_PC'
    lon = gimp.pc_lon
    lat = gimp.pc_lat

    nome = cartout + cub + '_h3p_col' + lab + '.eps'
    jirfu.stereoplot(lon, lat, col_c[cond], nomefi=nome, title='H3+ column [mol/cm2] (corr) ->' + cub, polo=polo,
                     min=0.2e12, max=2.9e12)
    nome = cartout + cub + '_h3p_temp' + lab + '.eps'
    jirfu.stereoplot(lon, lat, temp[cond], nomefi=nome, title='H3+ temperature [K] ->' + cub, polo=polo, min=700.,
                     max=1100)
    nome = cartout + cub + '_emiss_angle' + lab + '.eps'
    jirfu.stereoplot(lon, lat, gimp.emiss_angle, nomefi=nome, polo=polo, title='Emiss. ang [deg] ->' + cub, min=emin,
                     max=emax)
    nome = cartout + cub + '_offset' + lab + '.eps'
    jirfu.stereoplot(lon, lat, gimp.emiss_angle, nomefi=nome, polo=polo, title='Emiss. ang [deg] ->' + cub, min=omin,
                     max=omax)

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
    # fig.savefig(cartout+'col_h3p.eps', format='eps', dpi=150)
    # pl.close()
