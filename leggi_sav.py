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

cart = '/home/federico/Jiram/DATA/JM0003_SAV/'

nomeout = 'aur_JM0003_limb2000.pic'
thres = 0.005 # soglia per l'indice di H3+
mincu = 30 # minimo numero di pixel buoni nel cubo

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

with open(cart+'lista','r') as lista:
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
        for i in range(n_lin):
            for j in range(n_sam):
                corners = corners_
                pix = pix_
	        #if(geo[i,j,28] == 0):
	        #    print(geo[i,j,32])

                #cond = (geo[i,j,17] > 60.0 or geo[i,j,17] < -60.0) and (geo[i,j,28] == 1)# or geo[i,j,32] < 3000.0)
                #cond = geo[i,j,17] > 50.0 and (geo[i,j,28] == 1)# or geo[i,j,32] < 3000.0)
                cond = (geo[i,j,28] == 0 and geo[i,j,32]<2e7)
                geo[i,j,32]=1e-3*geo[i,j,32]
                #cond = geo[i,j,17] < -60.0 and (geo[i,j,28] == 1)# or geo[i,j,32] < 3000.0)
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

                    fu = jirfu.checkqual(wls,spe[i,j,:])

                    if(ind_h3p < thres): continue
                    if(fu == 0): continue

                    pixcu = np.append(pixcu, pix)
#
        pixcu = pixcu[1:]
        print(len(pixcu))
        if(len(pixcu) < mincu): continue
        pixtot = np.append(pixtot,pixcu)
        print(len(pixtot))

pixtot = pixtot[1:]
pickle.dump(pixtot,open(cart+nomeout,'w'))

#         print(np.shape(spe))
#         print(np.shape(geo))
#         #pl.imshow(spe[:,:,171], cmap='jet', interpolation="nearest")
#         # #pl.colorbar(orientation="horizontal")
#         #pl.show()
#         # pl.imshow(geo[:,:,25], cmap='jet', interpolation="nearest")
#         # pl.colorbar(orientation="horizontal")
#         # pl.show()
#         # sys.exit()

# pixtot = pickle.load(open(cart+'pix_nadir_N.pic','r'))
# pixs = pixtot.view(np.recarray)
#
# cond = (pixs.emiss_angle < 80) & (pixs.ind_h3p > 0.015)
# pl.scatter(pixs[cond].pc_lon,pixs[cond].pc_lat,c=pixs[cond].ind_h3p)
# pl.colorbar()
# pl.show()
# pl.close()
# #pl.plot(pixs[cond].wl[13], pixs[cond].spe[13])
# print(pixs[cond].ind_h3p[13])
# #pl.show()
# print('ok')
#
# print(len(pixs[cond]))
# n = len(pixs[cond])
#
# # for i in range(n):
# #     nome = cart2+'Spe_'+'{0:0=4d}'.format(i)+'.dat'
# #     print(nome)
# #     err = jirfu.findspi(wls,spe)
# #     jirfu.write_obs_JIR(pixs[cond].wl[i],pixs[cond].spe[i],nome,comment='emiss < 80 and ind_h3p > 0.015')
# #     # call(cart+'./fomi_mipas')
# #     # nomeout = cart+'output__mipas.dat'
# #     # alt_fomi, cr_fomi = sbm.leggioutfomi(nomeout)
#
# cartres = '/home/fede/Scrivania/Jiram/DATA/TEST_marisa/Res_prova/'
#
# ii,col, err_col = jirfu.read_res_jir_3(cartres+'CD-H3p.dat')
# ii,temp, err_temp = jirfu.read_res_jir_3(cartres+'VT-H3p.dat')
# ii,chi = jirfu.read_res_jir_4(cartres+'chisq.dat')
# ii,off = jirfu.read_res_jir_4(cartres+'offset_ok.dat')
#
# cond2 = ((chi > 4) | (col < 0.0) | (off < 0.0))
# print(len(col))
# print(type(col))
# print(len(col[cond2]))
# col[cond2] = float('NaN')
# temp[cond2] = float('NaN')
# chi[cond2] = float('NaN')
# off[cond2] = float('NaN')
#
# print(pixs[cond].emiss_angle)
#
# col_c = col*np.cos(np.pi*pixs[cond].emiss_angle/180.0)
#
# pl.scatter(pixs[cond].pc_lon,pixs[cond].pc_lat,c=col_c)
# pl.colorbar()
# pl.show()
#
# pl.scatter(pixs[cond].pc_lon,pixs[cond].pc_lat,c=temp)
# pl.colorbar()
# pl.show()
#
# pl.scatter(pixs[cond].pc_lon,pixs[cond].pc_lat,c=chi)
# pl.colorbar()
# pl.show()
#
# pl.scatter(pixs[cond].pc_lon,pixs[cond].pc_lat,c=off)
# pl.colorbar()
# pl.show()


# stepo = np.vectorize(jirfu.stereopolar)
# x,y = stepo(pixs[cond].pg_lon,pixs[cond].pg_lat)
# pl.scatter(x,y,c=col_c)
# pl.colorbar()
# pl.show()

# # Create a sphere
# r = 1.0
# pi = np.pi
# cos = np.cos
# sin = np.sin
# phi, theta = np.mgrid[0:pi:101j, 0:2 * pi:101j]
#
# x = r*sin(phi)*cos(theta)
# y = r*sin(phi)*sin(theta)
# z = r*cos(phi)
#
# mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 300))
# mlab.clf()
#
# spht = np.vectorize(sbm.sphtocart)
# print(type(pixs[cond].pc_lon))
# xss = np.zeros(len(pixs[cond].pc_lon))
# yss = np.zeros(len(pixs[cond].pc_lon))
# zss = np.zeros(len(pixs[cond].pc_lon))
# xss,yss,zss = spht(pixs[cond].pc_lon,pixs[cond].pc_lat)
#
# mlab.mesh(x , y , z, color=(0.0,0.5,0.5))
# mlab.points3d(xss, yss, zss, scale_factor=0.05)
#
# mlab.show()
