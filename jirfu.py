#!/usr/bin/python
# -*- coding: utf-8 -*-

import decimal
import numpy as np
import sys
import os.path
import matplotlib.pyplot as pl
import matplotlib.lines as lin
import math as m
from numpy import linalg as LA
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata


def stereoplot(lon,lat,col,nomefi,polo='N',min=0,max=0,title='',show=False):
    """
    Plots points on a stereographic map with colors.
    :return:
    """
    fig = pl.figure(figsize=(8, 6), dpi=150)
    pl.title(title)

    if polo == 'N':
        map = Basemap(projection='npstere',boundinglat=60,lon_0=180,resolution='l')
        map.drawparallels(np.arange(60,90,10))
    else:
        map = Basemap(projection='spstere',boundinglat=-60,lon_0=180,resolution='l')
        map.drawparallels(np.arange(-80,-50,10))

    aur_lon_0,aur_lat,aur_lon,aur_theta = leggi_map_aur(polo)

    map.drawmeridians(np.arange(0,360,30),labels=[1,1,0,1],fontsize=10)

    x, y = map(lon,lat)

    if(min == max):
        sca = map.scatter(x,y,c = col,s=5,edgecolors='none')
    else:
        sca = map.scatter(x,y,c = col,s=5,edgecolors='none',vmin=min,vmax=max)

    pl.colorbar()

    x, y = map(360-aur_lon,aur_lat)
    #map.scatter(x,y,color = 'white',edgecolors='black',s=15,marker = 'o')
    x = np.append(x,x[0])
    y = np.append(y,y[0])
    pl.plot(x,y,color='white',linewidth=2.0)
    pl.plot(x,y,color='black',linewidth=2.0,linestyle='--')
    if(show): pl.show()
    fig.savefig(nomefi, format='eps', dpi=150)
    pl.close()
    return


def stereomap(lon,lat,col,nomefi,polo='N',min=0,max=0,title='',show=False,lonlat=False,xres=50,lonres=180,latres=30,
              ncont=15,form=False,addpoints=False,divide=False):
    """
    Plots points on a stereographic map with colors.
    :return:
    """
    fig = pl.figure(figsize=(8, 6), dpi=150)
    pl.title(title)

    if polo == 'N':
        blat = 60
        map = Basemap(projection='npstere',boundinglat=blat,lon_0=180,resolution='l')
        map.drawparallels(np.arange(blat,90,10))
    else:
        blat = -60
        map = Basemap(projection='spstere',boundinglat=blat,lon_0=180,resolution='l')
        map.drawparallels(np.arange(-80,blat+10,10))

    aur_lon_0,aur_lat,aur_lon,aur_theta = leggi_map_aur(polo)

    map.drawmeridians(np.arange(0,360,30),labels=[1,1,0,1],fontsize=10)

    x, y = map(lon,lat)
    # print(x,y)
    # print(map(0.,-60.),map(90.,-60.),map(180.,-60.),map(270.,-60.))

    # Trovo i minimi e massimi delle coordinate del plot
    xu=np.zeros(4)
    yu=np.zeros(4)
    xu[0],yu[0] = map(0.,blat)
    xu[1],yu[1] = map(90.,blat)
    xu[2],yu[2] = map(180.,blat)
    xu[3],yu[3] = map(270.,blat)
    x0=np.min(xu)
    x1=np.max(xu)
    y0=np.min(yu)
    y1=np.max(yu)

    nx=xres
    ny=xres
    nstep=xres*1j
    xgri, ygri = np.mgrid[x0:x1:nstep, y0:y1:nstep]

    nsteplot=lonres*1j
    nsteplat=latres*1j
    if polo == 'N':
        grid_lat, grid_lon = np.mgrid[blat:90:nsteplat, 0:360:nsteplot]
    else:
        grid_lat, grid_lon = np.mgrid[-90:blat:nsteplat, 0:360:nsteplot]
    xg, yg = map(grid_lon, grid_lat)

    nlat = latres #np.shape(grid_lat)[0]
    nlon = lonres #np.shape(grid_lat)[1]
    steplo = 360/nlon
    stepla = 30/nlat

    # tolgo i NaN dai vettori
    cond2 = (~np.isnan(col))
    col2 = col[cond2]
    lon2 = lon[cond2]
    lat2 = lat[cond2]
    x2 = x[cond2]
    y2 = y[cond2]

    #grid_near = griddata((x,y), col, (xg, yg), method='nearest')
    #grid_near = griddata((lon,lat), col, (grid_lon, grid_lat), method='nearest')

    # Commento: i dati vanno grigliati a mano. Per ogni punto prendo quelli più vicini e li medio, buttando via i nan.

    lonlat=False

    if(lonlat):
        cols = -np.ones((nlat,nlon))
        num = np.zeros((nlat,nlon))
        for i in range(nlat):
            for j in range(nlon):
                glo = grid_lon[i,j]
                gla = grid_lat[i,j]
                cond = (lon2-glo >= -steplo/2) & (lon2-glo < steplo/2) & (lat2-gla >= -stepla/2) & (lat2-gla < stepla/2)
                if len(col2[cond]) > 0:
                    cols[i,j] = np.mean(col2[cond])
                num[i,j] = len(col2[cond])
    else:
        cols = -np.ones((nx,ny))
        num = np.zeros((nx,ny))
        steplo = xgri[1,0]-xgri[0,0]
        stepla = ygri[0,1]-ygri[0,0]
        for i in range(nx):
            for j in range(ny):
                glo = xgri[i,j]
                gla = ygri[i,j]
                cond = (x2-glo >= -steplo/2) & (x2-glo < steplo/2) & (y2-gla >= -stepla/2) & (y2-gla < stepla/2)
                if len(col2[cond]) > 0:
                    cols[i,j] = np.mean(col2[cond])
                num[i,j] = len(col2[cond])

    cols[(cols<0)]=float(np.nan)

    #pl.pcolormesh(xg, yg, cols)

    if lonlat:
        if(min == max):
            pl.contourf(xg, yg, cols,ncont)
        else:
            levels = np.linspace(min,max,ncont+1)
            cols[(cols > max)] = max
            cols[(cols < min)] = min
            pl.contourf(xg, yg, cols,levels=levels)
        #cs = pl.contour(xg, yg, cols,ncont)
        #pl.clabel(cs, inline=1,fontsize=10)#,manual=manual_locations)
    else:
        if(min == max):
            pl.contourf(xgri, ygri, cols,ncont)
        else:
            cols[(cols > max)] = max
            cols[(cols < min)] = min
            levels = np.linspace(min,max,ncont+1)
            if divide:
                pl.contourf(xgri, ygri, cols/np.sqrt(num),levels=levels)
            else:
                pl.contourf(xgri, ygri, cols,levels=levels)
        #cs = pl.contour(xgri, ygri, cols,ncont)
        #pl.clabel(cs, inline=1,fontsize=10)#,manual=manual_locations)

    if form:
        pl.colorbar(format='%.1e')
    else:
        pl.colorbar()

    if addpoints:
        if(min == max):
            print(np.min(cols[~np.isnan(cols)]),np.max(cols[~np.isnan(cols)]))
            sca = pl.scatter(x2,y2,c = col2,s=1,edgecolors='none',vmin=np.min(cols[~np.isnan(cols)]),vmax=np.max(cols[~np.isnan(cols)]))
        else:
            sca = pl.scatter(x2,y2,c = col2,s=1,edgecolors='none',vmin=min,vmax=max)

    #pl.colorbar()
    x, y = map(360-aur_lon,aur_lat)
    #map.scatter(x,y,color = 'white',edgecolors='black',s=15,marker = 'o')
    x = np.append(x,x[0])
    y = np.append(y,y[0])
    pl.plot(x,y,color='white',linewidth=2.0)
    pl.plot(x,y,color='black',linewidth=2.0,linestyle='--')
    if(show): pl.show()
    fig.savefig(nomefi, format='eps', dpi=150)
    pl.close()

    return


def stereopos(lon,lat,nomefi,color='black',marker='.',polo='N',title='',show=False):
    """
    Plots points on a stereographic map with colors.
    :return:
    """
    fig = pl.figure(figsize=(8, 6), dpi=150)
    pl.title(title)

    if polo == 'N':
        map = Basemap(projection='npstere',boundinglat=60,lon_0=180,resolution='l')
        map.drawparallels(np.arange(60,90,10))
    else:
        map = Basemap(projection='spstere',boundinglat=-60,lon_0=180,resolution='l')
        map.drawparallels(np.arange(-80,-50,10))

    aur_lon_0,aur_lat,aur_lon,aur_theta = leggi_map_aur(polo)

    map.drawmeridians(np.arange(0,360,30),labels=[1,1,0,1],fontsize=10)

    x, y = map(lon,lat)

    sca = map.scatter(x,y,s=1,marker=marker,color=color)
    # x = np.append(x,x[0])
    # y = np.append(y,y[0])
    # pl.plot(x,y)

    x, y = map(360-aur_lon,aur_lat)
    #map.scatter(x,y,color = 'white',edgecolors='black',s=15,marker = 'o')
    x = np.append(x,x[0])
    y = np.append(y,y[0])
    pl.plot(x,y,color='white',linewidth=2.0)
    pl.plot(x,y,color='black',linewidth=2.0,linestyle='--')
    if(show): pl.show()
    fig.savefig(nomefi, format='eps', dpi=150)
    pl.close()
    return


def findspi(wls,spe,thres=1.5,min=-2e-4):
    """
    Sets the mask at 0 if the point is a suspected spike. Checks line intensities on the H3+ line and each other point with more than thres the max of these lines is masked.
    :return: mask
    """

    lines = np.array([3315.0,3413,3531,3540,3665])
    vals = np.zeros(5)
    for line,i in zip(lines,range(5)):
        cond1 = (abs(wls-line) < 1.0)
        vals[i] = spe[cond1]

    cond = ((wls > 3200.) & (wls < 3800.)) & ((spe > thres*np.max(vals)) | (spe < min))
    mask = np.ones(len(spe), dtype = 'i4')
    mask[cond] = 0

    return mask


def checkqual(wls,spe,fondo,thres=-1e-3):
    """
    Checks the quality of the spectrum, counting elements lower than threshold with respect to the continuum.
    :return: 1 if ok, 0 in not
    """

    cond = (wls > 3250.) & (wls < 3700.) & (spe < fondo - thres)
    nu = sum(cond)

    mask = findspi(wls,spe)
    spi = len(mask[(mask == 0)])

    if(nu > 2 or spi > 2):
        ok = 0
    else:
        ok = 1

    return ok


def stereopolar(lon,lat,R=71000):
    """
    Calculates the stereographic projection, given lon, lat and radius.
    :return:
    """
    lamb = np.pi*0.5
    phi = np.pi
    lat = np.pi*lat/180.0
    lon = np.pi*lon/180.0

    k = 2*R/(1+m.sin(phi)*m.sin(lon)+m.sin(phi)*m.cos(lon)*m.cos(lat-lamb))
    x = k * m.cos(phi)*m.sin(lat-lamb)
    y = k * (m.cos(phi)*m.sin(lon)-m.sin(phi)*m.cos(lon)*m.cos(lat-lamb))

    return x,y


def leggi_map_aur(polo):
    """
    Legge le coordinate dell'ovale previsto dal modello VIP4.
    (Connerney, J.E.P., Açuna, M.H., Ness N.F., & Satoh, T. (1998).
    New models of Jupiter's magnetic field constrained by the Io Flux Tube footprint. J. Geophys. Res., 103, 11929 -11939)
    :param polo:
    :return:
    """
    if polo == 'S':
        filename = '/home/fede/Scrivania/Jiram/DATA/Model_VIP4/southR30table.txt'
    if polo == 'N':
        filename = '/home/fede/Scrivania/Jiram/DATA/Model_VIP4/northR30table.txt'

    infile = open(filename, 'r')
    infile.readline()
    infile.readline()
    data = [line.split() for line in infile]
    data_arr = np.array(data)
    aur_lon_0 = np.array([float(r) for r in data_arr[:, 0]])
    aur_lat = np.array([float(r) for r in data_arr[:, 1]])
    aur_lon = np.array([float(r) for r in data_arr[:, 2]])
    aur_theta = np.array([float(r) for r in data_arr[:, 3]])
    infile.close()

    return aur_lon_0,aur_lat,aur_lon,aur_theta


def fondojir(wl,spe):
    """
    Calculates the mean value of the background radiation in the spectrum.
    :param wl: Wavelengths
    :param spe: Spectrum
    :return:
    """
    cond = ((wl > 3475) & (wl < 3511)) | ((wl > 3558) & (wl < 3603)) | ((wl > 3735) & (wl < 3770))

    fondo = np.mean(spe[cond])

    return fondo


def ratio32_35(wl,spe,fondo):
    """
    Calculates ratio between H3+ lines.
    :param wl:
    :param spe:
    :return:
    """
    thres = 6e-4

    #cond1 = ((wl > 3195) & (wl < 3201))
    #cond2 = ((wl > 3518) & (wl < 3546))

    int1 = spe[133]-np.mean(spe[131:133])
    int1 = spe[133]+spe[134]-2*np.mean(spe[131:133])
    int2 = spe[171]+spe[170]-2*np.mean(spe[166:169])
    #int1 = spe[171] - fondo*2e3+(spe[170])
    int3 = spe[185]+spe[186]-2*np.mean(spe[182:185])
    ratio = int1/int2
    err = 3e-4/int1+3e-4/int2
    err = ratio*err
    if(int1 < thres or int2 < thres):
        ratio = float(np.nan)
        err = float(np.nan)

#    print(int1,int2,ratio,err)

    return ratio, err


def ratio35_37(wl,spe,fondo):
    """
    Calculates ratio between H3+ lines.
    :param wl:
    :param spe:
    :return:
    """
    thres = 6e-4

    #cond1 = ((wl > 3195) & (wl < 3201))
    #cond2 = ((wl > 3518) & (wl < 3546))

    int1 = spe[133]-np.mean(spe[131:133])
    int1 = spe[133]+spe[134]-2*np.mean(spe[131:133])
    int2 = spe[171]+spe[170]-2*np.mean(spe[166:169])
    #int1 = spe[171] - fondo*2e3+(spe[170])
    int3 = spe[185]+spe[186]-2*np.mean(spe[182:185])
    ratio = int2/int3
    err = 3e-4/int2+3e-4/int3
    err = ratio*err
    if(int2 < thres or int3 < thres):
        ratio = float(np.nan)
        err = float(np.nan)

#    print(int1,int2,ratio,err)

    return ratio, err


def ratio32_37(wl,spe,fondo):
    """
    Calculates ratio between H3+ lines.
    :param wl:
    :param spe:
    :return:
    """
    thres = 6e-4

    #cond1 = ((wl > 3195) & (wl < 3201))
    #cond2 = ((wl > 3518) & (wl < 3546))

    int1 = spe[133]-np.mean(spe[131:133])
    int1 = spe[133]+spe[134]-2*np.mean(spe[130:133])
    int2 = spe[171]+spe[170]-2*np.mean(spe[166:169])
    #int1 = spe[171] - fondo*2e3+(spe[170])
    int3 = spe[185]+spe[186]-2*np.mean(spe[182:185])
    ratio = int1/int3
    err = 3e-4/int2+3e-4/int3
    err = ratio*err
    if(int1 < thres or int3 < thres):
        ratio = float(np.nan)
        err = float(np.nan)

    return ratio, err


def ind_h3p(wl,spe,fondo):
    """
    Integrates the signal from H3+ lines.
    :param wl:
    :param spe:
    :return:
    """

    cond = ((wl > 3250) & (wl < 3275)) | ((wl > 3293) & (wl < 3330)) | ((wl > 3375) & (wl < 3463)) | ((wl > 3518) & (wl < 3546)) | ((wl > 3607) & (wl < 3630)) | ((wl > 3653) & (wl < 3681))

    intt = np.sum(spe[cond]-fondo)

    return intt


def integr_h3p(wl,spe,fondo,w1=3350,w2=3750):
    """
    Integrates the signal from H3+ lines.
    :param wl:
    :param spe:
    :return:
    """

    cond = (wl > w1) & (wl < w2)

    intt = np.trapz(spe[cond]*1e-3-fondo,x=wl[cond])

    return intt


def write_obs_JIR(freq,spe,mask,filename,comment=''):
    """
    Writes files of JIRAM observations. (JIRAM_MAP format)
    :return:
    """
    from datetime import datetime
    infile = open(filename, 'w')
    data = datetime.now()
    infile.write(comment+'\n')
    infile.write('\n')
    infile.write('Processed on: {}\n'.format(data))
    infile.write('\n')
    infile.write('Wavelength (nm), spectral data (W m^-2 um^-1 sr^-1), mask(0/1):\n')
    infile.write('{:1s}\n'.format('#'))

    for fr, ob, ma in zip(freq, spe, mask):
        if(np.isnan(ob)): ob = 0.0
        if(ma == 0): ob = 0.0
        infile.write('{:10.3f}{:15.5e}{:4d}\n'.format(fr,ob,ma))

    infile.close()
    return


def read_res_jir(filename):
    infile = open(filename, 'r')
    data = [line.split() for line in infile]
    data_arr = np.array(data)
    i = np.array([int(r) for r in data_arr[:, 0]])
    j = np.array([int(r) for r in data_arr[:, 1]])
    temp = np.array([float(r) for r in data_arr[:, 2]])
    err_t = np.array([float(r) for r in data_arr[:, 3]])
    infile.close()
    return i,j,temp,err_t


def read_res_jir_2(filename):
    infile = open(filename, 'r')
    data = [line.split() for line in infile]
    data_arr = np.array(data)
    i = np.array([int(r) for r in data_arr[:, 0]])
    j = np.array([int(r) for r in data_arr[:, 1]])
    temp = np.array([float(r) for r in data_arr[:, 2]])
    infile.close()
    return i,j,temp

def read_res_jir_3(filename):
    infile = open(filename, 'r')
    data = [line.split() for line in infile]
    data_arr = np.array(data)
    i = np.array([int(r) for r in data_arr[:, 0]])
    temp = np.array([float(r) for r in data_arr[:, 1]])
    err_t = np.array([float(r) for r in data_arr[:, 2]])
    infile.close()
    return i,temp,err_t


def read_res_jir_4(filename):
    infile = open(filename, 'r')
    data = [line.split() for line in infile]
    data_arr = np.array(data)
    i = np.array([int(r) for r in data_arr[:, 0]])
    temp = np.array([float(r) for r in data_arr[:, 1]])
    infile.close()
    return i,temp
