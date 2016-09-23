#!/usr/bin/python
# -*- coding: utf-8 -*-

import decimal
import numpy as np
import sys
import os.path
import matplotlib.pyplot as pl
import math as m
from numpy import linalg as LA


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


def checkqual(wls,spe,thres=-7e-4):
    """
    Checks the quality of the spectrum, counting elements lower than threshold.
    :return: 1 if ok, 0 in not
    """

    cond = (wls > 3250.) & (wls < 3700.) & (spe < thres)
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
