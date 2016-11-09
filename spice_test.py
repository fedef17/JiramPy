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
import spiceypy.spiceypy as spice
#from plotly.offline import download_plotlyjs, init_notebook_mode, iplot
#from plotly.graph_objs import *

#spice.tkvrsn("TOOLKIT")

spice.furnsh('/home/fede/Scrivania/Jiram/DATA/KERNELS_JIRAM/Kernels_jm0003/jm0003.mk')

print('ciao')
body = 'jupiter'

# questo legge i raggi del pianeta
n, rad = spice.bodvrd(body,'RADII',3)
print(rad)

# questo converte il tempo (la stringa start time dell'osservazione) in et:
tempo = '2016-08-27T14:55:58.817'
et = spice.str2et(tempo)
print(et)
tempo2 = '2016-08-27T14:55:59.817'
et2 = spice.str2et(tempo2)
print(et2)
#

targ = spice.bods2c('JUNO') # converte il nome in un numero intero identificativo
ref = 'IAU_JUPITER'
obs = spice.bods2c('jupiter')
pos,lt = spice.spkgps(targ, et, ref, obs)
print(pos)
pos2,lt2 = spice.spkgps(targ, et2, ref, obs)
print(pos2)

# Trova il boresight
et3 = spice.str2et('2016-08-26T17:23:47.817')
method = 'Ellipsoid'
target = 'JUPITER'
abcorr = 'CN+S' #"CN+S"
obs = 'JUNO'
dref = 'JUNO_JIRAM_S'
dvec = np.array([0.,0.,1.])
duduu = spice.sincpt(method, target, et3, ref, abcorr, obs, dref, dvec)
print(duduu)

# Per altri punti ti sposti lungo x di un ifov per ogni pixel (praticamente x è in gradi)
# il vettore dvec deve essere normalizzato a 1 (z serve a quello) ifov = 0.237767e-3 (mrad)
#init_notebook_mode()

