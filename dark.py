#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os.path
import matplotlib.pyplot as pl
import math as mt
import spect_base_module as sbm
from subprocess import call
import pickle
import scipy.io as io
from scipy.interpolate import PchipInterpolator as spline


lin1 = 30
lin2 = 48
sam1 = 6
sam2 = 30

cart = '/home/fede/Scrivania/Jiram/ANALISI/TEST_1/'

cubo = io.readsav(cart+'cubo_test.sav')
data = cubo.spe

print(type(data))

wew = io.readsav(cart+'wls.sav')
wl = wew.wls

#data[(data < 0)] = np.nan


mea = np.mean(data[:,lin1:lin2,sam1:sam2], axis = 1)
mea2 = np.mean(mea, axis = 1)



fig = pl.figure(figsize=(8, 6), dpi=150)
pl.title('Test data')
pl.imshow(data[172,20:55,50:240], cmap='jet', interpolation="nearest")
pl.colorbar(orientation="horizontal")
#pl.show()

fig.savefig(cart+'DATA.eps', format='eps', dpi=150)
pl.close()
