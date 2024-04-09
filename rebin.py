#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 15:01:28 2024

@author: aellien
"""

import os
import sys
import glob
import numpy as np

def rebin(im, xbin = 2, ybin = 2, btype = 'sum'):

    print('Input shape: %d x %d' %( im.shape[0], im.shape[1] ))
    print('XBIN = %d, YBIN = %d' %( xbin, ybin ))

    xedge = np.shape(im)[0]%xbin
    yedge = np.shape(im)[1]%ybin
    im = im[xedge:,yedge:]
    binim = np.reshape(im,(int(np.shape(im)[0]/xbin),xbin,int(np.shape(im)[1]/ybin),ybin))

    if btype == 'mean':
        binim = np.mean(binim,axis=3)
        binim = np.mean(binim,axis=1)
    elif btype == 'sum':
        binim = np.sum(binim,axis=3)
        binim = np.sum(binim,axis=1)
    print('New shape: %d x %d' %( binim.shape[0], binim.shape[1]  ))

    return binim