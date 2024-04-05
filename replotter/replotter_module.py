#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 15:03:07 2023

@author: kollar2


8-10-23 Making this able to load spec data too.

8-17-23 Found a bug that subtraction ind was not actually being set. Fixed it.



"""

import re
import random
# from scipy import *
import scipy
# import pylab
import matplotlib.pyplot as plt
import numpy as np
import time



import pickle
import datetime
import os
import sys

import copy

from scipy.ndimage import gaussian_filter





class replotter(object):
    def __init__(self, filepath,
                 vmin = -35,
                 vmax = -5,
                 cmap = 'YlGnBu_r',
                 div_cmap = 'bwr',
                 xlims = '',
                 ylims = '',
                 fig_xSize = 8,
                 fig_ySize = 6,
                 diff_cutoff = -35,
                 subtraction_ind = 0,
                 diff_colorBound = 10,
                 smoothe_data = False,
                 smoothing_sigma = 0.75,
                 verbose = False,
                 GHzVHz_correction = False):
        
        self.verbose = verbose
        
        #self.folder = folder
        #self.filename = filename
        self.filepath = filepath
        
        #color bar limits for transmission plots
        self.vmin = vmin
        self.vmax = vmax
        self.diff_cutoff = diff_cutoff
        self.diff_colorBound = diff_colorBound
        
        #settings for differential transmission
        self.subtraction_ind = subtraction_ind
        
        
        #color maps to use
        self.cmap = cmap
        self.div_cmap = div_cmap
        
        
        #smoothing settings
        self.smoothe_data = smoothe_data
        self.smoothing_sigma = smoothing_sigma
        
        #correction factor for old data that has frequncy axis in GHz
        self.GHzVHz_correction = GHzVHz_correction
        
        #load the data and determine the plot limits
        self.load_data()
        if xlims == '':
            #no limits given, determine automatically
            self.xlims = [self.freqs[0], self.freqs[-1]]
        else:
            self.xlims = xlims
            
        if ylims == '':
            #no limits given, determine automatically
            self.ylims = [self.voltages[0], self.voltages[-1]]
        else:
            self.ylims = ylims
            
        
        #default figure size and limits
        self.fig_xSize = fig_xSize
        self.fig_ySize = fig_ySize
        
        #trim the data and compute differential transmission
        self.compute_trimmed_data()
        self.compute_diff_transmission()
        
        if self.verbose:
            self.show_default_plots()
            
            
        
    def load_data(self):
        with open(self.filepath, 'rb') as file:
            pickledict = pickle.load(file)

        self.datadict = pickledict['Data']
        self.settings = pickledict['ExpSettings']
        
        if 'specdata' in self.datadict.keys():
            #this is spec data and I need to load differently
            specdict = self.datadict['specdata']
            
            self.voltages = self.datadict['voltages']
            self.labels = self.datadict['spec_labels']
            
            self.mags = specdict['mags']
            self.phases = specdict['phases']
            if self.GHzVHz_correction:
                self.freqs = specdict['xaxis']
            else:
                self.freqs = specdict['xaxis']/1e9
            
            self.mags = specdict['mags']
            self.phases = specdict['phases']
            #do the rolling background subtraction
            for vind in range(0, len(self.voltages)):
                self.mags[vind,:] = self.mags[vind,:] - np.mean(self.mags[vind,:])
                self.phases[vind,:] = self.phases[vind,:] - np.mean(self.phases[vind,:])
            
            
            
        else:
            self.data = self.datadict['full_data']
            self.voltages = self.datadict['voltages']
            # filename = datadict['filename']
            self.labels = self.datadict['labels']
            
            self.mags = self.data['mags']
            self.phases = self.data['phases']
            if self.GHzVHz_correction:
                self.freqs = self.data['xaxis']
            else:
                self.freqs = self.data['xaxis']/1e9
        
        if self.smoothe_data:
            raw_mags = copy.deepcopy(self.mags)
            self.mags = gaussian_filter(self.mags, sigma = self.smoothing_sigma)
            
        return
        
    def compute_trimmed_data(self):
        #zeroout the noisy stuff at low power
        self.trimmedMat = copy.deepcopy(self.mags)
        noise_region = np.where(self.trimmedMat < self.diff_cutoff)
        self.trimmedMat[noise_region] = self.diff_cutoff
        return
    
    def compute_diff_transmission(self):
        #take a reference trace
        self.ref_trace = self.trimmedMat[self.subtraction_ind,:]
        

        self.deltaMat = copy.deepcopy(self.trimmedMat) - self.ref_trace
        
        self.delta_max = np.max(self.deltaMat)
        self.delta_min = - np.max(self.deltaMat) #np.min(deltaMat)
        return
    
    def show_baseline_plot(self, fignum = 1):
        ####plot the raw data, and the cutoff for differential transmission
        fig1 = plt.figure(fignum)
        plt.clf()
        
        ax = plt.subplot(1,2,1)
        plt.pcolormesh(self.freqs, self.voltages, self.mags, cmap = 'hot',shading = 'auto')
        plt.xlabel(self.labels[0])
        plt.ylabel(self.labels[1])
        plt.title('Raw Mag Data')
        plt.colorbar()
        
        ax = plt.subplot(1,2,2)
        plt.plot(self.freqs, self.mags[-1,:], color = 'mediumblue', label = 'data')
        plt.plot(self.freqs, self.diff_cutoff * np.ones(len(self.freqs)), color = 'firebrick', label = 'cutoff')
        plt.xlabel(self.labels[0])
        plt.ylabel('S21 Log Mag (dB)')
        plt.title('Single Trace')
        ax.legend(loc = 'lower right')
        
        fig1.set_size_inches([14, 6.5])
        plt.tight_layout()
        plt.show()
        # save_name = filename + '_Replot.png'
        # fig1.savefig(save_name, dpi = 400)
        ######################## 
        return
    
    def show_transmission_plot(self, fignum = 2, savefig = False):
        fig2 = plt.figure(fignum)
        plt.clf()
        
        ax = plt.subplot(1,1,1)
        plt.pcolormesh(self.freqs, self.voltages, self.mags, 
                       cmap = self.cmap,
                       vmin = self.vmin,
                       vmax = self.vmax,
                       shading = 'auto')
        plt.xlabel(self.labels[0])
        plt.ylabel(self.labels[1])
        plt.title('Transmission v Applied Flux')
        plt.colorbar()
        ax.set_xlim(self.xlims)
        ax.set_ylim(self.ylims)
        
        plt.suptitle(self.filename)
        
        fig2.set_size_inches([self.fig_xSize, self.fig_ySize])
        plt.tight_layout()
        plt.show()
        
        if savefig:
            save_name = self.filename + '_' + self.cmap + '.png'
            fig2.savefig(save_name, dpi = 500)
            # save_name = filename + '_' + cmap + '.pdf'  #WARNING - PDF save seems to crash. Maybe map is too big
            # fig22.savefig(save_name)
        
        return
    
    def show_differential_plot(self, fignum = 3, savefig = True):
        #plot the subtracted data, but using deltMat wich zeros
        #out parts of the data that are too low to have signal
        fig3 = plt.figure(fignum)
        plt.clf()
        
        ax = plt.subplot(1,1,1)
        plt.pcolormesh(self.freqs, self.voltages, self.deltaMat, 
                       cmap = self.div_cmap, 
                       vmax = self.diff_colorBound, 
                       vmin = -self.diff_colorBound,
                       shading = 'auto')
        plt.xlabel(self.labels[0])
        plt.ylabel(self.labels[1])
        plt.title('Differential Transmission v Applied Flux (cutoff, log mode)')
        plt.colorbar()
        ax.set_xlim(self.xlims)
        ax.set_ylim(self.ylims)
        
        plt.suptitle(self.filename)
        
        # fig23.set_size_inches([7.5, 2.75])
        fig3.set_size_inches([self.fig_xSize, self.fig_ySize])
        plt.tight_layout()
        plt.show()
        
        if savefig:
            save_name = self.filename + '_Differential.png'
            fig3.savefig(save_name, dpi = 500)
            # save_name = filename + '_Differential.pdf' #WARNING - PDF save seems to crash. Maybe map is too big
            # fig23.savefig(save_name)
        
        return
        
    
    def show_default_plots(self, savefig = False):
        
        self.show_baseline_plot(fignum =1)
        
        self.show_transmission_plot(fignum = 2, savefig = savefig)       
        
        self.show_differential_plot(fignum = 3, savefig = savefig)
        
        return
    
    
    
    
    
    
    
    