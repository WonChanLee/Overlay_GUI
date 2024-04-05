#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 11:37:12 2023

@author: kollar2
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, CheckButtons
import scqubits as scq

import copy




def stitch_data_files(ax, replot_list,
                     convert_to_flux=False,
                     offset=0, 
                     volts_per_flux=1,
                     use_phase_data = False,
                     freq_on_x_axis = True):
    '''
    Takes a list of files and combines them into on big colormap plot
    The code will subtract the mean from every row of the data to smooth out 
    variations in the transmission plot
    The xaxis can be volts or flux quanta
    Params:
        ax: axis of the figure where plots should be made
        files: list of replot classes for the data files
        (ordered by the user to get best results) of the colormaps we want to plot
        saveDir: directory containing all the files
        convert_to_flux: boolean to toggle the x axis
        offset: offset voltage corresponding to integer flux
        volts_per_flux: conversion between voltage and flux quanta
        
        if phase = True, then use the flux data
        
    This is modified from stitch_colormaps in the old fluxonium gui 
    New version is intended for lattice project data, which has a replotter
    class which stores lots of useful plot settings associated with each data file
    Will make use of that.    
        
    '''
    for curr_scan in replot_list:
        #print(curr_scan.filename)
        
        voltages = curr_scan.voltages
        
        # specdata = data['specdata']
        mags  = curr_scan.mags
        phases = curr_scan.phases
        if use_phase_data:
            thedata = phases
        else:
            thedata = mags
            
        freqs = curr_scan.freqs
        
        # for ind, mag in enumerate(thedata):
        #     thedata[ind] = mag - np.mean(mag)
        
        
        if convert_to_flux:
            #print('flux')
            voltaxis = volt_to_flux(voltages, offset, volts_per_flux)
        else:
            voltaxis = voltages
            
        #configure settings to view with voltage on x or y depending on aspect ratio of the data
        if not freq_on_x_axis:
            thedata = np.transpose(thedata)
            xaxis = voltaxis
            yaxis = freqs
            if convert_to_flux:
                xlabel = 'Flux (phi_ext/phi_0)'
            else:
                xlabel = 'Voltage (V)'
            ylabel = 'Frequency (GHz)' 
        else:
            yaxis = voltaxis
            xaxis = freqs
            if convert_to_flux:
                ylabel = 'Flux (phi_ext/phi_0)'
            else:
                ylabel = 'Voltage (V)'
            xlabel = 'Frequency (GHz)'
            
        
        #rescale all the data so that it can be plotted on the same color scale
        vmax = curr_scan.vmax
        vmin = curr_scan.vmin
        
        rescaledData = copy.deepcopy(thedata)
        rescaledData = (rescaledData-vmin)/(vmax-vmin) #hopefull this is now scaled to run between 0 and 1
            
        im=ax.pcolormesh(xaxis[:rescaledData.shape[1]], yaxis, 
                         rescaledData, 
                         cmap=curr_scan.cmap,
                         vmin=0, 
                         vmax=1,
                         shading = 'auto')

        ax.set_ylabel(ylabel)

        ax.set_xlabel(xlabel)
    return
        


def volt_to_flux(volts, offset, volts_per_flux):
    '''
    Helper function to convert volts to flux
    Params:
        volts: applied external voltage (double or array of doubles)
        offset: voltage where we have integer flux (double)
        volts_per_flux: conversion from volts to flux quantum (double)
    Returns:
        flux or flux array rounded to 6 places
        
        changed to 6 decimals because of ultrfine lattice spec data AK 8-10-23
    '''
        
    return np.round((volts-offset)/volts_per_flux, 6)

def flux_to_volt(flux, offset, volts_per_flux):
    '''
    Helper function to convert volts to flux
    Params:
        flux: external flux (double or array of doubles)
        offset: voltage where we have integer flux (double)
        volts_per_flux: conversion from volts to flux quantum (double)
    Returns:
        volts or volts array rounded to 6 places
        
        changed to 6 decimals because of ultrfine lattice spec data AK 8-10-23
    '''
    return np.round(flux*volts_per_flux + offset, 6)





#############
##stuff for calculating transmon properties


q_e = 1.602e-19  # electron charge in C
h = 6.626e-34  # Planck constant in Js
hbar = h / (2 * np.pi)  # Planck constant / 2p in Js
phi_0 = hbar / (2 * q_e)  # flux quantum / 2p in Wb



def fluxonium_energy_single_flux(phi, EC, EJ, EL, n=3):
    fx = scq.Fluxonium(EJ = EJ, EC = EC, EL = EL, flux = phi, cutoff = 110)
    return fx.eigenvals(evals_count = n)


def fluxonium_energies(EC, EJ, EL,
                      min_flux = -0.6, 
                      max_flux = 0.6, 
                      flux_points = 51,
                      n = 3):
    '''
    Computes the spectrum of the transmon specified by EC, EJ1, EJ2 over the 
    flux range [min_flux, max_flux] with flux_points points
    Params:
        EC (GHz)
        EJ_sum (GHz)
        EJ_diff (GHz)
        min_flux: min ext flux/2pi 
        max_flux: max ext flux/2pi
        flux_points: positive int (should be odd to avoid weird aliasing)
        
        n: number of transmon levels kept
    Returns:
        flux: np array of fluxes
        energies: matrix holding the computed eigen vals over flux range
    '''
    flux = np.linspace(min_flux, max_flux, flux_points)
    # (evals_s, _) = flux_sweep(flux*2*np.pi, EC, EJ, EL)
    energies = np.zeros((n, len(flux)))
    for find in range(0, len(flux)):
        phi = flux[find]
        energies[:,find] = fluxonium_energy_single_flux(phi = phi, EC = EC, EJ = EJ, EL = EL, n=n)
    
    return energies, flux





##########
#functions for the overlay gui


mastercolors = {}
mastercolors['e1'] = 'red'
mastercolors['e2'] = 'gold'#'goldenrod'

# mastercolors['g2'] = 'dodgerblue'
# mastercolors['e1'] = 'limegreen'
# mastercolors['e2'] = 'orange'#'goldenrod'
# mastercolors['e3'] = 'red'#'indianred'#'cyan'
# mastercolors['e4'] = 'fuchsia'#'mediumorchid'#hotpink'


def plot_line(ax, flux, energy,
              linewidth=3,
              alpha=0.5,
              color='red',
              linestyle='-',
              marker='o',
              label='default',
              freq_on_x_axis = True):
    '''
    Standard line plotting function to keep everything consistent and uncluttered late in the code
    
    added option to flip x and y axes of the data depending on which aspect ratio is best for viewing AK 8-10-23
    
    '''
    if freq_on_x_axis:
        xvals = energy
        yvals = flux
    else:
        xvals = flux
        yvals = energy
    
    
    
    [line, ] = ax.plot(
        xvals,
        yvals,
        linewidth=linewidth,
        alpha=alpha,
        color=color,
        linestyle=linestyle, 
        marker=marker,
        label = label,
        markersize=5
        )
    return line

def bare_transitions(ax, lines, full_energies, flux, freq_on_x_axis = True):
    '''
    Plot the n lowest energy levels of the transmon with the groundstate subtracted
    '''
    labels = ['e1', 'e2']

    colors = [mastercolors['e1'], mastercolors['e2']]
    marker = 'o'
    energies = full_energies[1:] - full_energies[0]
    for energy, label, color in zip(energies, labels, colors):
        line = plot_line(ax, flux, energy, color=color, label=label, marker=marker,
                         freq_on_x_axis = freq_on_x_axis)
        lines.append(line)
    ax.legend(loc='lower right')

def two_photon(ax, lines, full_energies, flux, freq_on_x_axis = True):
    '''
    Plot two photon transitions (essentially higher energy/2)
    '''
    labels = ['e1/2', 'e2/2']
    colors = [mastercolors['e1'],mastercolors['e2']]
    marker = 'd'
    linestyle = '--'
    energies = (full_energies[1:] - full_energies[0])/2
    for energy, label, color in zip(energies, labels, colors):
        line = plot_line(ax, flux, energy, color=color, label=label, marker=marker, linestyle = linestyle,
                         freq_on_x_axis = freq_on_x_axis)
        lines.append(line)
    ax.legend(loc='lower right')


def cavity_assisted(ax, lines, full_energies, flux, cavity_freq, freq_on_x_axis = True):
    '''
    Plot the cavity assisted transitions
    '''
    labels = ['e2 cav assist']
    colors = [mastercolors['e2']]
    marker = 's'
    linestyle = '-.'
    energies = np.abs((full_energies[2:] - full_energies[0]) - cavity_freq)
    for energy, label, color in zip(energies, labels, colors):
        line = plot_line(ax, flux, energy, color=color, label=label, marker=marker, linestyle=linestyle,
                         freq_on_x_axis = freq_on_x_axis)
        lines.append(line)
    ax.legend(loc='lower right')
    
def cavity_line(ax, lines, flux, cav_freq, freq_on_x_axis = True):
    '''
    Plot the cavity line itself
    '''
    line = plot_line(ax, flux, cav_freq*np.ones(len(flux)), 
                     color='k', 
                     label='cavity', 
                     marker='', 
                     linewidth = 2,
                     alpha = 0.5,
                     freq_on_x_axis = freq_on_x_axis)
    lines.append(line)
    ax.legend(loc='lower right')




# check box
def fluxonium_plot_levels(ax, lines, energies, flux, cav_freq,
                         bare=True, 
                         twoPhoton=False, 
                         cavAssist=False,
                         freq_on_x_axis = True):
    '''
    Plot all the user specified lines from provided spectral data:
        lines: list of lines to keep track of so we can remove them before plotting new ones
        energies: matrix of the n lowest energies that are being considered
        flux: np array of flux values at which the energies were computed
        cav_freq: frequency of the cavity in GHz
        bare: boolean indicating if the bare transistions should be plotted
        twoPhoton: boolean indicating if the two photon transistions should be plotted
        cavAssist: boolean indicating if the cavity assisted transistions should be plotted
        
        freq_on_x_axis: boolean to set the orientation of the plot

    '''
    #for line in lines:
    #    try:
    #        line.remove()
    #    except:
    #        continue
    #lines = []
    
    if bare:
        bare_transitions(ax, lines, energies, flux,
                         freq_on_x_axis = freq_on_x_axis)
    if twoPhoton:
        two_photon(ax, lines, energies, flux,
                         freq_on_x_axis = freq_on_x_axis)
    if cavAssist:
        cavity_assisted(ax, lines, energies, flux, cav_freq,
                         freq_on_x_axis = freq_on_x_axis)
    
    cavity_line(ax, lines, flux, cav_freq,
                         freq_on_x_axis = freq_on_x_axis)
    return lines






#### the main overlay gui itself
######################

def overlay_GUI_fluxonium(fig, ax, lines, 
                         EC, EJ, EL, 
                         cav_freq = 40, 
                         min_flux=-0.5, 
                         max_flux=0.5, 
                         flux_points=51, 
                         bare=True, 
                         twoPhoton=False, 
                         cavAssist=False,
                         freq_on_x_axis = True,
                         ):
    '''
    GUI allowing the user to tweak the transmon values to get the best data fit by eye
    There are also options to change which transitions are plotted to minimize 
    clutter in the early fitting stages
    Params:
        fig: reference figure
        ax: reference axes
        UI: list containing pointers to all the widget objects (important to keep them active)
        
        EC (GHz)
        EJ (GHz)
        
        cav_freq (GHz)
        
        min_flux (phi/2pi)
        max_flux (phi/2pi)
        flux_points (int)
        
        bare, twoPhoton, cavAssist: booleans controlling 
        which transitions are plotted
        
        freq_on_x_axis: boolean to set the orientation of the plot
    '''
    # Get the starting energy spectrum and plot the lines
    # lines, energies and flux are variables used throughout the code so need to be global (I think?)
    # global lines
    global energies
    global flux
    #global lines
    
    #lines = []
    #for line in lines:
    #    lines.remove()
    energies, flux = fluxonium_energies(EC, EJ, EL, min_flux, max_flux, flux_points)
    print(energies)
    lines = fluxonium_plot_levels(ax, lines, energies, flux, 
                                 cav_freq, 
                                 bare = bare,
                                 twoPhoton = twoPhoton,
                                 cavAssist = cavAssist,
                                 freq_on_x_axis = freq_on_x_axis)
    
    # #draw in the cavity
    # line = plot_line(ax, flux, cav_freq*np.ones(len(flux)), 
    #                  color='k', 
    #                  label='cavity', 
    #                  marker='', 
    #                  linewidth = 2,
    #                  alpha = 0.75,
    #                  freq_on_x_axis = freq_on_x_axis)
    
    # Adjust the subplots region to leave some space for the sliders and buttons
    #fig.subplots_adjust(left=0.1, bottom=0.35)
    axis_color = 'lightgoldenrodyellow'
    #fig.canvas.update()
    fig.canvas.draw_idle()