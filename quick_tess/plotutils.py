#!/usr/bin/env python

import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib
import numpy as np

#Dom Rowan 2023

desc="""
Utility plotting functions
"""

colors = ["#3696ff", "#f70065", "#011a7c", "#761954", "#8800b2"]

def plotparams(ax, labelsize=15):
    '''
    Basic plot params

    :param ax: axes to modify

    :type ax: matplotlib axes object

    :returns: modified matplotlib axes object
    '''
    ax.minorticks_on()
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(direction='in', which='both', labelsize=labelsize)
    ax.tick_params('both', length=8, width=1.8, which='major')
    ax.tick_params('both', length=4, width=1, which='minor')
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.5)
    return ax

def plotparams_cbar(cbar):

    cbar.ax.tick_params(direction='out', which='both', labelsize=15)
    cbar.ax.tick_params('y', length=8, width=1.8, which='major')
    cbar.ax.tick_params('y', length=4, width=1, which='minor')

    for axis in ['top', 'bottom', 'left', 'right']:
        cbar.ax.spines[axis].set_linewidth(1.5)

    return cbar

def plt_return(created_fig, fig, ax, savefig, dpi=300):
    if created_fig:
        if savefig is not None:
            fig.savefig(savefig, dpi=dpi)
            return 0
        else:
            plt.show()
            return 0
    else:
        return ax

def fig_init(ax=None, use_plotparams=True, figsize=(12,6),**kwargs):

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize, **kwargs)
        created_fig=True
    else:
        created_fig=False
        fig=None

    if use_plotparams and not has_twin(ax):
        ax = plotparams(ax)

    return fig, ax, created_fig

def has_twin(ax):
    for other_ax in ax.figure.axes:
        if other_ax is ax:
            continue
        if other_ax.bbox.bounds == ax.bbox.bounds:
            return True
    return False

def format_latex_label(label):

    return label.replace('_', ' ')

def many_colors():
    
    colors = [(25,25,25),(0,92,49),(43,206,72),(255,204,153),
              (148,255,181),(143,124,0),(157,204,0),
              (255,0,16),(94,241,242),(0,153,143),(224,255,102),(116,10,255),
              (153,0,0),(255,80,5),
              (194,0,136),(0,51,128),(255,164,5),(66,102,0),
              (240,163,255),(0,117,220),(153,63,0),(76,0,92)][::-1]

    return colors

def plotparams_bax(bax):

    for ax in bax.axs:
        ax.spines['top'].set_visible(True)
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(direction='in', which='both', labelsize=15)
        ax.spines['top'].set_linewidth(1.5)
        ax.spines['bottom'].set_linewidth(1.5)
    bax.axs[0].spines['left'].set_linewidth(1.5)
    bax.axs[-1].spines['right'].set_linewidth(1.5)
    bax.axs[-1].spines['right'].set_visible(True)
    bax.axs[-1].spines['right'].set_linewidth(1.5)
    bax.axs[-1].yaxis.set_ticks_position('right')
    bax.axs[-1].yaxis.set_ticklabels([])

    """
    for i in range(len(bax.axs)):
        bax.axs[i].spines['top'].set_visible(True)
        bax.axs[i].xaxis.set_ticks_position('both')
        if
    """


    for ax in bax.axs:
        ax.minorticks_on()
        ax.tick_params('both', length=8, width=1.8, which='major')
        ax.tick_params('both', length=4, width=1, which='minor')
        for l in ax.lines:
            if len(l._x) == 2:
                l.set_linewidth(1.5)

    for i in range(len(bax.axs)): 
        if (i != 0) and (i != len(bax.axs)-1): 
            bax.axs[i].tick_params(axis='y', which='minor', left=False) 
         
        bounds = bax.axs[i].get_position().bounds 
        size = bax.fig.get_size_inches() 
        ylen = bax.d*np.sin(bax.tilt*np.pi/180)*size[0]/size[1] 
        xlen = bax.d*np.cos(bax.tilt*np.pi/180) 
 
        d_kwargs=dict(transform=bax.fig.transFigure,  
                      color=bax.diag_color, clip_on=False, lw=1.5) 
 
        if i != len(bax.axs)-1: 
            xpos = bounds[0]+bounds[2] 
            ypos = bounds[1]+bounds[3] 
            bax.axs[i].plot((xpos-xlen, xpos+xlen), (ypos-ylen, ypos+ylen), **d_kwargs) 
 
        if i != 0: 
            xpos = bounds[0] 
            ypos = bounds[1]+bounds[3] 
            bax.axs[i].plot((xpos-xlen, xpos+xlen), (ypos-ylen, ypos+ylen), **d_kwargs) 
 
 
    return bax
