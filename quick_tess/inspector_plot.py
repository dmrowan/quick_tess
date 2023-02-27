#!/usr/bin/env python

import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.gridspec import GridSpec
import numpy as np
import os
import pandas as pd

from . import data_parser
from . import variablestar
from . import plotutils
from . import tessutils

#Dom Rowan 2023

desc="""
Plot the sector-by-sector TESS LC
"""

def plot(path, parser, savefig=None, parser_args=None, 
         period=None, **kwargs):
    
    if not tessutils.check_iter(path):
        path = [ path ]

    vs = variablestar.VariableStar(path, parser=parser, parser_args=parser_args)
    if period is not None:
        vs.period = period
    else:
        vs.ls_periodogram(**kwargs)

    if len(path) == 1:
        
        fig, ax = plt.subplots(1, 1, figsize=(12, 10))
        fig.subplots_adjust(top=.98, right=.98)
        ax[0] = vs.plot_full(time_offset=2.45e6, ax=ax[0], plot_kwargs=dict(rasterized=True))
        ax[1] = vs.plot_lc(nphase=2, ax=ax[1], plot_kwargs=dict(rasterized=True))
        ax[1].text(.95, 1.0, f'Sector {vs.df.sector.iloc[0]}',
                   ha='right', va='center', fontsize=20,
                   transform=ax[1].transAxes,
                   bbox=dict(edgecolor='black', facecolor='white'))

    else:
        
        fig = plt.Figure(figsize=(12, 2+4*(2+len(path))))
        fig.subplots_adjust(top=.98, right=.98, bottom=0.05)

        gs = GridSpec(2+len(path), 1)

        bax = vs.plot_broken(fig=fig, gs=gs[0], time_offset=2.45e6, 
                             plot_kwargs=dict(rasterized=True))
        bax.axs[0].set_ylabel('TESS T Mag', fontsize=20)
        ax1 = plt.Subplot(fig, gs[1])
        ax1 = vs.plot_lc(nphase=2, ax=ax1, plot_kwargs=dict(rasterized=True))
        fig.add_subplot(ax1)

        ax1.text(.95, 1.0, f'All Sectors',
                 ha='right', va='center', fontsize=20,
                 transform=ax1.transAxes,
                 bbox=dict(edgecolor='black', facecolor='white'))


        for i in range(len(path)):
            vsi = variablestar.VariableStar(path[i], parser=parser, parser_args=parser_args)

            if period is not None:
                vsi.period = period
            else:
                vsi.ls_periodogram(**kwargs)

            axi = plt.Subplot(fig, gs[i+2])
            axi = vsi.plot_lc(nphase=2, ax=axi, plot_kwargs=dict(rasterized=True))
            axi.text(.95, 1.0, f'Sector {vsi.df.sector.iloc[0]}',
                     ha='right', va='center', fontsize=20,
                     transform=axi.transAxes,
                     bbox=dict(edgecolor='black', facecolor='white'))
            fig.add_subplot(axi)

    if savefig is None:
        plt.show()
    else:
        fig.savefig(savefig, dpi=300)
