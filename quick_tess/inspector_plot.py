#!/usr/bin/env python

import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.gridspec import GridSpec
import numpy as np
import os
import pandas as pd
import pickle

from . import data_parser
from . import variablestar
from . import plotutils
from . import tessutils

#Dom Rowan 2023

desc="""
Plot the sector-by-sector TESS LC
"""

def plot(path, parser, parser_args=None, savefig=None, savelc=None,
         period=None, double_period=True, period_pickle=None,
         **kwargs):

    '''
    Make the inspector plot for a TESS lc

    path (str or iterable) : path to TESS LC
    parser : function to convert to "standard" form
    parser_args : optional kwargs to pass to parser function
    savefig : output fname
    savelc: output LC in a format suitable for Period04
    period : pass a period (skips Lomb Scargle period search)
    double_period : double the period input by the LS (does nothing if period is not None)
    period_pickle : output name to save periods to

    kwargs : optional parameters to pass to periodogram call

    returns list of periods
    '''
    
    if not tessutils.check_iter(path):
        path = [ path ]

    if double_period:
        n = 2
    else:
        n = 1

    period_dict = {}

    vs = variablestar.VariableStar(path, parser=parser, parser_args=parser_args)
    if vs.df.empty:
        return
    if period is not None:
        vs.period = period
    else:
        vs.ls_periodogram(**kwargs)
        vs.phase_fold(n=n)

    if len(path) == 1:
        
        fig, ax = plt.subplots(2, 1, figsize=(12, 10))
        fig.subplots_adjust(top=.98, right=.98)
        ax[0] = vs.plot_full(time_offset=2.45e6, ax=ax[0], plot_kwargs=dict(rasterized=True))
        ax[1] = vs.plot_lc(nphase=2, ax=ax[1], plot_kwargs=dict(rasterized=True))
        ax[1].text(.95, 1.0, f'Sector {vs.df.sector.iloc[0]}',
                   ha='right', va='center', fontsize=20,
                   transform=ax[1].transAxes,
                   bbox=dict(edgecolor='black', facecolor='white'))

        period_dict[vs.df.sector.iloc[0]] = vs.period

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

        period_dict['all'] = vs.period

        for i in range(len(path)):
            vsi = variablestar.VariableStar(path[i], parser=parser, parser_args=parser_args)

            if vsi.df.empty:
                continue

            if period is not None:
                vsi.period = period
            else:
                vsi.ls_periodogram(**kwargs)
                vsi.phase_fold(n=n)

            period_dict[vsi.df.sector.iloc[i]] = vsi.period

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

    if period_pickle is not None:
        
        with open(period_pickle, 'wb') as p:
            pickle.dump(period_dict, p)
     
    return period_dict

    if savelc is not None:
                  vs.to_period04_format(savelc)
