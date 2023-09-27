#!/usr/bin/env python

from collections import namedtuple
import os
import numba
import numpy as np
import pandas as pd
from tqdm import tqdm
import sys


#Dom Rowan 2023

desc="""
Utility functions for quick_tess 
"""

path_tup = namedtuple('path_tup', ['path', 'tic', 'sector', 'lctype', 'gaia'])

#Determine LCtype, sector tic from lc path
def parse_path(lc_path):

    if 'hlsp_qlp_tess' in lc_path:
        lctype='qlp'
        tic = int(os.path.split(lc_path)[-1].split('_')[4].split('-')[1])
        sector = int(os.path.split(
                lc_path)[-1].split('_')[4].split('-')[0].lstrip('s'))
        gaia = None
    elif 'hlsp_tess-spoc' in lc_path:
        lctype='spoc'
        tic = int(os.path.split(lc_path)[-1].split('_')[4].split('-')[0])
        sector = int(os.path.split(
                lc_path)[-1].split('_')[4].split('-')[1].lstrip('s'))
        gaia = None
    elif 'hlsp_tglc_tess' in lc_path:
        lctype = 'tglc'
        tic = None
        gaia = os.path.split(lc_path)[-1].split('_')[4].split('-')[2]
        sector = int(os.path.split(lc_path)[-1].split('_')[4].split('-')[2].lstrip('s'))
    else:
        raise ValueError(f'invalid lc path format {lc_path}')

    return path_tup(lc_path, tic, sector, lctype, gaia)

def sort_paths(paths):

    if not check_iter(paths):
        return paths

    sectors = [ parse_path(p).sector for p in paths ]

    df = pd.DataFrame({'path':paths, 'sector':sectors})
    df = df.sort_values(by='sector', ascending=True).reset_index(drop=True)

    return df.path.to_numpy()


#Check if list tuple, np array, etc
def check_iter(var):
    if isinstance(var, str):
        return False
    elif hasattr(var, '__iter__'):
        return True
    else:
        return False

#utility function for reading in df from various extensions
def pd_read(table_path, low_memory=False):

    if (not isinstance(table_path, pd.core.frame.DataFrame)
        and check_iter(table_path)):
        df0 = pd_read(table_path[0])
        for i in range(1, len(table_path)):
            df0 = pd.concat([df0, pd_read(table_path[i])])
        return df0
    else:
        if type(table_path) == pd.core.frame.DataFrame:
            return table_path
        elif table_path.endswith('.csv') or table_path.endswith('.dat'):
            return pd.read_csv(table_path, low_memory=low_memory)
        elif table_path.endswith('.pickle'):
            return pd.read_pickle(table_path)
        else:
            raise TypeError("invalid extension")

def pd_write(df, table_path, index=False):
    
    if table_path.endswith('.csv'):
        df.to_csv(table_path, index=index)
    elif table_path.endswith('.pickle'):
        df.to_pickle(table_path, index=index)
    else:
        raise TypeError("invalid extension for ellutils pd write")

#decorator to create a multiprocessing list to store return vals from function 
def manager_list_wrapper(func, L, *args, **kwargs):
    return_vals = func(*args, **kwargs)
    print(return_vals)
    L.append(return_vals)
    return return_vals

def manager_list_wrapper_silent(func, L, *args, **kwargs):
    return_vals = func(*args, **kwargs)
    L.append(return_vals)
    return return_vals

class HiddenPrints:
    '''
    From https://stackoverflow.com/questions/8391411/how-to-block-calls-to-print
    '''
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


#Used for plot_lc legend
def sort_column(arr):

    try:
        arr.sort()
        return arr
    except:
        arr_numeric = [
                float(''.join([char for char in s if char.isnumeric()]))
                if any([char.isnumeric() for char in s])
                else s for s in arr ]
        df = pd.DataFrame({'original':arr, 'numeric':arr_numeric})
        df = df.sort_values(by='numeric')
        return list(df.original)

def extend_phase(xvals, yvals, nphase, xtype='phase'):

    if xtype == 'phase':
        extended_xvals = np.concatenate([xvals+i for i in range(nphase)])
    elif xtype == 'phi':
        extended_xvals = np.concatenate([xvals+2*np.pi*i for i in range(nphase)])
    else:
        raise ValueError('xtype must be phase or phi')
                
    if len(yvals.shape) == 1:
        extended_yvals = np.tile(yvals, nphase)
    else:
        extended_yvals = np.zeros((yvals.shape[0], yvals.shape[1]*2), dtype=float)
        for i in range(yvals.shape[0]):
            extended_yvals[i] = np.tile(yvals[i],2)

    return extended_xvals, extended_yvals

def filter_format(s):

    if ':' in s:
        s = s.split(':')[1]

    if '_' not in s:
        return s
    else:
        return r'$'+s+r'$'

#requires nparrays
def compute_chi2(yvals, ymodel, yerr):
    return np.sum( (yvals-ymodel)**2 / (yerr**2))

def find_runs(lst):
    runs = []
    run = []
    for i in range(len(lst)):
        if i == 0 or lst[i] != lst[i-1] + 1:
            if len(run) > 1:
                runs.append(run)
            run = []
        run.append(i)
    if len(run) > 1:
        runs.append(run)
    return runs

@numba.jit(nopython=True)
def pdm_evaluate(times, yvals, period):
    
    phase = np.mod(times, period)/period

    bins = np.arange(-0.05, 1.05, 0.05)

    bin_numbers = np.digitize(phase, bins)
    bin_stds = np.zeros(len(bins)-1)

    for i in range(1, len(bins)):
        bin_data = yvals[bin_numbers == i]
        if len(bin_data):
            bin_stds[i-1] = np.std(bin_data)
        else:
            bin_stds[i-1] = 0

    variance = np.sum(np.power(bin_stds,2))

    return variance

def pdm(times, yvals, period_arr, progress=True):
    
    t0 = np.min(times)
    times = times-t0

    variance = np.zeros(len(period_arr))
    if progress:
        iterable = tqdm(range(len(period_arr)))
    else:
        iterable = range(len(period_arr))
    for i in iterable:
        variance[i] = pdm_evaluate(times, yvals, period_arr[i])

    return variance

