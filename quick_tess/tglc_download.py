#!/usr/bin/env python

from astropy import units as u
from astropy.io import fits
from astropy.table import Table
import os

from tglc.quick_lc import tglc_lc
from . import tessutils


#Dom Rowan 2023

desc="""
Download data for a TIC using TGLC
"""

def download(TIC, verbose=True, savedir=None, **kwargs):

    
    target = f'TIC {TIC}'
    
    if savedir is None:
        savedir = os.environ.get('QUICKTESS_DIR', None)
        if savedir is None:
            print("Set the QUICKTESS_DIR environment variable to path to\
                   the downloaded data")
            return
    local_directory = os.path.join(savedir, f'{TIC}_tglc/')
    os.makedirs(local_directory, exist_ok=True)

    kwargs.setdefault('size', 90)
    kwargs.setdefault('save_aper', False)
    kwargs.setdefault('limit_mag', 16)
    kwargs.setdefault('get_all_lc', False)
    kwargs.setdefault('first_sector_only', False)
    kwargs.setdefault('sector', None)
    kwargs.setdefault('prior', None)

    tglc_lc(target=target, 
            local_directory=local_directory, 
            **kwargs)

    paths = [ os.path.join(local_directory, 'lc', f) for f in 
              os.listdir(os.path.join(local_directory, 'lc')) ]

    paths = tessutils.sort_paths(paths)

    if len(paths) > 1 and (kwargs['sector'] is not None):
        return [ p for p in paths if tessutils.parse_path(p).sector == kwargs['sector'] ][0]
    else:
        return paths
        

