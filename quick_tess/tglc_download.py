#!/usr/bin/env python

from astropy import units as u
from astropy.io import fits
from astropy.table import Table
import os

from tglc.quick_lc import tglc_lc
from . import tessutils

data_path = os.environ.get('QUICKTESS_DIR', None)

if data_path is None:
    print("Please set the QUICKTESS_DIR environment variable to path ot the downloaded data")

#Dom Rowan 2023

desc="""
Download data for a TIC using TGLC
"""

def download(TIC, verbose=True, **kwargs):

    
    target = f'TIC {TIC}'
    local_directory = os.path.join(data_path, f'{TIC}_tglc/')
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
        

