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

def download(TIC, verbose=True):

    
    target = f'TIC {TIC}'
    local_directory = os.path.join(data_path, f'{TIC}_tglc/')
    os.makedirs(local_directory, exist_ok=True)

    tglc_lc(target=target, 
            local_directory=local_directory, 
            size=90, 
            save_aper=False,
            limit_mag=16,
            get_all_lc=False,
            first_sector_only=False,
            sector=None,
            prior=None)  

    paths = [ os.path.join(local_directory, 'lc', f) for f in 
              os.listdir(os.path.join(local_directory, 'lc')) ]

    return tessutils.sort_paths(paths)

