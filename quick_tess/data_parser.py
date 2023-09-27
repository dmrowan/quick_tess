#!/usr/bin/env python

from astropy import log
from astropy.io import fits
from astropy.table import Table
from astropy.units import UnitsWarning
import numpy as np
import pandas as pd
import os

import warnings
warnings.filterwarnings("ignore", message="invalid value encountered in log10")
warnings.simplefilter('ignore', category=UnitsWarning)

from . import tessutils

desc="""
Parser for different types of data that can be read into VariableStar
"""

#TESS SPOC LCs from MAST
def tess_spoc(datafile_list, pdc=False):

    if not tessutils.check_iter(datafile_list):
        datafile_list = [datafile_list]

    df_list = []
    for datafile in datafile_list:

        tab = Table.read(datafile, hdu=1)
        df = tab.to_pandas()
        df['hjd'] = df.TIME + 2457000

        if pdc:
            df['mag'] = -2.5*np.log10(df.PDCSAP_FLUX)
            df['mag_err'] = abs((-2.5*np.log10(np.e))*(df.PDCSAP_FLUX_ERR/df.PDCSAP_FLUX))
        else:
            df['mag'] = -2.5*np.log10(df.SAP_FLUX)
            df['mag_err'] = abs((-2.5*np.log10(np.e))*(df.SAP_FLUX_ERR/df.SAP_FLUX))
        df['sector'] = fits.open(datafile)[0].header['SECTOR']
        df = df.rename({'QUALITY':'quality'}, axis=1)
        df = df[~np.isnan(df['mag'])]
        df = df[df['quality'] == 0]
        tess_mag = fits.open(datafile)[0].header['TESSMAG']
        df['mag'] = df['mag'] + (tess_mag - df.mag.median())
        df['filter'] = 'T'
        df_list.append(df)
        #Add orbit column
        idx = np.argmax(df.hjd.diff())
        orbit = np.full(len(df), df.sector.iloc[0]*2+7, dtype=float)
        orbit[idx:] += 0.5
        df['orbit'] = orbit

    df_out = pd.concat(df_list)
    df_out = df_out.reset_index(drop=True)

    return df_out

#TESS QLP LCs from MAST
def tess_qlp(datafile_list, ksp=False):

    if not tessutils.check_iter(datafile_list):
        datafile_list = [datafile_list]

    df_list = []
    for datafile in datafile_list:
    
        #might need to add quality filter
        tab = Table.read(datafile, hdu=1)
        df = tab.to_pandas()
        df['hjd'] = df.TIME + 2457000


        tess_mag = fits.open(datafile)[0].header['TESSMAG']
        if ksp:
            df['mag'] = -2.5*np.log10(df.KSPSAP_FLUX)
            df['mag'] = df['mag'] + (tess_mag - df.mag.median())
            df = df[~np.isnan(df['mag'])]
            df['mag_err'] = abs((-2.5*np.log10(np.e))*(df.KSPSAP_FLUX_ERR/df.KSPSAP_FLUX))
        else:
            try:
                df['mag'] = -2.5*np.log10(df.SAP_FLUX)
                df['mag'] = df['mag'] + (tess_mag - df.mag.median())
                df = df[~np.isnan(df['mag'])]
                if 'KSPSAP_FLUX_ERR' in df.columns:
                    df['mag_err'] = abs((-2.5*np.log10(np.e))*(df.KSPSAP_FLUX_ERR/df.SAP_FLUX))
                else:
                    df['mag_err'] = 1e-4 #new QLP light curves dont have KSP column
            except:
                continue

        df['orbit'] = df.ORBITID
        df['sector'] = tessutils.parse_path(datafile).sector
        df['filter'] = 'T'
        df = df.rename(columns={'QUALITY':'quality'})
        df = df[df['quality'] == 0]
        if not df.empty:
            if np.all(df.orbit.to_numpy() == 0):
                idx = np.argmax(df.hjd.diff())
                orbit = np.full(len(df), df.sector.iloc[0]*2+7, dtype=float)
                orbit[idx:] += 0.5
                df['orbit'] = orbit
        df['mag'] = df.mag.replace([-np.inf, np.inf], np.nan)
        df = df[df.mag.notna()]
        df_list.append(df)


    if len(df_list):
        df_out = pd.concat(df_list)
        df_out = df_out.reset_index(drop=True)
    else:
        df_out = pd.DataFrame()

    return df_out

def tess_tglc(datafile_list):
    
    if not tessutils.check_iter(datafile_list):
        datafile_list = [datafile_list]

    df_list = []
    for datafile in datafile_list:
        
        lc = fits.open(datafile)
        df = Table(lc[1].data).to_pandas()
        tess_mag = lc[0].header['TESSMAG']

        df['sector'] = lc[0].header['SECTOR']
        df['filter'] = 'T'
        df['hjd'] = df.time + 2457000

        df['mag'] = -2.5*np.log10(df.cal_psf_flux)
        df['mag'] = df['mag'] + (tess_mag - df.mag.median())
        df = df[~np.isnan(df.mag)].reset_index(drop=True)
        df['mag_err'] = 0.01
        df['mag'] = df.mag.replace([-np.inf, np.inf], np.nan)
        df = df[df.mag.notna()]
        df_list.append(df)

    df_out = pd.concat(df_list)
    
    return df_out.reset_index(drop=True)

    
