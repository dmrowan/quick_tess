#!/usr/bin/env python

from astropy import log
from astropy.io import fits
from astroquery.mast import Observations
from astropy import units
import numpy as np
import os

desc="""
Download tess data from MAST for QLP and SPOC pipelines
"""

from . import tessutils

data_path = os.environ.get('QUICKTESS_DIR', None)

if data_path is None:
    print("Please set the QUICKTESS_DIR environment variable to path ot the downloaded data")


def download(TIC=None, ra=None, dec=None,
             pipeline='QLP',
             download_dir=data_path,
             verbose=False):

    full_path = None
    if TIC is None:
        assert(ra is not None)
        assert(dec is not None)

        obs_table = Observations.query_region(
                f'{ra} {dec}',
                radius=5.0*units.arcsecond, verbose=verbose)
        idx = np.where(obs_table['provenance_name'] == pipeline)[0]
        full_path = []
        if len(idx) != 0:
            TIC_arr = np.array(obs_table['target_name'][idx])
            if len(np.unique(TIC_arr)) != 1:
                idx = [i for i in idx 
                       if obs_table['target_name'][i]==obs_table['target_name'][idx[0]] ]
            
            TIC = obs_table['target_name'][idx[0]]
            for i in range(len(idx)):
                print(obs_table['obs_id'][idx[i]])
                full_path.append(
                        os.path.join(
                                asassn_config.tessdir.tess_cache, 
                                'mastDownload/HLSP', 
                                obs_table['obs_id'][idx[i]],
                                obs_table['obs_id'][idx[i]]+'.fits'))
            if pipeline == 'TESS-SPOC':
                full_path = [ x.replace('_tp.fits', '_lc.fits') 
                              for x in full_path ]
        else:
            log.info(f'No TESS data on MAST for {ra} {dec}')
            return None

    obsTable = Observations.query_criteria(provenance_name=pipeline,
                                           target_name=TIC, verbose=verbose)
        
    if len(obsTable) == 0:
        return None

    data = Observations.get_product_list(obsTable, verbose=verbose)
    with tessutils.HiddenPrints():
        download_lc = Observations.download_products(
                data, download_dir=download_dir)
        
    if full_path is None:
        if pipeline == 'QLP':
            full_path = [ os.path.join(download_dir, 'mastDownload/HLSP',
                                       x, x+'.fits') 
                          for x in obsTable['obs_id'] ]
        elif pipeline == 'TESS-SPOC':
            full_path = [ os.path.join(download_dir, 'mastDownload/HLSP',
                                       x, x.replace('_tp', '_lc.fits'))
                          for x in obsTable['obs_id'] ]
        else:
            raise ValueError(
                    f'pipeline name {pipeline} not recognized')

    return full_path



