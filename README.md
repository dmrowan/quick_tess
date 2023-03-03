# quick_tess
Download and inspect TESS light curves

Set the QUICKTESS_DIR environmental variable to store where LCs are downloaded to.

Example usage:
```
from quick_tess import mast_download, inspector_plot, data_parser

path = mast_download.download(302105449) #enter TIC
inspector_plot.plot(path, parser=data_parser.tess_qlp,
                    savefig='temp.pdf', samples_per_peak=20)
```

The default pipeline is the Quick Look Pipeline QLP. You can search for and download SPOC data with the `pipeline` argument:

```
from quick_tess import mast_download, inspector_plot, data_parser

path = mast_download.download(302105449, pipeline='TESS-SPOC') #enter TIC
inspector_plot.plot(path, parser=data_parser.tess_spoc, #notice the parser argument changed
                    savefig='temp.pdf', samples_per_peak=20)
```

Not all TICs have QLP or SPOC light curves, even if they have been observed. We can use the TGLC tool to get LCs from the FFIs. Note that this is fairly slow:

```
from quick_tess import tglc_download, inspector_plot, data_parser
path = tglc_download.download(302105449)
inspector_plot.plot(path, parser=data_parser.tess_tglc)
```

You can also pass any of the optional arguments for tglc.quick_lc.tglc_lc directly to tglc_download.download
