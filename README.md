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
