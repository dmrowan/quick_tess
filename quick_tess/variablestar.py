#/usr/bin/env python

from astropy import units as u
from astropy.timeseries import LombScargle
from astropy import log
from brokenaxes import brokenaxes
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import pandas as pd
import pickle
import random
from scipy import ndimage
from scipy import signal
from scipy import stats

from . import data_parser
from . import tessutils
from . import plotutils

#Dom Rowan 2023

desc="""
Basic variable star class for quick_tess
"""

colors = ['#DE1509', '#C9189E', '#AA2AFA', '#3524D1', '#1D73F2']

class VariableStar:

    def __init__(self, datafile,
                 parser=None, 
                 parser_args=None, 
                 verbose=True):

        if parser_args is None:
            parser_args = {}

        if parser is not None:
            self.df = parser(datafile, **parser_args)
        else:
            self.df = datafile

        if type(self.df) != pd.core.frame.DataFrame:
            raise TypeError('parser must return pd.core.frame.DataFrame')

        #Initialize params
        self.datafile = datafile
        self.verbose = verbose
        self._sigma_clipped = False

        try:
            self.timecol = 'hjd'
        except:
            self.timecol = self.df.columns[0]

        self._t0 = None
        self.df_clipped = pd.DataFrame(columns=self.df.columns)
        self.df_clipped['clipped_by'] = None
        self.scatter_color = 'black'

    @property
    def verbose(self):
        return self._verbose

    @verbose.setter
    def verbose(self, value):
        self._verbose=bool(value)

    @property
    def timecol(self):
        return self._timecol

    @timecol.setter
    def timecol(self, value):
    
        if len(self.df.columns) == 0:
            self._timecol = value
        elif value in self.df.columns:
            self._timecol = value
        else:
            raise ValueError('timecol not found in VariableStar.df')

    @property
    def t0(self):
        return self._t0

    @t0.setter
    def t0(self, value):
        if value is None:
            pass
        elif value < 0:
            raise ValueError('t0 must be postiive')
        self._t0 = value

    @property
    def sigma_clipped(self):
        return self._sigma_clipped

    @sigma_clipped.setter
    def sigma_clipped(self, value):
        self._sigma_clipped = bool(value)

    def ls_periodogram(self, maximum_frequency=10, minimum_frequency=None,
                       yvals='mag', samples_per_peak=5,
                       plot=False, ax=None, savefig=None,
                       plot_kwargs=None):

        #Create astropy LS periodogram

        ls = LombScargle(self.df[self.timecol], self.df[yvals])
        fal = ls.false_alarm_level(1e-5)

        #Use autofrequeny grid over specified range
        self.ls_freq, self.ls_power = ls.autopower(
                minimum_frequency=minimum_frequency,
                maximum_frequency=maximum_frequency,
                samples_per_peak=samples_per_peak)
        self.ls = ls


        df_ls = pd.DataFrame({'freq':self.ls_freq,
                              'power':self.ls_power})
        df_ls = df_ls.sort_values(by='freq', ascending=True, 
                                  ignore_index=True)

        self.df_ls = df_ls

        #Select the period corresponding to the highest peak in the periodogram
        self.period = np.power(df_ls.freq.iloc[np.argmax(df_ls.power)], -1)

        if (plot) or (savefig is not None) or (ax is not None):
            
            fig, ax, created_fig = plotutils.fig_init(ax=ax)

            if plot_kwargs is None:
                plot_kwargs = {}
            plot_kwargs.setdefault('color', 'black')
            plot_kwargs.setdefault('zorder', 1)

            ax = plotutils.plotparams(ax)
            ax.set_xlabel('Frequency (1/d)', fontsize=20)
            ax.set_ylabel('Power', fontsize=20)
            ax.axhline(fal, color='gray', label=r'FAL=${}$'.format(1e-5),
                       lw=3, zorder=2)

            ax.axvline(1/self.period, color='black', lw=1.5, alpha=0.4)

            ax.plot(self.ls_freq, self.ls_power, **plot_kwargs)

            return plotutils.plt_return(created_fig, fig, ax, savefig)

        return self.period

    def pdm(self, prange=0.1, period=None, niters=10000,
            plot=False, ax=None, savefig=None, progress=True):
        
        if period is None:
            period = self.period

        period_arr = np.linspace(period-prange, period+prange, niters)

        variances = tessutils.pdm(self.df[self.timecol].to_numpy(), 
                                  self.df.mag.to_numpy(),
                                  period_arr, progress=progress)

        self.period = period_arr[np.argmin(variances)]

        if (plot) or (ax is not None) or (savefig is not None):
            
            fig, ax, created_fig = plotutils.fig_init(ax=ax, figsize=(8,4))
            if created_fig:
                fig.subplots_adjust(top=.98, right=.98)
            ax.plot(period_arr, variances, color='black', ls='-', lw=2, zorder=2)
            ax.set_xlabel('Period (d)', fontsize=20)
            ax.set_ylabel(r'$\sigma_m^2$ (mag)', fontsize=20)
            ax.axvline(self.period, ls='--', lw=1, color='xkcd:red', zorder=1)

            return plotutils.plt_return(created_fig, fig, ax, savefig)

    @property
    def is_phase_folded(self):
        return 'phase' in self.df.columns

    #Phase fold data
    def phase_fold(self, t0=None, n=1):
        
        if not hasattr(self, 'period'):
            if self._verbose:
                log.info('Computing period with LS periodogram')
            self.ls_periodogram()

        self.period_n = n
        self.period = self.period*self.period_n

        if len(self.df) != 0:
            if t0 is not None:
                self.t0 = t0

            if self.t0 is None:
                self.t0 = min(self.df[self.timecol])

            self.df['phase'] = ((self.df[self.timecol]-self.t0)%self.period)/self.period
            self.df['phi'] = self.df.phase*2*np.pi
        else:
            self.df['phase'] = []
            self.df['phi'] = []

    def _fold_clipped(self):
        
        t0 = self.t0
        period = self.period

        if (('phase' not in self.df_clipped.columns) or
            (np.any(self.df_clipped.phase.isna()))):

            if not self.df_clipped.empty:
                self.df_clipped['phase'] = ((self.df_clipped[self.timecol]-t0)%period)/period
                self.df_clipped['phi'] = self.df_clipped.phase*2*np.pi
            else:
                self.df_clipped['phase'] = None
                self.df_clipped['phi'] = None

    def plot_full(self, ax=None, savefig=None, 
                  time_offset=0,color_column=None, color_dict=None, 
                  legend=False, plot_kwargs=None,
                  yvals='mag'):

        if plot_kwargs is None:
            plot_kwargs = {}

        plot_kwargs.setdefault('marker', '.')
        plot_kwargs.setdefault('edgecolor', 'none')
        plot_kwargs.setdefault('alpha',0.6)

        if yvals not in list(self.df.columns):
            raise ValueError(f'column {yvals} not found in vs')

        #Create figure
        fig, ax, created_fig = plotutils.fig_init(ax=ax)

        if time_offset == 0:
            ax.set_xlabel(plotutils.format_latex_label(self.timecol.upper()), 
                          fontsize=20)
        else:
            ax.set_xlabel(
                    plotutils.format_latex_label(
                            self.timecol.upper())+r'$-{}$'.format(time_offset),
                    fontsize=20)

        ax.set_ylabel(plotutils.format_latex_label(yvals), fontsize=20)

        if color_column is None:
            if 'color' not in list(plot_kwargs.keys()):
                plot_kwargs['color'] = self.scatter_color

            ax.scatter(self.df[self.timecol]-time_offset, 
                       self.df[yvals],
                       **plot_kwargs)
        else:
            if color_dict is None:
                color_dict = dict(zip(list(set(self.df[color_column])), 
                                           plotutils.many_colors()))

            ax.scatter(self.df[self.timecol]-time_offset, self.df[yvals],
                       color=[color_dict[c] for c in self.df[color_column]],
                       **plot_kwargs)

        if not ax.yaxis_inverted() and 'mag' in yvals:
            ax.invert_yaxis()

        if legend:
            current_xlim = ax.get_xlim()
            current_ylim = ax.get_ylim()
            for value in tessutils.sort_column(list(
                    set(self.df[color_column]))):
                ax.scatter([-1], [0], color=color_dict[value], 
                           alpha=0.6, label=value)
            ax.set_xlim(current_xlim)
            ax.set_ylim(current_ylim)
            ax.legend(loc='upper right', edgecolor='black', fontsize=15)

        return plotutils.plt_return(created_fig, fig, ax, savefig)

    def plot_broken(self, fig=None, gs=None, savefig=None, xlims=None,
                    combine_adjacent_sectors=False,
                    plot_kwargs=None, yvals='mag', time_offset=0):

        if plot_kwargs is None:
            plot_kwargs = {}

        plot_kwargs.setdefault('color', self.scatter_color)
        plot_kwargs.setdefault('marker', '.')
        plot_kwargs.setdefault('edgecolor', 'none')

        #Determine breaks
        if xlims is None:
            if 'sector' in self.df.columns:
                xlims = [ 
                        (self.df[self.df.sector==sector][
                                self.timecol].min()-time_offset-1,
                         self.df[self.df.sector==sector][
                                self.timecol].max()-time_offset+1)
                        for sector in np.sort(
                                np.unique(self.df.sector.to_numpy())) ]
            else:
                xlims=None

        if combine_adjacent_sectors:
            sector_list = np.sort(np.unique(self.df.sector.to_numpy()))

            runs = tessutils.find_runs(sector_list)

            xlims_new = []
            for i in range(len(xlims)):
                
                if i not in [ item for sublist in runs for item in sublist ]:
                    xlims_new.append(xlims[i])
                else:
                    run = [r for r in runs if i in r][0]
                    if i == run[0]:
                        xlims_new.append( (xlims[i][0], xlims[run[-1]][1]) )
                    else:
                        continue
            
            xlims = xlims_new

        
        if fig is None:
            created_fig = True
            fig = plt.Figure(figsize=(12, 6))
            bax = brokenaxes(fig=fig, xlims=xlims)
        else:
            created_fig = False
            bax = brokenaxes(fig=fig, subplot_spec=gs, xlims=xlims)

        bax.scatter(self.df[self.timecol]-time_offset, 
                    self.df[yvals], **plot_kwargs)

        if not bax.yaxis_inverted()[0]:
            bax.invert_yaxis()

        if time_offset == 0:
            bax.big_ax.set_xlabel(
                    plotutils.format_latex_label(self.timecol.upper()), 
                    fontsize=20, labelpad=25)
        else:
            bax.big_ax.set_xlabel(
                    plotutils.format_latex_label(
                            self.timecol.upper())+r'$-{}$'.format(time_offset),
                    fontsize=20, labelpad=25)
        bax.axs[0].set_ylabel(plotutils.format_latex_label(yvals), fontsize=20)

        bax = plotutils.plotparams_bax(bax)

        if created_fig:
            if savefig is not None:
                fig.savefig(savefig)
            else:
                plt.show()
        else:
            return bax

    def plot_lc(self, nphase=2, ax=None, savefig=None, 
                phase_column='phase', plot_binned=True,
                color_column=None, color_dict=None,
                legend=False, title=False, plot_kwargs=None,
                label_period=True, plot_clipped=True,
                errorbar=False, yvals='mag', yerr='mag_err'):

        if plot_kwargs is None:
            plot_kwargs = {}

        plot_kwargs.setdefault('marker', '.')
        plot_kwargs.setdefault('alpha', 0.6)
        if errorbar:
            plot_kwargs.setdefault('ls', '')
            plot_kwargs.setdefault('markeredgecolor', 'none')
        else:
            plot_kwargs.setdefault('edgecolor', 'none')
        if 'phase' not in self.df.columns:
            self.phase_fold()

        if phase_column not in self.df.columns:
            raise ValueError(f'column {phase_column} not found in vs')

        if yvals not in list(self.df.columns):
            raise ValueError(f'column {yvals} not found in vs')

        if errorbar and (yerr not in list(self.df.columns)):
            raise ValueError(f'column {yerr} not found in vs')

        #Create figure
        fig, ax, created_fig = plotutils.fig_init(ax=ax)

        if created_fig:
            fig.subplots_adjust(top=.98, right=.98)

        if color_column is None:
            if errorbar:
                if 'c' not in list(plot_kwargs.keys()):
                    plot_kwargs.setdefault('color', self.scatter_color)
                ax.errorbar(
                        *tessutils.extend_phase(self.df[phase_column], self.df[yvals], nphase),
                        yerr=tessutils.extend_phase(
                                self.df[phase_column], self.df[yerr], nphase)[1],
                         **plot_kwargs)
            else:
                if 'c' not in list(plot_kwargs.keys()):
                    plot_kwargs.setdefault('color', self.scatter_color)
                ax.scatter(
                        *tessutils.extend_phase(
                                self.df[phase_column], 
                                self.df[yvals], nphase),
                        **plot_kwargs)
        else:
            if color_dict is None:
                colors = plotutils.many_colors()
                random.shuffle(colors)
                color_dict = dict(zip(list(set(self.df[color_column])), colors))

            plot_kwargs.pop('color', None)
            ax.scatter(*tessutils.extend_phase(self.df[phase_column], self.df[yvals], nphase),
                       color=[ color_dict[c] for c in np.tile(self.df[color_column], nphase) ],
                       **plot_kwargs)
                       #**{i:d[i] for i in d if i not in ['ecolor', 'ls']})
            if errorbar:
                ax.errorbar(
                        *tessutils.extend_phase(self.df[phase_column], self.df[yvals], nphase),
                        yerr=tessutils.extend_phase(
                                self.df[phase_column], self.df[yerr], nphase)[1],
                        marker=None,
                        color=[ color_dict[c] for c in np.tile(self.df[color_column], nphase) ],
                        **plot_kwargs)

        #Plot binned LC (from Tharindu)
        if plot_binned:
            bins_ph=np.arange(-0.05,1.05,0.05) #bin the light curve in phase
            bin_means, bin_edges, binnumber = stats.binned_statistic(
                    self.df.phase,self.df[yvals], statistic='median', bins=bins_ph) 
            bin_std, bin_edges2, binnumber2 = stats.binned_statistic(
                    self.df.phase,self.df[yvals], statistic='std', bins=bins_ph)
            bin_width = (bin_edges[1] - bin_edges[0])
            bin_centers = bin_edges[1:] - bin_width/2
            lcbin=pd.DataFrame({'phase':bin_centers,'yvals':bin_means,'errs':bin_std})

            ax.errorbar(*tessutils.extend_phase(lcbin.phase, lcbin.yvals, nphase),
                        yerr=tessutils.extend_phase(lcbin.phase, lcbin.errs, nphase)[1],
                        marker='o', color='red', ls='')

        #Plot clipped points
        if not self.df_clipped.empty and plot_clipped:

            if errorbar:
                ax.errorbar(
                        *tessutils.extend_phase(
                                self.df_clipped[phase_column],
                                self.df_clipped[yvals], nphase),
                        yerr=tessutils.extend_phase(
                                self.df_clipped[phase_column], self.df_clipped[yerr], nphase)[1],
                        color='gray', label='clipped', alpha=0.8, ls='', 
                        marker=plot_kwargs['marker'])
            else:
                ax.scatter(*tessutils.extend_phase(
                        self.df_clipped[phase_column], self.df_clipped[yvals],nphase),
                           color='gray', alpha=0.8, marker=plot_kwargs['marker'])

        ax.set_xlabel('Phase', fontsize=20)

        if (not self.df.empty) and ('filter' in self.df.columns):
            filt = tessutils.filter_format(self.df['filter'].value_counts().index[0])
            ax.set_ylabel(f'{filt} Mag', fontsize=20)
        else:
            ax.set_ylabel('Mag', fontsize=20)

        if not ax.yaxis_inverted() and 'flux' not in yvals:
            ax.invert_yaxis()

        if legend:
            if color_column != 'camera':
                current_xlim = ax.get_xlim()
                current_ylim = ax.get_ylim()
                for value in tessutils.sort_column(list(set(self.df[color_column]))):
                    ax.scatter([-1], [0], color=color_dict[value], alpha=0.6, label=value)
                ax.set_xlim(current_xlim)
                ax.set_ylim(current_ylim)
            ax.legend(loc='upper right', edgecolor='black', fontsize=15)

        elif label_period:
            current_ylim = ax.get_ylim()
            ylim_range = current_ylim[0]-current_ylim[1]
            ax.set_ylim(current_ylim[0], current_ylim[1]-0.08*ylim_range)
            ax.text(.05, .95, r'$P={}$d'.format(round(self.period, 2)),
                    fontsize=18, transform=ax.transAxes, 
                    ha='left', va='top',
                    bbox=dict(facecolor='white', edgecolor='none',
                              alpha=0.5))

        return plotutils.plt_return(created_fig, fig, ax, savefig)

    #Sigma clip using std and median
    def sigma_clip(self, nsigma=5):

        #Compute limits
        std = self.df.mag.std()
        
        if not tessutils.check_iter(nsigma):
            nsigma = [nsigma, nsigma]

        #order of nsigma is (faint, brighter)
        lower = self.df.mag.median()-nsigma[1]*std
        upper = self.df.mag.median()+nsigma[0]*std
        to_clip = np.array([ not (lower <= m <= upper) for m in self.df.mag])

        if any(to_clip):
            idx_clip = np.where(to_clip)[0]
            idx_keep = np.where(~to_clip)[0]

            if self._verbose:
                log.info(f"{len(idx_clip)} points sigma-clipped")

            df_clipped = self.df.iloc[idx_clip,].copy()
            df_clipped['clipped_by'] = 'sigma_clip'
            self.df_clipped = pd.concat([self.df_clipped, df_clipped]).reset_index(
                    drop=True)
            self.df = self.df.iloc[idx_keep,].reset_index(drop=True)
        else:
            self.df_clipped=pd.DataFrame(columns=self.df.columns)

        #update sigma clipping flag
        self._sigma_clipped=True
        self.nsigma=nsigma

    @property
    def clipped_fraction(self):
        
        return len(self.df_clipped)/(len(self.df_clipped)+len(self.df))

    #Clip using median filter
    def median_filter(self, nsigma=5, window=9, rolling_cutoff=0.4,
                      use_rolling=False,
                      plot=False, savefig=None, ax=None,
                      wide_ylim=True, no_clip=False):

        if not (0 < rolling_cutoff <= 1):
            raise ValueError("Rolling cutoff must be between 0 and 1")

        if isinstance(nsigma, str):
            if nsigma not in self.df.columns:
                raise ValueError(f"nsigma column {nsigma} not found in vs.df")
            nsigma = self.df[nsigma]

        #Sigma clipping can only be performed once
        if self._sigma_clipped:
            if self._verbose:
                log.info("Sigma clipping already performed")
            return self.df

        #Fold LC first
        if 'phase' not in self.df.columns:
            self.phase_fold()

        if all([np.isnan(e) for e in self.df.mag_err]):
            use_rolling=True

        #Sort df by phase
        df_sorted = self.df.sort_values(by='phase')
        df_sorted = df_sorted.reset_index(drop=True)

        #Run median filter with wrapping
        df_sorted['smooth'] = ndimage.median_filter(
                df_sorted.mag, window, mode='wrap')

        df_sorted['mf_upper'] = df_sorted.smooth+nsigma*df_sorted.mag_err
        df_sorted['mf_lower'] = df_sorted.smooth-nsigma*df_sorted.mag_err

        if no_clip:
            self.df = df_sorted
            return self.df['smooth']

        to_clip = (np.abs(df_sorted.mag - df_sorted.smooth) > 
                   nsigma*df_sorted.mag_err)

        #See if points meet clipping criteria
        if any(to_clip) or use_rolling:
            #Identify relevant indicies
            idx_clip = np.where(to_clip)[0]
            idx_keep = np.where(~to_clip)[0]

            #Alternative clip for high amplitude pulsators (ex: M)
            if (len(idx_clip)/len(df_sorted) >= rolling_cutoff) or use_rolling:

                #Take the median filter array
                padded = np.concatenate(
                        (np.array(df_sorted.smooth.iloc[-window//2+1:]),
                         np.array(df_sorted.smooth),
                         np.array(df_sorted.smooth.iloc[:window//2])))
                ts = pd.Series(padded)
                rolling_std_padded = ts.rolling(window=window, center=True).std()
                rolling_std = np.array(rolling_std_padded[~np.isnan(rolling_std_padded)])

                df_sorted['rolling_std'] = rolling_std

                to_clip = (np.abs(df_sorted.mag - df_sorted.smooth) >
                           nsigma*df_sorted.rolling_std)
                        
                idx_clip = np.where(to_clip)[0]
                idx_keep = np.where(~to_clip)[0]

                self._used_rolling_std = True

                #Define columns for data frame
                df_sorted['mf_upper'] = df_sorted.smooth+nsigma*df_sorted.rolling_std
                df_sorted['mf_lower'] = df_sorted.smooth-nsigma*df_sorted.rolling_std
            else:
                self._used_rolling_std = False

            if self._verbose:
                log.info(f"{len(idx_clip)} points clipped using median filter")

            #Applying clipping to data
            df_clipped = df_sorted.iloc[idx_clip,].copy()
            df_clipped['clipped_by'] = 'mf'
            self.df_clipped = pd.concat([self.df_clipped, df_clipped]).reset_index(drop=True)
            self.df = df_sorted.iloc[idx_keep,].reset_index(drop=True)

        else: #creat empty dataframe
            self.df = df_sorted.reset_index(drop=True)
            self._used_rolling_std=False

        #Update sigma clipping parameters
        self.nsigma = nsigma
        self.median_window = window
        self._sigma_clipped=True

        #Plotting options
        if (plot) or (savefig is not None) or (ax is not None):
            
            fig, ax, created_fig = plotutils.fig_init(ax=ax)

            #Plot data, median lc, and clipped points
            ax.plot(self.df.phase, self.df.smooth, color=colors[4],
                    lw=4, zorder=2, label='Median LC')
            ax.scatter(self.df.phase, self.df.mag, color=self.scatter_color,
                       zorder=1, label='_nolegend_')
            ylim = ax.get_ylim()
            if not self.df_clipped.empty:
                ax.scatter(self.df_clipped.phase, self.df_clipped.mag, 
                           color='gray', zorder=2,
                           label='clipped')
            if not wide_ylim:
                ax.set_ylim(ylim)

            #Shade region where points are kept
            fill_color = colors[0]
            ax.fill_between(self.df.phase, self.df.smooth, self.df.mf_upper,
                            color=fill_color, alpha=0.2,
                            zorder=0.5, label='_nolegend_')
            ax.fill_between(self.df.phase, self.df.smooth, self.df.mf_lower,
                            color=fill_color, alpha=0.2,
                            zorder=0.5, label='_nolegend_')

            #Misc plotting params
            ax = plotutils.plotparams(ax)
            ax.set_ylim(ax.get_ylim()[::-1])
            ax.set_xlabel("Phase", fontsize=20)
            ax.set_ylabel("Mag", fontsize=20)
            ax.legend(fontsize=15, edgecolor='black', loc='upper right')

            return plotutils.plt_return(created_fig, fig, ax, savefig)

    #manual outlier rejection
    def manual_clip(self, npoints=1, clip_bright=False, clip_faint=False, timeout=30):

        fig, ax, created_fig = plotutils.fig_init()

        ax.scatter(self.df.phase, self.df.mag, color=self.scatter_color)
        ax.set_xlabel('Phase', fontsize=20)
        ax.set_ylabel('Mag', fontsize=20)
        ax.set_title("Clipping Points Manually", color='red', fontsize=15)
        ax.set_ylim(ax.get_ylim()[::-1])

        selected = plt.ginput(npoints, timeout=timeout)

        plt.close()

        matched_indicies = []
        for s in selected:

            distances = np.sqrt(np.power(s[0]-self.df.phase, 2)+np.power(s[1]-self.df.mag,2))
            matched_idx = np.where(distances == min(distances))[0][0]
            matched_indicies.append(matched_idx)

        idx_keep = [ i for i in range(len(self.df)) if i not in matched_indicies ]

        df_clipped = self.df.iloc[matched_indicies].copy()
        df_clipped['clipped_by'] = 'manual_clip'
        self.df_clipped = pd.concat([self.df_clipped, df_clipped]).reset_index(drop=True)

        self.df = self.df.iloc[idx_keep,].reset_index(drop=True)

        self._sigma_clipped=True

    def _apply_smooth(self, smooth_func, yvals, column_tag):

        new_column = f'{column_tag}_{yvals}'
        if new_column in self.df.columns:
            log.warning(f'Overwriting column {new_column}')

        df_sorted = self.df.sort_values(by='phase', ascending=True)
        df_sorted[new_column] = smooth_func(df_sorted[yvals])
        self.df = df_sorted.sort_index()

    def _plot_smooth(self, smooth_column, ax=None, savefig=None, 
                     plot_kwargs=None, plot_lc_kwargs=None):

        if plot_kwargs is None:
            plot_kwargs = {}

        if plot_lc_kwargs is None:
            plot_lc_kwargs = {}

        fig, ax, created_fig = plotutils.fig_init(ax=ax)

        plot_kwargs.setdefault('color', 'xkcd:green')
        plot_kwargs.setdefault('lw', 3)
        plot_kwargs.setdefault('ls', '-')
        plot_kwargs.setdefault('alpha', 0.7)

        ax = self.plot_lc(nphase=1, ax=ax, plot_kwargs=plot_lc_kwargs)
        df_plot = self.df.sort_values(by='phase', ascending=True)
        ax.plot(df_plot.phase, df_plot[smooth_column], **plot_kwargs)

        return plotutils.plt_return(created_fig, fig, ax, savefig)

    def savgol_filter(self, window, order, mode='wrap', yvals='mag',
                      plot=False, ax=None, savefig=None, 
                      plot_lc_kwargs=None,
                      plot_kwargs=None):

        if plot_lc_kwargs is None:
            plot_lc_kwargs = {}

        if plot_kwargs is None:
            plot_kwargs = {}

        #Define smooth func
        smooth_func = lambda x: signal.savgol_filter(
                x, window, order, mode=mode)

        #Apply to axis
        self._apply_smooth(smooth_func, yvals, 'savgol')

        #Plot results
        if (plot) or (ax is not None) or (savefig is not None):
            plot_kwargs.setdefault('label', 'Savitzky-Golay Filter')

        if (plot) or (ax is not None) or (savefig is not None):

            plot_kwargs.setdefault('label', 'Savitzky-Golay Filter')
            return self._plot_smooth(
                    f'savgol_{yvals}', ax=ax, savefig=savefig,
                    plot_kwargs=plot_kwargs, plot_lc_kwargs=plot_lc_kwargs)

    def gauss_smooth(self, window, order=0, mode='wrap', yvals='mag',
                        plot=False, ax=None, savefig=None,
                        plot_lc_kwargs=None,
                        plot_kwargs=None):

        if plot_lc_kwargs is None:
            plot_lc_kwargs = {}

        if plot_kwargs is None:
            plot_kwargs = {}

        #Define smooth function
        smooth_func = lambda x: ndimage.gaussian_filter(x, window, order=order, mode=mode)

        #Apply to axis
        self._apply_smooth(smooth_func, yvals, 'gauss')

        #Plot results
        if (plot) or (ax is not None) or (savefig is not None):
            plot_kwargs.setdefault('label', 'Gaussian Filter')
            return self._plot_smooth(
                    f'gauss_{yvals}', ax=ax, savefig=savefig, 
                    plot_kwargs=plot_kwargs, plot_lc_kwargs=plot_lc_kwargs)
        
    #Fit a polynomial to unfolded LC
    def fit_poly(self, n=1, plot=False, savefig=None, ax=None,
                 plot_full_kwargs=None, show_median=False, 
                 color=colors[0]):

        if plot_full_kwargs is None:
            plot_full_kwargs = {}

        if 'time_offset' not in list(plot_full_kwargs.keys()):
            time_offset=0
        else:
            time_offset = plot_full_kwargs['time_offset']

        #Fit a polynomial of degree N
        if not tessutils.check_iter(n):
            n = [n]

        self.pf = [np.polyfit(self.df[self.timecol]-time_offset, self.df.mag, x) for x in n]
        #Compute the chi2
        self.poly_chi2 = [ tessutils.compute_chi2(
                self.df.mag, np.poly1d(x)(self.df[self.timecol]-time_offset), 
                self.df.mag_err) for x in self.pf ]

        if len(self.pf) == 1:
            self.pf = [self.pf[0]]

        if len(self.poly_chi2) == 1:
            self.poly_chi2 = self.poly_chi2[0]

        #Plot the fit over unfolded light curve
        if (plot) or (ax is not None) or (savefig is not None):
            fig, ax, created_fig = plotutils.fig_init(ax=ax)
            ax = self.plot_full(ax=ax, **plot_full_kwargs)

            xvals = np.linspace(self.df.hjd.min()-time_offset, self.df.hjd.max()-time_offset)

            yvals = [np.poly1d(x)(xvals) for x in self.pf]
            if not tessutils.check_iter(color):
                color = [color]
            if len(color) < len(yvals):
                raise ValueError('pass colors as list')

            for i in range(len(yvals)):
                ax.plot(xvals, yvals[i], color=color[i], lw=3)

            if show_median:
                ax.plot(xvals, [self.df.mag.mean()]*len(xvals), lw=2, color='gray', ls='--')

            return plotutils.plt_return(created_fig, fig, ax, savefig)
        else:
            return self.pf, self.poly_chi2

    def phase_mask(self, phase_range, phase_column='phase'):
        
        if phase_column not in list(self.df.columns):
            raise ValueError(f'phase column {phase_column} not found')

        if len(np.array(phase_range).shape) == 1:
            phase_range = [phase_range]

        for i in range(len(phase_range)):
            
            idx_mask = np.where( (self.df[phase_column] > phase_range[i][0]) &
                                 (self.df[phase_column] < phase_range[i][1]) )[0]
            if len(idx_mask):
                self.df = self.df[~self.df.index.isin(idx_mask)]
                self.df = self.df.reset_index(drop=True)

        return self.df

    def to_period04_format(self, outfile):
        
        df_out = self.df[[self.timecol, 'mag']].copy()
        if outfile is not None:
            df_out.to_csv(outfile, index=False, sep='\t', header=False)

        return df_out

    def output_full_lightcurve(self, outfile):
        
        if outfile is not None:
            self.df.to_csv(outfile, index=False)

        return self.df
		
    #Save to pickle
    def to_pickle(self, out):
        with open(out, 'wb') as p:
            pickle.dump(self, p)

    def __repr__(self):
        outs = f"VariableStar: {self.name}\n"
        outs += self.df.head().to_string() 
        outs += f"\nN obs:{len(self.df)}"
        return outs
