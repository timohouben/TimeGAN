# -*- coding: utf-8 -*-
# ------------------------------------------------------------------------------
# Created by : Timo Houben
# Created on : December 2020
# Edited on : March 5, 2021
# ------------------------------------------------------------------------------
__author__ = "Timo Houben"
__copyright__ = "Copyright (c) 2021, Timo Houben, application"
__credits__ = ["Timo Houben"]
__license__ = "MIT"
__version__ = "0.0.0"
__maintainer__ = "Timo Houben"
__email__ = "timo.houben@gufz.de"
__status__ = "development"
# ------------------------------------------------------------------------------
# Description
'''
Class for spectral analysis of data.
'''
# ------------------------------------------------------------------------------
#
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack as fftpack
from scipy import signal
import os

class SpectralAnalysis(object):
    
    def __init__(self, csvpathlist):
        self.csvpathlist = csvpathlist
        if type(self.csvpathlist) != list:
            self.csvpathlist = list(self.csvpathlist)
        self.data = pd.read_csv(csvpathlist[0])
        # concatinate dataframes
        # missing: CHECK FOR SAME HEADER
        for i, csvpath in enumerate(self.csvpathlist[1:]):
            df_temp = pd.read_csv(csvpath)
            header = list(df_temp.columns)
            new_header = [h+"_"+str(i+1) for h in header]
            new_header = dict(zip(header, new_header))
            df_temp.rename(columns = new_header, inplace=True)
            self.data = pd.concat([self.data, df_temp], sort=False, axis=1)
        
        self.rootdir = os.path.dirname(self.csvpathlist[0])

    def get_all_power_spectra(self, time_step_size, method="scipyffthalf"):     
        self.all_power_spectra = pd.DataFrame(columns=["frequency"])
        for column in self.data:
            freq, spec = self.power_spectrum(self.data[column], time_step_size, method=method)
            df_temp = pd.DataFrame({"frequency":freq,"spectrum"+"_"+column:spec})
            self.all_power_spectra = pd.merge(self.all_power_spectra, df_temp, on="frequency", how="outer")
        return self.all_power_spectra


    def power_spectrum(self, series, time_step_size, method="scipyffthalf"):
        """
        Plots the spectra for each time series for whole period.

        This script computes the power spectral density estimate of a time series.
        You can choose between three methods. 'scipyffthalf' and 'scipyperio'
        reveals almost exactly the same results. 'scipywelch' computes a smoothed
        periodogram.
        
        Parameters
        ----------

        series : 1D array, list
            Time series of and input process of e.g. a LTI system. If considering an
            aquifer as filter of the LTI, the input signal would be equal to the
            recharge time series of the aquifer.

        time_step_size : string
            Size of the time step in seconds to obtain Hz on abszisse.
        
        method : string, Default: 'scipyffthalf'
            Method which will be used to derive the spectrum.
            'scipyffthalf'
            # ======================================================================
            # method 1: Periodogram: Power Spectral Density: abs(X(w))^2
            #           http://staff.utia.cas.cz/barunik/files/QFII/04%20-%20Seminar/04-qf.html
            # ======================================================================
            'scipywelch'
            # ======================================================================
            # method 2: scipy.signal.welch
            #           https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.signal.welch.html#r145
            # ======================================================================
            'scipyperio'
            # ======================================================================
            # method 3: Scipy.signal.periodogram
            #           https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.signal.periodogram.html
            # ======================================================================
        
        Yields
        ------
        frequency : 1D array
            Corresponding frequencies of the Fourier Transform.

        power_spectrum : 1D array
            Power spectrum of time series.

        """
        len_series = len(series)
        # define the sampling frequency/time step
        sampling_frequency = 1.0 / time_step_size  # [Hz] second: 1, day: 1.1574074074074E-5
        # calculate spectrum
        if method == "scipyffthalf":
            # calulate fast fourier transform
            fourier = fftpack.fft(series)
            # calculate the one sided power spectrum
            power_spectrum = abs(fourier[: int(round(len(fourier) / 2))]) ** 2 
            # calculate the corresponding frequencies
            frequency = (
            abs(fftpack.fftfreq(len_series, time_step_size))[
                : int(round(len_series / 2))
            ]
        )

        # smoothed spectrum by Welch et al. (????)
        elif method == "scipywelch":
            # define smoothing window size
            nperseg = int(round(len_series / 10))
            frequency, power_spectrum = signal.welch(
                series, sampling_frequency, nperseg=nperseg, window="hamming"
            )

        else: 
            print("method must be either 'scipyffthalf' or 'scipywelch'")

        return frequency, power_spectrum

    
    def plot_power_spectrum(self, frequency, spectrum):
        pass    




    def plot_all_power_spectra(self, comment, ylims=(1e-7, 1e15), xlims=(1e-8, 1e-3), withmarkers=True):
        
        year = 1/(365*24*60*60)
        month = 1/(30*24*60*60)
        week = 1/(7*24*60*60)
        day = 1/(24*60*60)

        cols = int(len(self.csvpathlist))
        rows = int((len(self.all_power_spectra.columns) - 1) / cols)

        fig, axes = plt.subplots(rows, cols, sharex=True, sharey=True, squeeze=False, figsize=(cols*4+1, rows*2+1))

        for column, ax in zip((list(self.all_power_spectra)[1:]), np.ravel(axes, order="F")):
            ax.loglog(self.all_power_spectra["frequency"], self.all_power_spectra[column], label=column)
            ax.set_title(column)
            ax.set_ylim(ylims)
            ax.set_xlim(xlims)
            if withmarkers == True:
                for vline, vlabel in zip([year, month, week, day], ["year", "month", "week", "day"]):
                    ax.vlines(vline, ymin=ylims[0], ymax=ylims[1], linestyle="--", alpha=0.3, linewidth=0.5)
                    ax.text(vline, pow(10, np.log10(ylims[0]) + 0.2), str(vlabel))

        plt.tight_layout()
        plt.savefig(os.path.join(self.rootdir, "all_power_spectra_" + comment + ".png"), dpi=300)