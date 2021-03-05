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
Script to execute spectralanalysis class on climate data from Mohit. 
'''
# ------------------------------------------------------------------------------
#
from spec import SpectralAnalysis

obs = "/Users/houben/phd/studies/timegan/data/cllimate_data/clim_obs.csv"
sim1 = "/Users/houben/phd/studies/timegan/data/cllimate_data/clim_sim_1.csv"
sim2 = "/Users/houben/phd/studies/timegan/data/cllimate_data/clim_sim_2.csv"

spectral = SpectralAnalysis([obs, sim1, sim2])
spectral.data

all_power_spectra = spectral.get_all_power_spectra(time_step_size=3600, method="scipyffthalf")
spectral.plot_all_power_spectra(comment="scipyffthalf")

all_power_spectra = spectral.get_all_power_spectra(time_step_size=3600, method="scipywelch")
spectral.plot_all_power_spectra(comment="scipywelch")
