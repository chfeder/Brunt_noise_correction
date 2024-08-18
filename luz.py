#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Christoph Federrath 2024

import numpy as np
from cfpack import stop, print, matplotlibrc
import cfpack as cfp


# fit function
def fitfunc(x, a, b):
    return np.log( x**a / (x**a + b) )


# ===== the following applies in case we are running this in script mode =====
if __name__ == "__main__":

    # read data
    dat = cfp.read_ascii("luz.txt")
    x = dat["SNR"]
    y = dat["sigma_ratio"]

    # plot data
    cfp.plot(x=x, y=y, marker='o', linestyle="", label="data")

    # fit data
    params = {"a": [0.1, 1, 3.0], "b": [0.1, 1, 3.0]}
    fitres = cfp.fit(fitfunc, x, np.log(y), params=params)

    # plot fit
    xfit = np.logspace(-2, 2, 500)
    yfit = np.exp(fitfunc(xfit, *fitres.popt))
    cfp.plot(x=xfit, y=yfit, label="fit")

    # final plot command to save formatted figure
    cfp.plot(xlabel=r"SNR", ylabel=r"$\sigma_\mathrm{\rho,\,\infty}\,/\,\sigma_\mathrm{\rho,\,SNR}$",
             xlog=True, ylog=True, xlim=[5e-3,2e2], ylim=[1e-3,2], save="luz.pdf")

