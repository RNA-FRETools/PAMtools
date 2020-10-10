#!/usr/bin/env python3

import numpy as np
import json
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.optimize
import scipy.special
import pandas as pd
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
import naturalcolors.colorpalette as ncp

# seaborn settings
sns.set_style('white')
sns.set_context("notebook")
sns.set(font='Arial')

#cols = [[0.31, 0.45, 0.56], [0.6, 0.6, 0.6], [0.75, 0.51, 0.38]]
cols = [[0.1, 0.1, 0.1], [0.9, 0.9, 0.9]]
cmap_list, cmap_linseg = ncp.make_colormap(cols, 'grayscale')
#colors = ncp.get_colors(cmap_linseg, 5, scramble=True)
colors = ncp.get_colors(cmap_linseg, 5)
bwo = ncp.get_cmap('bluewhiteorange')


def set_ticksStyle(x_size=4, y_size=4, x_dir='in', y_dir='in'):
    sns.set_style('ticks', {'xtick.major.size': x_size, 'ytick.major.size': y_size, 'xtick.direction': x_dir, 'ytick.direction': y_dir})


def fit(fun, x_data, y_data, p0, bounds=(-np.inf, np.inf)):
    """
    Wrapper for the curve_fit function of the scipy.optimize module.
    The curve_fit optimizes the mean and sigma parameters
    while the nnls optimizes the relative weights of the Gaussians.

    Parameters
    ----------
    fun : callable
          The model function f(x,...) taking x values as a first argument followed by the function parameters
    x_data : array_like
             array of the independent variable
    y_data : array_like
             array of the dependent variable
    p0 : array_like
         start values for the fit model
    bounds : 2-tuple of float or 2-tuple of array_like, optional
             lower and upper bounds for each parameter in p0. Can be either a tuple of two scalars
             (same bound for all parameters) or a tuple of array_like with the same length as p0.
             To deactivate parameter bounds set: `bounds=(-np.inf, np.inf)`

    Returns
    -------
    p : ndarray
        optimized fit parameters
    p_std : ndarray
            standard deviation of optimized fit parameters
    """
    p, cov = scipy.optimize.curve_fit(fun, x_data, y_data, p0, bounds=bounds)
    p_std = np.sqrt(np.diag(cov))
    return p, p_std


class Hist2d:

    def __init__(self, data, verbose=False):
        """
        2D histogram class

        Attributes
        ----------
        data : dict
               dictionary of 2D raw data ordered by contour, image, hex and scatter keys
        """
        self.contour = {}
        self.image = {}
        self.hex = {}
        self.scatter = {}
        self.XY1D = {}
        self.XY1DfitSum = {}
        self.XY1DfitComp = {}
        self.limits = {}
        self.XY1DfitParam_PAMtools = {}
        self.XY1DfitParamStd_PAMtools = {}
        self.XY1DfitSum_PAMtools = {}
        self.XY1DfitComp_PAMtools = {}
        self.parameters = {}
        additional_parameters = ['photons_per_window', 'crosstalk', 'direct_excitation', 'gamma_factor', 'donor_lifetime']
        for p in additional_parameters:
            try:
                self.parameters[p] = data[p]
            except KeyError:
                if verbose:
                    print('\"{}\" is not specified.'.format(p))
        for i in ['x', 'y', 'z']:
            self.contour[i] = np.array(data['contour'][i]) if data['contour'][i] is not None else None
            self.image[i] = np.array(data['image'][i]) if data['image'][i] is not None else None
            self.hex[i] = np.array(data['hex'][i]) if data['hex'][i] is not None else None
            if i != 'z':
                self.scatter[i] = data['scatter'][i] if data['scatter'][i] is not None else None
                self.XY1D[('X', i)] = np.array(data['X1D'][i])
                self.XY1D[('Y', i)] = np.array(data['Y1D'][i])
                self.XY1DfitSum[('X', i)] = np.array(data['X1DfitSum'][i])
                self.XY1DfitSum[('Y', i)] = np.array(data['Y1DfitSum'][i])
                self.XY1DfitComp[('X', i)] = np.array(data['X1DfitComp'][i])
                self.XY1DfitComp[('Y', i)] = np.array(data['Y1DfitComp'][i])
                self.limits[i] = data['limits'][i]
        self.contour['levels'] = np.array(data['contour']['levels']) if data['contour']['levels'] is not None else None
        self.cmapLimits = data['cmapLimits']

        if self.hex['x'] is not None:
            self.hex['patches'] = self.makeHex()

        self.sigma_shotnoise()


    def makeHex(self):
        """
        Make hexagons patches

        Returns
        -------
        patches : list of matplotlib.patches.Polygon
        """
        patches = []
        for i in range(self.hex['x'].shape[1]):
            patches.append(Polygon(np.vstack((self.hex['x'][:, i], self.hex['y'][:, i])).T))
        return patches

    def plot2Dhist(self, style='contour', axis_labels=('FRET', 'Stoichiometry'), cmap='RdGy_r', imgOffset=0, PAM_fit=True, PAMtools_fit=False, label=None, hist_color=[0.7, 0.7, 0.7], show_components=True, edgecolor=None, linewidth=0):
        """
        Display a 2D contour / image / hex or scatter plot

        Parameters
        ----------
        style : str (optional)
                style of the plot ('contour', 'image', 'hex', 'scatter'). Default='contour'
        axis_labels : tuple of str
                      labels of the x/y-axis
        cmap : colormap or str
               colormap object or name of a registered colormap
        imgOffset : int (optional)
                    Lower clipping offset (in %) for image plots
        PAM_fit : bool (default=True)
        PAMtools_fit : bool (default=False)
        label : str
        hist_color : rgb array
        """
        imgOff_counts = imgOffset * self.image['z'].max() / 100

        with sns.axes_style('ticks'):
            set_ticksStyle()
            f, ax = plt.subplots(nrows=2, ncols=2, figsize=(3, 3), sharex=False, sharey=False, squeeze=False, gridspec_kw={'width_ratios': [1, 0.3], 'height_ratios': [0.3, 1]})
            plt.subplots_adjust(wspace=0.05)
            plt.subplots_adjust(hspace=0.05)
            if style == 'image':
                img = np.ma.masked_where(self.image['z'] <= imgOff_counts, self.image['z'])  # zeros are always masked
                ax[1, 0].imshow(img, cmap=cmap, origin='lower', extent=[min(self.image['x']), max(self.image['x']), min(self.image['y']), max(self.image['y'])])
            elif style == 'hex':
                p = PatchCollection(self.hex['patches'], cmap=cmap, edgecolor=None)
                p.set_array(self.hex['z'])
                p.set_edgecolor('None')
                ax[1, 0].add_collection(p)
            elif style == 'scatter':
                ax[1, 0].plot(self.scatter['x'], self.scatter['y'], '.', markersize=1, color='black')
            else:
                ax[1, 0].contourf(self.contour['x'], self.contour['y'], self.contour['z'], levels=self.contour['levels'], cmap=cmap)
            ax[0, 0].bar(self.XY1D[('X', 'x')], self.XY1D[('X', 'y')], width=self.XY1D[('X', 'x')][1] - self.XY1D[('X', 'x')][0], color=hist_color, edgecolor=edgecolor, linewidth=linewidth)
            ax[1, 1].barh(self.XY1D[('Y', 'x')], self.XY1D[('Y', 'y')], height=self.XY1D[('Y', 'x')][1] - self.XY1D[('Y', 'x')][0], color=hist_color, edgecolor=edgecolor, linewidth=linewidth)

            if PAM_fit:
                for i in range(self.XY1DfitComp[('X', 'x')].shape[0]):
                    ax[0, 0].plot(self.XY1DfitComp[('X', 'x')][i, :], self.XY1DfitComp[('X', 'y')][i, :], color=colors[i], linestyle='--')
                    ax[1, 1].plot(self.XY1DfitComp[('Y', 'y')][i, :], self.XY1DfitComp[('Y', 'x')][i, :], color=colors[i], linestyle='--')
                ax[0, 0].plot(self.XY1DfitSum[('X', 'x')], self.XY1DfitSum[('X', 'y')], color='black')
                ax[1, 1].plot(self.XY1DfitSum[('Y', 'y')], self.XY1DfitSum[('Y', 'x')], color='black')

            if label:
                ax[1, 0].text(self.limits['x'][0] + np.diff(self.limits['x']) / 10, self.limits['y'][1] - np.diff(self.limits['y']) / 7, label, horizontalalignment='left')

            if PAMtools_fit and 'X' in self.XY1DfitParam_PAMtools:
                n = len(self.XY1DfitParam_PAMtools['X']['mu'])
                if n > 1 and show_components:
                    for i in range(n):
                        ax[0, 0].plot(self.XY1DfitSum_PAMtools[('X', 'x')], self.XY1DfitComp_PAMtools[('X', 'y', i)], ':', color=colors[i], linewidth=1.5)
                ax[0, 0].plot(self.XY1DfitSum_PAMtools[('X', 'x')], self.XY1DfitSum_PAMtools[('X', 'y')], color='black')

            if PAMtools_fit and 'Y' in self.XY1DfitParam_PAMtools:
                m = len(self.XY1DfitParam_PAMtools['Y']['mu'])
                if m > 1 and show_components:
                    for i in range(m):
                        ax[1, 1].plot(self.XY1DfitComp_PAMtools[('Y', 'y', i)], self.XY1DfitSum_PAMtools[('Y', 'x')], ':', color=colors[i], linewidth=1.5)
                ax[1, 1].plot(self.XY1DfitSum_PAMtools[('Y', 'y')], self.XY1DfitSum_PAMtools[('Y', 'x')], color='black')

            ax[0, 1].set_axis_off()
            ax[0, 0].set_axis_off()
            ax[1, 1].set_axis_off()
            ax[1, 0].set_xlim(self.limits['x'])
            ax[1, 0].set_ylim(self.limits['y'])
            ax[0, 0].set_xlim(self.limits['x'])
            ax[1, 1].set_ylim(self.limits['y'])
            ax[1, 0].set_xlabel(axis_labels[0])
            ax[1, 0].set_ylabel(axis_labels[1])
        self.ax = ax

    def plot1Dhist(self, axis_labels=('FRET', 'occupancy'), PAM_fit=True, PAMtools_fit=False, label=None, hist_color=[0.7, 0.7, 0.7], show_components=True):
        """
        Display a FRET histogram
        """
        with sns.axes_style('ticks'):
            set_ticksStyle()
            f, ax = plt.subplots(nrows=1, ncols=1, figsize=(2.5,2), sharex=False, sharey=False, squeeze=False)
            ax[0, 0].bar(self.XY1D[('X', 'x')], self.XY1D[('X', 'y')], width=(self.XY1D[('X', 'x')][1] - self.XY1D[('X', 'x')][0]), color=hist_color, linewidth=0)

            if PAM_fit:
                for i in range(self.XY1DfitComp[('X', 'x')].shape[0]):
                    ax[0, 0].plot(self.XY1DfitComp[('X', 'x')][i, :], self.XY1DfitComp[('X', 'y')][i, :], color=colors[i], linestyle='--')
                ax[0, 0].plot(self.XY1DfitSum[('X', 'x')], self.XY1DfitSum[('X', 'y')], color='black')

            if label:
                ax[1, 0].text(self.limits['x'][0] + np.diff(self.limits['x']) / 10, self.limits['y'][1] - np.diff(self.limits['y']) / 7, label, horizontalalignment='left')

            if PAMtools_fit and 'X' in self.XY1DfitParam_PAMtools:
                n = len(self.XY1DfitParam_PAMtools['X']['mu'])
                if n > 1 and show_components:
                    for i in range(n):
                        ax[0, 0].plot(self.XY1DfitSum_PAMtools[('X', 'x')], self.XY1DfitComp_PAMtools[('X', 'y', i)], ':', color=colors[i], linewidth=1.5)
                ax[0, 0].plot(self.XY1DfitSum_PAMtools[('X', 'x')], self.XY1DfitSum_PAMtools[('X', 'y')], color='black')

            ax[0, 0].set_xlim(self.limits['x'])
            ax[0, 0].set_xlabel(axis_labels[0])
            ax[0, 0].set_ylabel(axis_labels[1])
        self.ax = ax

    def plot_BVA(self, x_axis, x_axis_label=None, style='contour', cmap='RdGy_r', imgOffset=0, PAM_fit=False, PAMtools_fit=False, label=None, hist_color=[0.7, 0.7, 0.7], line_color='black'):
        """
        Display a 2D contour / image / hex or scatter plot of the Burst Variance Analysis

        x_axis : str
                 display FRET or proximity ratio on the x-axis('FRET', 'PR')
        x_axis_label : str
        style : str (optional)
                style of the plot ('contour', 'image', 'hex', 'scatter'). Default='contour'
        cmap : colormap or str
               colormap object or name of a registered colormap
        imgOffset : int (optional)
                    Lower clipping offset (in %) for image plots
        PAM_fit : bool (default=False)
        PAMtools_fit : bool (default=False)
        label : str
        hist_color : rgb array or str
        line_color : rgb array or str
        """
        if x_axis_label is None:
            x_axis_label = x_axis
        self.plot2Dhist(style=style, axis_labels=(x_axis_label, 'BVA st.dev.'), cmap=cmap, imgOffset=imgOffset, PAM_fit=PAM_fit, PAMtools_fit=PAMtools_fit, label=label, hist_color=hist_color)
        try: 
            self.ax[1, 0].plot(self.BVA[x_axis], self.BVA['sigma'], color=line_color)
        except ValueError:
            print('Parameters are missing. Line from BVA can not be drawn.')


    def fit_histogram(self, axis, mu0=[0.5], sigma0=[0.1], mu_bounds=(0, np.inf), sigma_bounds=(0, np.inf), verbose=True, fit_function='gaussian'):
        """
        Wrapper for the curve_fit function of the scipy.optimize module
        The curve_fit optimizes the Gaussian parameters (mu, sigma)
        while the nnls optimizes the relative weights of the Gaussians.

        Parameters
        ----------
        fun : callable
              The model function f(x,...) taking x values as a first argument followed by the function parameters
        axis : int
               Axis to fit (0: x-axis, 1: y-axis)
        mu0 : array_like
             start values for the Gaussian center
        sigma0 : array_like
                 start values for the standard deviation of the Gaussians
        mu_bounds : 2-tuple of float or 2-tuple of array_like, optional
                 lower and upper bounds for each parameter in mu0. Can be either a tuple of two scalars
                 (same bound for all parameters) or a tuple of array_like with the same length as mu0.
                 To deactivate parameter bounds set: `bounds=(-np.inf, np.inf)`
        sigma_bounds : 2-tuple of float or 2-tuple of array_like, optional
                 lower and upper bounds for each parameter in sigma0. See also mu_bounds
        """
        if axis == 0:
            a = 'X'
        else:
            a = 'Y'
        p0 = []
        gauss_bounds = ([], [])
        for i, (m, s) in enumerate(zip(mu0, sigma0)):
            p0.append(m)
            p0.append(s)
            if np.isscalar(mu_bounds[0]):
                mbl = mu_bounds[0]
            else:
                mbl = mu_bounds[0][i]

            if np.isscalar(mu_bounds[1]):
                mbu = mu_bounds[1]
            else:
                mbu = mu_bounds[1][i]

            if np.isscalar(sigma_bounds[0]):
                sbl = sigma_bounds[0]
            else:
                sbl = sigma_bounds[0][i]

            if np.isscalar(sigma_bounds[1]):
                sbu = sigma_bounds[1]
            else:
                sbu = sigma_bounds[1][i]

            gauss_bounds[0].append(mbl)
            gauss_bounds[0].append(sbl)
            gauss_bounds[1].append(mbu)
            gauss_bounds[1].append(sbu)


        if fit_function == 'beta':
            y_data = self.XY1D[a, 'y'][(self.XY1D[a, 'x']>0) & (self.XY1D[a, 'x']<1)]
            x_data = self.XY1D[a, 'x'][(self.XY1D[a, 'x']>0) & (self.XY1D[a, 'x']<1)]
        else:
            y_data = self.XY1D[a, 'y']
            x_data = self.XY1D[a, 'x']

        self._fitfunction = fit_function
        self._y_data = y_data
        p, p_std = fit(self._model_func, x_data, y_data, p0, bounds=gauss_bounds)
        A, x, y = self.nnls_convol_irfexp(x_data, p)
        self.XY1DfitParam_PAMtools[a] = {'ampl': x / sum(x), 'mu': p[0::2], 'sigma': p[1::2]}
        self.XY1DfitParamStd_PAMtools[a] = {'mu': p_std[0::2], 'sigma': p_std[1::2]}
        halfbinwidth = (x_data[1]-x_data[0])/2
        self.XY1DfitSum_PAMtools[a, 'x'] = np.linspace(x_data[0]-halfbinwidth, x_data[-1]+halfbinwidth, 200)
        A_fit = self.prepare_distributions(self.XY1DfitSum_PAMtools[a, 'x'], p)
        for i in range(len(x)):
            self.XY1DfitComp_PAMtools[a, 'y', i] = np.dot(A_fit[:, i], x[i])
        self.XY1DfitSum_PAMtools[a, 'y'] = np.dot(A_fit, np.array(x))

    @staticmethod
    def gauss(x_data, mu, sigma):
        """
        Calculate a Gaussian PDF

        Parameters
        ----------
        x_data : array_like
                 array of the independent variable
        mu : float
             mean of the Gaussian distribution
        sigma : float, optional
                standard deviation of the Gaussian distribution

        Returns
        -------
        gaussian : ndarray
        """
        gaussian = 1/(sigma*np.sqrt(2*np.pi))*np.exp(-(x_data - mu)**2 / (2 * sigma**2)).T
        return gaussian

    @staticmethod
    def beta(x_data, mu, sigma):
        """
        Calculate a beta PDF

        Parameters
        ----------
        x_data : array_like
                 array of the independent variable
        mu : float
             mean of the beta distribution
        sigma : float, optional
                standard deviation of the beta distribution

        Returns
        -------
        beta : ndarray

        References
        ----------
        Gopich and Szabo, Theory of Single-Molecule FRET Efficiency Histograms, Wiley (2011)

        """
        A = mu**2*(1-mu)/sigma**2-mu
        D = mu*(1-mu)**2/sigma**2-1+mu
        beta = scipy.special.gamma(A+D)/(scipy.special.gamma(A)*scipy.special.gamma(D))*x_data**(A-1)*(1-x_data)**(D-1)
        return beta

    def prepare_distributions(self, x_data, p0):
        distributions = []
        for k in range(0, len(p0), 2):
            if self._fitfunction == 'beta':
                distributions.append(self.beta(x_data, *p0[k:k + 2]))
            else:
                distributions.append(self.gauss(x_data, *p0[k:k + 2]))
            A = np.array(distributions).T
        return A

    def nnls_convol_irfexp(self, x_data, p0):
        """
        Solve non-negative least squares for series of IRF-convolved single-exponential decays.
        First, the IRF is shifted, then convolved with each exponential decay individually (decays 1,...,n),
        merged into an m x n array (=A) and finally plugged into scipy.optimize.nnls(A, experimental y-data) to which
        compute `argmin_x || Ax - y ||_2`.

        Parameters
        ----------
        x_data : array_like
                 array of the independent variable
        p0 : array_like
             start values for the fit model

        Returns
        -------
        A : ndarray
            matrix containing irf-convoluted single-exponential decays in the first n columns
            and ones in the last column (background counts)
        x : ndarray
            vector that minimizes `|| Ax - y ||_2`
        y : ndarray
            fit vector computed as `y = Ax`

        """
        A = self.prepare_distributions(x_data, p0)
        x, rnorm = scipy.optimize.nnls(A, self._y_data)
        y = np.dot(A, np.array(x))
        return A, x, y

    def _model_func(self, x_data, *p0):
        """
        Wrapper function for nnls_irfshift_convol

        Parameters
        ----------
        x_data : array_like
                 array of the independent variable
        p0 : array_like
             start values for the fit model

        See also
        --------
        nnls_convol_irfexp : Calculate non-linear least squares of IRF-convolved single-exponential decays

        Returns
        -------
        y : ndarray
            fit vector computed as `y = Ax`
        """
        A, x, y = self.nnls_convol_irfexp(x_data, p0)
        return y

    def PR2FRET(self, proximity_ratio, verbose=False):
        """
        Convert proximity ratio into absolute FRET efficiency

        Parameters
        ----------
        proximity_ratio : array_like
        alpha : float
                donor-acceptor spectral crosstalk (percentage of donor emission into acceptor detection channel)
        gamma : float
                difference in quantum yield and detection efficieny of donor and acceptor
        delta : float
                percentage of direct acceptor excitation by the green laser

        Returns
        -------
        FRET : array_like
        
        Reference
        ---------
        """
        try:
            FRET = (1-(1+self.parameters['crosstalk']+self.parameters['direct_excitation'])*(1-proximity_ratio))/(1-(1+self.parameters['crosstalk']-self.parameters['gamma_factor'])*(1-proximity_ratio))
        except KeyError:
            FRET = None
            if verbose:
                print('Parameter \"crosstalk\", \"direct_excitation\" or \"gamma_factor\" is missing.')
        return FRET


    def sigma_shotnoise(self, verbose=False):
        """
        Calculate shot noise limit from Burst Variance Analysis (BVA)
        """
        self.BVA = {}
        self.BVA['PR'] = np.linspace(0,1,1000)
        self.BVA['FRET']  = self.PR2FRET(self.BVA['PR'])
        try:
            self.BVA['sigma'] = np.sqrt(self.BVA['PR']*(1-self.BVA['PR'])/self.parameters['photons_per_window'])
        except KeyError:
            self.BVA['sigma'] = None
            if verbose:
                print('Parameter \"photons_per_window\" is missing. Expected proximity ratio variance can not be calculated.')

    def read_fit(self, filename):
        """
        Read PAM fit results from file

        Parameters
        ----------
        filename : str
        """
        self.XY1DfitParam = pd.read_csv(filename, header=2, sep='\t')

    @classmethod
    def fromfile(cls, filename, verbose=False):
        """
        Alternative constructor for the pamplot.Hist2d class

        Parameters
        ---------
        filename : str
        """
        with open(filename) as f:
            data = json.load(f)
        return cls(data, verbose)

class FCS:

    def __init__(self, data, parameters, verbose=False):
        """
        FCS class

        Attributes
        ----------
        data : pandas.DataFrame
               Dataframe with columns ['time', 'data', 'error', 'fit', 'res']
        """
        self.data = data
        self.parameters = parameters
        try:
            self.calcBrightness()
        except KeyError:
            if verbose:
                print('\"average counts\" (kHz) is not present in parameter file.')
        self.calcConcentration()
        self.parameters['Rh(298K)_nm'] = self.hydrodynamic_radius(parameters['D'])

    @classmethod
    def fromfile(cls, filename, verbose=False):
        """
        Alternative constructor for the pamplot.FCS class

        Parameters
        ---------
        filename : str
        verbose : boolean
        """
        data = pd.read_csv(filename, sep='\t')
        parameters = pd.read_csv('{}_param.txt'.format(filename[:-4]), sep='\t', header=2).iloc[0]
        if np.isnan(parameters['y0']):
            parameters['y0'] = 0
        return cls(data, parameters, verbose)

    @staticmethod
    def hydrodynamic_radius(D, T=298.15, eta=8.9*10**-4):
        """
        Calculate the hydrodnamic radius in nanometer

        Parameters
        ----------
        T : float
            temperature in Kelvin   
        eta : float
              viscosity of the solvent in Pa*s
        """
        return 1.38*10**-23*T/(6*np.pi*eta*D*10**-12)*10**9

    def calcConcentration(self):
        """
        Calculate the concentration in nanomolar
        """
        V_eff = np.pi**(3/2)*self.parameters['w_xy']**2*self.parameters['w_z']*10**-15
        N_A = 6.022*10**23
        self.parameters['c_nM'] = self.parameters['N']/(V_eff*N_A*10**-9)
        return '{:0.2f} nM'.format(self.parameters['c_nM'])

    def calcBrightness(self):
        """
        Calculate the molecular brightness in kHz per molecule
        """
        self.parameters['B'] = self.parameters['average_counts']/self.parameters['N']



    def plotFCS(self, color=ncp.get_colors(bwo, 1), time_limits=(10**-7, 0.1), legend_strings=None, normalization='None', normalization_time=10**-6, linecolor='k', figsize=(2.25, 2)):
        """
        Plot FCS curve

        Parameters
        ----------
        color : array_like 
                colors of FCS curve
        time_limits : tuple
        legend_strings : list of str
        normalization : str
                        normalization method {'None', 'G0', 'time'} (default='None')
                        G0 normalizes to fit
                        time normalzes to time point closest to reference time point given in normalization_time
        normalization_time : float
                             reference time point for normalization method 'time'
        """
        ax = plotFCS([self], color_list=color, time_limits=time_limits, legend_strings=legend_strings, normalization=normalization, linecolor=linecolor, figsize=figsize)
        return ax



def plotFCS(fcs_list, color_list=None, time_limits=(10**-7, 0.1), legend_strings=None, normalization='None', normalization_time=10**-6, linecolor='k', figsize=(2.25, 2)):
    """
    Plot FCS curve

    Parameters
    ----------
    fcs_list : list of pp.pamplot.FCS class
    color_list : array_like 
                 colors of FCS curves
    time_limits : tuple
                  time limits in seconds
    legend_strings : list of str
    normalization : str
                    normalization method {'None', 'G0', 'time'} (default='None')
                    G0 normalizes to fit
                    time normalzes to time point closest to reference time point given in normalization_time
    normalization_time : float
                         reference time point for normalization method 'time'
    """

    if color_list is None:
        color_list = ncp.get_colors(bwo, len(fcs_list))

    with sns.axes_style('ticks'):
        set_ticksStyle()
        f, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize, sharex=True, sharey=True, squeeze=False)
        for i,fcs in enumerate(fcs_list):
            time = fcs.data.loc[(fcs.data['time']>time_limits[0]) & (fcs.data['time']<time_limits[1]), 'time']
            data = fcs.data.loc[(fcs.data['time']>time_limits[0]) & (fcs.data['time']<time_limits[1]), 'data']
            fit = fcs.data.loc[(fcs.data['time']>time_limits[0]) & (fcs.data['time']<time_limits[1]), 'fit']
            if normalization == 'G0':
                first = fit.iloc[0]
                data = (data-fcs.parameters['y0'])/first
                fit = (fit-fcs.parameters['y0'])/first
                ax[0,0].set_ylabel('norm. G($\\tau$)')
            elif normalization == 'time':
                first = data.iloc[np.argmin(abs(time-normalization_time))]
                data = (data-fcs.parameters['y0'])/first
                fit = (fit-fcs.parameters['y0'])/first
                ax[0,0].set_ylabel('norm. G($\\tau$)')
            else:
                ax[0,0].set_ylabel('G($\\tau$)')
            ax[0,0].semilogx(time, data, '.', color=color_list[i])
            ax[0,0].semilogx(time, fit, '-', color=linecolor)
        
        ax[0,0].set_xlabel('lag time $\\tau$ (s)')
        ax[0,0].set_xlim(time_limits)
        if legend_strings is not None:
            ax[0,0].legend(legend_strings, frameon=False)
        return ax

def diffusion(tau, tau_D, N, w_xy, w_z):
    s = w_z/w_xy
    return 1/N*(1+tau/tau_D)**-1*(1+tau/(s**2*tau_D))**-(1/2)

def diffusion_triplet(tau, tau_D, N, w_xy, w_z):
    s = w_z/w_xy
    return 1/N*(1+T/(1-T)*np.exp(-tau/tau_T))*(1+tau/tau_D)**-1*(1+tau/(s**2*tau_D))**-(1/2)

class Timetrace:
    """
    FCS class

    Attributes
    ----------
    data : dict
           dictionary with keys ['time', 'counts']
    """
    def __init__(self, GG, GR, verbose=False):
        self.GG = pd.DataFrame(GG)
        self.GR =  pd.DataFrame(GR)

    @classmethod
    def fromfiles(cls, filename_GG, filename_GR, verbose=False):
        """
        Alternative constructor for the pamplot.Hist2d class

        Parameters
        ---------
        filename : list of str
        """
        with open(filename_GG) as f:
            GG = json.load(f)
        with open(filename_GR) as f:
            GR = json.load(f)
        return cls(GG, GR, verbose)


    def plotTrace(self, time_limits=[0,1], count_limits=[-50,50], figsize=(5, 2)):
        """
        Plot time trace

        Parameters
        ----------
        time_limits : tuple
                      time limits in seconds
        """
        with sns.axes_style('ticks'):
            set_ticksStyle()
            colors = [[0.21, 0.59, 0.32], [0.75, 0.25, 0.16]]
            f, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize, sharex=False, sharey=False, squeeze=False)
            ax[0,0].plot(self.GG.loc[(self.GG['time']>time_limits[0]) & (self.GG['time']<time_limits[1]), 'time'], self.GG.loc[(self.GG['time']>time_limits[0]) & (self.GG['time']<time_limits[1]), 'counts'], color=colors[0])
            ax[0,0].plot(self.GR.loc[(self.GR['time']>time_limits[0]) & (self.GR['time']<time_limits[1]), 'time'], -self.GR.loc[(self.GR['time']>time_limits[0]) & (self.GR['time']<time_limits[1]), 'counts'], color=colors[1])
            ax[0,0].set_ylabel('photon counts')
            ax[0,0].set_xlabel('time (s)')
            ax[0,0].set_xlim(time_limits)
            ax[0,0].set_ylim(count_limits)
        self.ax = ax
        