#!/usr/bin/env python3

import numpy as np
import json
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
import naturalcolors.colorpalette as ncp

# seaborn settings
sns.set_style('white')
sns.set_context("notebook")
sns.set(font='Arial')

cols = [[0.31, 0.45, 0.56], [0.6, 0.6, 0.6], [0.75, 0.51, 0.38]]
cmap_list, cmap_linseg = ncp.make_colormap(cols, 'ColorsOfNature')
colors = ncp.get_colors(cmap_linseg, 6, scramble=True)


def set_ticksStyle(x_size=4, y_size=4, x_dir='in', y_dir='in'):
    sns.set_style('ticks', {'xtick.major.size': x_size, 'ytick.major.size': y_size, 'xtick.direction': x_dir, 'ytick.direction': y_dir})


class Hist2d():

    def __init__(self, data):
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
        self.X1D = {}
        self.Y1D = {}
        self.X1DfitSum = {}
        self.Y1DfitSum = {}
        self.X1DfitComp = {}
        self.Y1DfitComp = {}
        self.limits = {}
        for i in ['x', 'y', 'z']:
            self.contour[i] = np.array(data['contour'][i]) if data['contour'][i] is not None else None
            self.image[i] = np.array(data['image'][i]) if data['image'][i] is not None else None
            self.hex[i] = np.array(data['hex'][i]) if data['hex'][i] is not None else None
            if i != 'z':
                self.scatter[i] = data['scatter'][i] if data['scatter'][i] is not None else None
                self.X1D[i] = np.array(data['X1D'][i])
                self.Y1D[i] = np.array(data['Y1D'][i])
                self.X1DfitSum[i] = np.array(data['X1DfitSum'][i])
                self.Y1DfitSum[i] = np.array(data['Y1DfitSum'][i])
                self.X1DfitComp[i] = np.array(data['X1DfitComp'][i])
                self.Y1DfitComp[i] = np.array(data['Y1DfitComp'][i])
                self.limits[i] = data['limits'][i]
        self.contour['levels'] = np.array(data['contour']['levels']) if data['contour']['levels'] is not None else None
        self.cmapLimits = data['cmapLimits']

        if self.hex['x'] is not None:
            self.hex['patches'] = self.makeHex()

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

    def plot2Dhist(self, style='contour', axis_labels=('FRET', 'Stoichiometry'), cmap='Blues_r', imgOffset=0):
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
            ax[0, 0].bar(self.X1D['x'], self.X1D['y'], width=self.X1D['x'][1] - self.X1D['x'][0], color=[0.7, 0.7, 0.7], linewidth=0)
            ax[1, 1].barh(self.Y1D['x'], self.Y1D['y'], height=self.Y1D['x'][1] - self.Y1D['x'][0], color=[0.7, 0.7, 0.7], linewidth=0)
            for i in range(self.X1DfitComp['x'].shape[0]):
                ax[0, 0].plot(self.X1DfitComp['x'][i, :], self.X1DfitComp['y'][i, :], color=colors[i], linestyle='--')
                ax[1, 1].plot(self.Y1DfitComp['y'][i, :], self.Y1DfitComp['x'][i, :], color=colors[i], linestyle='--')
            ax[0, 0].plot(self.X1DfitSum['x'], self.X1DfitSum['y'], color='black')
            ax[1, 1].plot(self.Y1DfitSum['y'], self.Y1DfitSum['x'], color='black')

            ax[0, 1].set_axis_off()
            ax[0, 0].set_axis_off()
            ax[1, 1].set_axis_off()
            ax[1, 0].set_xlim(self.limits['x'])
            ax[1, 0].set_ylim(self.limits['y'])
            ax[0, 0].set_xlim(self.limits['x'])
            ax[1, 1].set_ylim(self.limits['y'])
            ax[1, 0].set_xlabel(axis_labels[0])
            ax[1, 0].set_ylabel(axis_labels[1])

    @classmethod
    def fromfile(cls, filename):
        """
        Alternative constructor for the pamplot.Hist2d class

        Parameters
        ---------
        filename : str
        """
        with open(filename) as f:
            data = json.load(f)
        return cls(data)
