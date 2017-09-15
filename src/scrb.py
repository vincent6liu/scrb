#!/usr/local/bin/python3

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from functools import reduce, partial
from mpl_toolkits.mplot3d import Axes3D  # necessary for 3D graph
import os
import platform
import pandas as pd
import tkinter as tk
import tkinter.font as font
import numpy as np
from tkinter import filedialog, ttk
from itertools import *
import csv
import matplotlib
import warnings
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import seaborn as sns
import re
from copy import deepcopy
from matplotlib.font_manager import FontProperties

# set plotting defaults
with warnings.catch_warnings():
    warnings.simplefilter('ignore')  # catch experimental ipython widget warning
    sns.set(context="paper", style='ticks', font_scale=1.5, font='Bitstream Vera Sans')

matplotlib.rcParams['image.cmap'] = 'viridis'
size = 12


def qualitative_colors(n):
    """ Generalte list of colors
    :param n: Number of colors
    """
    return sns.color_palette('Set1', n)


def get_fig(fig=None, ax=None, figsize=[6.5, 6.5]):
    """fills in any missing axis or figure with the currently active one
    :param ax: matplotlib Axis object
    :param fig: matplotlib Figure object
    """
    if not fig:
        fig = plt.figure(figsize=figsize)
    if not ax:
        ax = plt.gca()
    return fig, ax


def density_2d(x, y):
    """return x and y and their density z, sorted by their density (smallest to largest)

    :param x:
    :param y:
    :return:
    """
    xy = np.vstack([np.ravel(x), np.ravel(y)])
    z = gaussian_kde(xy)(xy)
    i = np.argsort(z)
    return np.ravel(x)[i], np.ravel(y)[i], np.arcsinh(z[i])



class SCRBGui(tk.Tk):
    def __init__(self, parent):
        tk.Tk.__init__(self, parent)
        self._parent = parent

        self.menubar = tk.Menu(self)
        self.fileMenu = tk.Menu(self.menubar, tearoff=0)

        self.currentPlot = None
        self.data = {'matrix': None, 'tsne': None, 'cluster': None, 'clusterlab':None, 'genelist': None}

        self.countmatrix_file = None
        self.cluster_file = None
        self.tsne_file = None
        self.genelist_file = None

        self.initialize()

    # updated
    def initialize(self):
        self.grid()

        # set up menu bar
        self.menubar.add_cascade(label="File", menu=self.fileMenu)
        self.fileMenu.add_command(label="Save Plot", state='disabled', command=self.save_plot)
        self.fileMenu.add_command(label="Exit", command=self.quit_scrb)

        self.config(menu=self.menubar)

        # intro screen
        tk.Label(self, text="SCRB", font=('Helvetica', 48), fg="black", bg="white", padx=100, pady=15).grid(row=0)
        tk.Label(self, text="Single Cell RNA Browser", font=('Helvetica', 25), fg="black",
                 bg="white", padx=100, pady=0).grid(row=1)

        self.helv20 = font.Font(family='Helvetica', size=20)
        tk.Button(self, text="Load Files", font=self.helv20, command=self.load_files, height=3, width=10).grid(row=2)

        # update
        self.protocol('WM_DELETE_WINDOW', self.quit_scrb)
        self.grid_columnconfigure(0, weight=1)
        self.resizable(False, False)
        self.update()
        self.geometry(self.geometry())
        self.focus_force()

    def load_files(self):
        self.import_files = tk.Toplevel()
        self.import_files.resizable(False, False)
        self.import_files.title('Import files')

        # get file path of count matrix
        matrixContainer = tk.Frame(self.import_files)
        matrixContainer.grid(column=0, row=0, sticky='w')
        tk.Label(matrixContainer, text="Count matrix file: ").grid(column=0, row=0, sticky='w')
        tk.Button(matrixContainer, text="Load", command=lambda: self.get_filename('count matrix')).grid(column=1,
                                                                                                        row=0,
                                                                                                        sticky='w')
        self.matrixfileVar = tk.StringVar()
        self.matrixfileVar.set("None selected")
        tk.Label(matrixContainer, textvariable=self.matrixfileVar).grid(column=2, row=0, sticky='e')

        # get file path of cluster
        clusterContainer = tk.Frame(self.import_files)
        clusterContainer.grid(column=0, row=1, sticky='w')
        tk.Label(clusterContainer, text="Cluster file: ").grid(column=0, row=0, sticky='w')
        tk.Button(clusterContainer, text="Load", command=lambda: self.get_filename('cluster')).grid(column=1, row=0,
                                                                                                    sticky='w')
        self.clusterfileVar = tk.StringVar()
        self.clusterfileVar.set("None selected")
        tk.Label(clusterContainer, textvariable=self.clusterfileVar).grid(column=2, row=0, sticky='e')

        # get file path of tsne
        tsneContainer = tk.Frame(self.import_files)
        tsneContainer.grid(column=0, row=2, sticky='w')
        tk.Label(tsneContainer, text="tSNE file: ").grid(column=0, row=0, sticky='w')
        tk.Button(tsneContainer, text="Load", command=lambda: self.get_filename('tsne')).grid(column=1, row=0,
                                                                                              sticky='w')

        self.tsnefileVar = tk.StringVar()
        self.tsnefileVar.set("None selected")
        tk.Label(tsneContainer, textvariable=self.tsnefileVar).grid(column=2, row=0, sticky='e')

        # get file path of gene list
        genelistContainer = tk.Frame(self.import_files)
        genelistContainer.grid(column=0, row=3, sticky='w')
        tk.Label(genelistContainer, text="Gene list file: ").grid(column=0, row=0, sticky='w')
        tk.Button(genelistContainer, text="Load", command=lambda: self.get_filename('gene list')).grid(column=1, row=0,
                                                                                                       sticky='w')
        self.genelistfileVar = tk.StringVar()
        self.genelistfileVar.set("None selected")
        tk.Label(genelistContainer, textvariable=self.genelistfileVar).grid(column=2, row=0, sticky='e')

        # final buttons
        finalButtonContainer = tk.Frame(self.import_files)
        finalButtonContainer.grid(column=0, row=4, sticky='w')
        tk.Button(finalButtonContainer, text="Cancel", command=self.import_files.destroy).grid(column=0, row=0, padx=35)
        tk.Button(finalButtonContainer, text="Process", command=self.process_data).grid(column=1, row=0, padx=35)

        self.wait_window(self.import_files)

    def get_filename(self, type: str):
        if type == "count matrix":
            self.countmatrix_file = filedialog.askopenfilename(title='Load count matrix file',
                                                               initialdir='~/.magic/data')
            self.matrixfileVar.set(self.countmatrix_file)
        elif type == "cluster":
            self.cluster_file = filedialog.askopenfilename(title='Load cluster file', initialdir='~/.magic/data')
            self.clusterfileVar.set(self.cluster_file)
        elif type == "tsne":
            self.tsne_file = filedialog.askopenfilename(title='Load tsne file', initialdir='~/.magic/data')
            self.tsnefileVar.set(self.tsne_file)
        elif type == "gene list":
            self.genelist_file = filedialog.askopenfilename(title='Load gene list file', initialdir='~/.magic/data')
            self.genelistfileVar.set(self.genelist_file)

    def process_data(self):

        # make sure all the data are loaded
        if None in {self.countmatrix_file, self.cluster_file, self.tsne_file, self.genelist_file}:
            warning_window = tk.Toplevel()
            tk.Label(warning_window, text="Warning: not all files are loaded!").grid(column=0, row=0)
            tk.Button(warning_window, text="Ok", command=warning_window.destroy).grid(column=0, row=1)
            return

        # read the data matrix and filter out genes with 0 reads
        matrix = pd.DataFrame.from_csv(self.countmatrix_file)
        sums = matrix.sum(axis=0)
        to_keep = np.where(sums > 0)[0]
        matrix = matrix.iloc[:, to_keep]
        self.data['matrix'] = matrix

        # read the cluster information
        with open(self.cluster_file) as cluster_file:
            clusters = cluster_file.readline()
            clusters = clusters.split(',')
            clusters[-1] = clusters[-1][0]
            clusters = pd.Series(clusters, index=matrix.index)
            self.data['cluster'] = clusters
            labels = []
            for line in cluster_file:
                newlab = re.split('[\t ,]', line)
                newlab[1] = newlab[1][:-1]
                newlab = tuple(newlab)
                labels.append(newlab)
            self.data['clusterlab'] = labels

        # read tsne information
        tsne = pd.DataFrame.from_csv(self.tsne_file)
        self.data['tsne'] = tsne

        # read gene list file
        gene_list = pd.DataFrame.from_csv(self.genelist_file)
        self.data['genelist'] = gene_list

        # close the file loading window
        self.import_files.destroy()

        # construct the main window
        for item in self.grid_slaves():
            item.grid_forget()

        # list of genes ranked by p-value
        self.genes_list = ttk.Treeview(height=34)
        self.genes_list.heading('#0', text='Genes')
        self.genes_list.grid(column=0, row=0, rowspan=6, sticky='NSEW')
        ysb = ttk.Scrollbar(orient=tk.VERTICAL, command=self.genes_list.yview)
        xsb = ttk.Scrollbar(orient=tk.HORIZONTAL, command=self.genes_list.xview)
        self.genes_list.configure(yscroll=ysb.set, xscroll=xsb.set)

        # option to visualize gene expression
        self.visual_button = tk.Button(text="Select gene(s)", command=self.exp_visual,
                                       font=self.helv20, height=5, width=30, wraplength=60)
        self.visual_button.grid(column=0, row=8, sticky='NSEW')

        self.notebook = ttk.Notebook(height=700, width=700)
        self.notebook.grid(column=1, row=0, rowspan=14, columnspan=4, sticky='NSEW')
        self.tabs = []

        # update
        self.geometry('1100x750')

        # visualize cluster on tSNE projection upon loading
        self._visualizeCluster()

        # add the gene list
        for index, row in gene_list.iterrows():
            entry = str(index).upper() + ' (p-value = ' + str(row['p-value']) + ')'
            self.genes_list.insert('', 'end', text=entry, open=True)

    def exp_visual(self):
        self.get_genes = tk.Toplevel()
        self.get_genes.resizable(False, False)
        self.get_genes.title("Select genes")

        entry = []
        for key in self.genes_list.selection():
            curgene = self.genes_list.item(key, 'text').split(' (')[0]
            entry.append(curgene)
        entry = ', '.join(entry)

        geneEntryContainer = tk.Frame(self.get_genes)
        geneEntryContainer.grid(column=0, row=0, sticky='w')
        self.geneVar = tk.StringVar()
        self.geneVar.set(entry)
        tk.Label(geneEntryContainer, text="Select one or more genes, separated by commas:").grid(column=0,
                                                                                                 row=0, sticky='w')
        tk.Entry(geneEntryContainer, textvariable=self.geneVar).grid(column=1, row=0, sticky='w')

        # final buttons
        finalButtonContainer = tk.Frame(self.get_genes)
        finalButtonContainer.grid(column=0, row=1, sticky='w')
        tk.Button(finalButtonContainer, text="Cancel", command=self.get_genes.destroy).grid(column=0, row=0, padx=35)
        tk.Button(finalButtonContainer, text="Process", command=self._exp_visual).grid(column=1, row=0, padx=35)

        self.wait_window(self.get_genes)

    def _exp_visual(self):
        self.get_genes.destroy()

        genes = re.split('[, ]+', self.geneVar.get().upper())
        num_genes, side = len(genes), 2
        if num_genes == 0:
            return

        fontP = FontProperties()
        fontP.set_size('xx-small')

        # get the correct side length
        while num_genes+1 > side**2:
            side += 1

        matrix = self.data['matrix']
        tsnedata = self.data['tsne']
        communities = self.data['cluster']

        fig, axarr = plt.subplots(side, side)
        fig.set_size_inches(7, 7)
        sc = axarr[0, 0].scatter(tsnedata['tSNE1'], tsnedata['tSNE2'], s=size,
                            c=communities.values, edgecolors='none', cmap='rainbow')
        axarr[0, 0].set(adjustable='box-forced')
        axarr[0, 0].set_title('Cluster', fontsize=12)

        """
        lp = lambda i: plt.plot([], color=sc.cmap(sc.norm(i)), ms=np.sqrt(size), mec="none",
                                label="Cluster {:g}".format(i), ls="", marker="o")[0]
        handles = [lp(int(i)) for i in np.unique(communities)]
        plt.legend(handles=handles, prop=fontP, bbox_to_anchor=(0, 1), loc='upper left', ncol=1).set_frame_on(True)
        """


        for i in range(side):
            for j in range(side):
                if i == 0 and j == 0:
                    for tick in axarr[i, j].xaxis.get_major_ticks():
                        tick.label.set_fontsize(10)
                    for tick in axarr[i, j].yaxis.get_major_ticks():
                        tick.label.set_fontsize(10)
                elif len(genes) != 0:
                    curgene = genes.pop(0)
                    expression = matrix[curgene]
                    expression = np.log10(expression + 0.1)
                    lgmatrix = np.log10(matrix.as_matrix().flatten() + 0.1)
                    sc = axarr[i, j].scatter(tsnedata['tSNE1'], tsnedata['tSNE2'], s=size, c=expression.values,
                                             edgecolors='none', cmap='coolwarm',
                                             vmin=min(lgmatrix), vmax=np.percentile(lgmatrix, 97))
                    axarr[i, j].set(adjustable='box-forced')
                    axarr[i, j].set_title(curgene, fontsize=12)

                    for tick in axarr[i, j].xaxis.get_major_ticks():
                        tick.label.set_fontsize(10)
                    for tick in axarr[i, j].yaxis.get_major_ticks():
                        tick.label.set_fontsize(10)

                    fig.colorbar(sc, ax=axarr[i, j], orientation='vertical')

                else:
                    axarr[i, j].set(adjustable='box-forced')
                    # turn off ones that do not display genes
                    axarr[i, j].axis('off')

        plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
        plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
        if side > 2:
            plt.setp([a.get_yticklabels() for a in axarr[:, 2]], visible=False)

        plt.subplots_adjust(wspace=None, hspace=0.3)

        self.tabs.append([tk.Frame(self.notebook), fig])
        self.notebook.add(self.tabs[len(self.tabs) - 1][0], text="Expression")

        self.canvas = FigureCanvasTkAgg(fig, self.tabs[len(self.tabs) - 1][0])
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=10, columnspan=4, sticky='NSEW')

        self.currentPlot = 'tsne'


    def _visualizeCluster(self):
        tsnedata = self.data['tsne']
        communities = self.data['cluster']

        if min(set(communities)) == 0:
            communities = [x+1 for x in communities]
        color = communities

        self.fig = plt.figure(figsize=[7, 7])
        gs = gridspec.GridSpec(1, 1)
        self.ax = self.fig.add_subplot(gs[0, 0])

        self.plot_tsne(tsnedata, self.fig, self.ax, color=color)
        self.ax.set_title('Cluster Visualization')
        self.ax.set_xlabel('tSNE1')
        self.ax.set_ylabel('tSNE2')

        gs.tight_layout(self.fig)

        self.tabs.append([tk.Frame(self.notebook), self.fig])
        self.notebook.add(self.tabs[len(self.tabs) - 1][0], text="Cluster")

        self.canvas = FigureCanvasTkAgg(self.fig, self.tabs[len(self.tabs) - 1][0])
        self.canvas.show()
        self.canvas.get_tk_widget().grid(column=1, row=1, rowspan=10, columnspan=4, sticky='NSEW')

        self.currentPlot = 'tsne'

    @staticmethod
    def plot_tsne(tsne, fig=None, ax=None, density=False, color=None, ge=False, title='tSNE projection'):
        """Plot tSNE projections of the data
        Must make sure the object being operated contains tSNE data
        :param tsne: pd.Dataframe that contains tsne data
        :param fig: matplotlib Figure object
        :param ax: matplotlib Axis object
        :param title: Title for the plot
        """
        fontP = FontProperties()
        fontP.set_size('xx-small')

        fig, ax = get_fig(fig=fig, ax=ax)
        if isinstance(color, pd.Series) and ge:
            sc = plt.scatter(tsne['tSNE1'], tsne['tSNE2'], s=size,
                             c=color.values, edgecolors='none', cmap='Oranges')
        elif isinstance(color, pd.Series) and not ge:  # cluster visualization
            sc = plt.scatter(tsne['tSNE1'], tsne['tSNE2'], s=size,
                             c=color.values, edgecolors='none', cmap='rainbow')
            lp = lambda i: plt.plot([], color=sc.cmap(sc.norm(i)), ms=np.sqrt(size), mec="none",
                                    label="Cluster {:g}".format(i), ls="", marker="o")[0]
            handles = [lp(int(i)) for i in np.unique(color)]
            plt.legend(handles=handles, prop=fontP, loc='upper right').set_frame_on(True)
        elif density:
            # Calculate the point density
            xy = np.vstack([tsne['tSNE1'], tsne['tSNE2']])
            z = gaussian_kde(xy)(xy)

            # Sort the points by density, so that the densest points are plotted last
            idx = z.argsort()
            x, y, z = tsne['tSNE1'][idx], tsne['tSNE2'][idx], z[idx]

            plt.scatter(x, y, s=size, c=z, edgecolors='none')
            plt.colorbar()
        else:
            plt.scatter(tsne['tSNE1'], tsne['tSNE2'], s=size, edgecolors='none',
                        color=qualitative_colors(2)[1] if color is None else color)

        ax.set_title(title)
        #plt.axis('tight')
        plt.tight_layout()
        return fig, ax

    def save_plot(self):
        tab = self.notebook.index(self.notebook.select())
        default_name = self.notebook.tab(self.notebook.select(), "text")

        plotFileName = filedialog.asksaveasfilename(title='Save Plot', defaultextension='.png',
                                                    initialfile=default_name)
        if plotFileName is not None:
            self.tabs[tab][1].savefig(plotFileName)

    def quit_scrb(self):
        self.quit()
        self.destroy()


def launch():
    app = SCRBGui(None)
    if platform.system() == 'Darwin':
        app.focus_force()
    elif platform.system() == 'Windows':
        app.lift()
        app.call('wm', 'attributes', '.', '-topmost', True)
        app.after_idle(app.call, 'wm', 'attributes', '.', '-topmost', False)
    elif platform.system() == 'Linux':
        app.focus_force()

    app.title('SCRB')

    while True:
        try:
            app.mainloop()
            break
        except UnicodeDecodeError:
            pass


if __name__ == "__main__":
    launch()
