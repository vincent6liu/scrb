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
        self.fileMenu.add_command(label="Exit", command=self.quit_scrb())

        self.config(menu=self.menubar)

        # intro screen
        tk.Label(self, text="SCRB", font=('Helvetica', 48), fg="black", bg="white", padx=100, pady=15).grid(row=0)
        tk.Label(self, text="Single Cell RNA Browser", font=('Helvetica', 25), fg="black",
                 bg="white", padx=100, pady=0).grid(row=1)

        self.helv20 = font.Font(family='Helvetica', size=20)
        tk.Button(self, text="Load Files", font=self.helv20, command=self.load_files, height=3, width=10).grid(row=2)

        # update
        self.protocol('WM_DELETE_WINDOW', self.quit_scrb())
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
                newlab = line.split('\t')
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
        self.genes_list = ttk.Treeview(height=30)
        self.genes_list.heading('#0', text='Genes')
        self.genes_list.grid(column=0, row=0, rowspan=6, sticky='NSEW')
        ysb = ttk.Scrollbar(orient=tk.VERTICAL, command=self.genes_list.yview)
        xsb = ttk.Scrollbar(orient=tk.HORIZONTAL, command=self.genes_list.xview)
        self.genes_list.configure(yscroll=ysb.set, xscroll=xsb.set)

        # option to visualize gene expression
        self.visual_button = tk.Button(text="Select gene(s)", command=self.exp_visual,
                                       font=self.helv20, height=5, width=30)
        self.visual_button.grid(column=0, row=8, sticky='NSEW')

        self.notebook = ttk.Notebook(height=600, width=600)
        self.notebook.grid(column=1, row=0, rowspan=14, columnspan=4, sticky='NSEW')
        self.tabs = []

        # update
        self.geometry('1000x650')

    def exp_visual(self):
        pass

    def save_plot(self):
        pass  # to be implemented

    def quit_scrb(self):
        pass  # to be implemented


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
