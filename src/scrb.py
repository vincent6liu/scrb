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
import csv


class SCRBGui(tk.Tk):
    def __init__(self, parent):
        tk.Tk.__init__(self, parent)
        self._parent = parent

        self.menubar = tk.Menu(self)
        self.fileMenu = tk.Menu(self.menubar, tearoff=0)

        self.currentPlot = None
        self.data = {}

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

        helv20 = font.Font(family='Helvetica', size=20)
        tk.Button(self, text="Load Files", font=helv20, command=self.load_files, height=3, width=10).grid(row=2)

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
        pass  # to be implemented

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
