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
