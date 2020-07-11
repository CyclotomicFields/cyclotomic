#!/usr/bin/env python3
import numpy as np
import matplotlib, glob, sys, itertools, math
import csv

# Required pre plt import configuration
matplotlib.rcParams['pgf.rcfonts'] = False
matplotlib.rcParams['pgf.texsystem'] = 'pdflatex'
matplotlib.rcParams['text.usetex'] = True
matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rcParams['figure.figsize'] = "4,4"

import matplotlib.pyplot as plt

from argparse import ArgumentParser
parser = ArgumentParser(description="Generate graphs from csv benchmark data")
parser.add_argument("-l", "--log", help="use a logarithmic scale for the y axis", action="store_true")
parser.add_argument("-p", "--prefix", help="prefix to use for the files generated", default="default_prefix")
parser.add_argument("-x", "--x", help="column to use as x axis and label", nargs=2, metavar=("column_index", "label"), required=True)
parser.add_argument("-y", "--y", help="column to use as y axis and label", nargs=2, metavar=("column_index", "label"), required=True)
parser.add_argument("-f", "--file", metavar=("label", "file.csv"), help="label and csv files to graph", nargs=2, action="append", required=True)
args = parser.parse_args()

NEXT_FIGURE = 1
log = args.log

def do_plot_multiple(results_mult, xindex, yindex, xlabel, ylabel, title, fname, labels):
    global NEXT_FIGURE, log
    marker_color = itertools.cycle(itertools.product(".x+", "rgbcmy"))
    plt.figure(NEXT_FIGURE)
    NEXT_FIGURE += 1
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    handles = []
    for i in range(len(results_mult)):
        to_plot = np.array(results_mult[i]).T
        if log:
            for j in range(len(to_plot[yindex])):
                to_plot[yindex][j] = math.log10(to_plot[yindex][j])
        marker, color = next(marker_color)
        handle, = plt.plot(to_plot[xindex], to_plot[yindex], label=labels[i], marker=marker, linestyle="None", color=color)
        handles.append(handle)
    lgd = plt.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    fig = plt.gcf()
    fig.savefig(fname+".pgf", bbox_extra_artists=[lgd], bbox_inches='tight')
    fig.savefig(fname+".png", bbox_extra_artists=[lgd], bbox_inches='tight')
    print(fname)

def read_pts(fname):
    data = []

    with open(fname, "r", newline="") as in_csv:
        reader = csv.reader(in_csv, delimiter=",")
        return np.array([[float(s) for s in row] for row in reader])

# plots multiple files on the same graph, for comparison
def plot_multiple(fnames, out_fname, header):
    y_label = header[1]
    if log:
        y_label = f"log({y_label})"
    ptss = [read_pts(fname[1]) for fname in fnames]
    do_plot_multiple(ptss, 0, 1, header[0], y_label, f'', f'{out_fname}', [fname[0] for fname in fnames])

def tolatex(data, header, output_file):
    output = ""
    output += '\\begin{tabular}{|'
    output += 'n{4}{2}|'*len(data[0])
    output += '}\n'
    output += "\\hline\n"
    output += " & ".join(['{{'+str(col)+'}}' for col in header]) + "\\\\\n"
    output += "\\hline\n"
    output += " \\\\\n".join([" & ".join([str(col) for col in row]) for row in data])
    output += " \\\\\n\\hline"
    output += '\n\\end{tabular}\n'
    open(output_file, "w").write(output)

header = [args.x[1], args.y[1]]

fnames = args.file
plot_multiple(fnames, args.prefix, header)

for (_, input_csv) in fnames:
    with open(input_csv, "r", newline="") as input_file:
        csv_reader = csv.reader(input_file, delimiter=",")
        results = [[int(x) for x in row] for row in csv_reader]
        tolatex(results, header, input_csv + ".tex")
