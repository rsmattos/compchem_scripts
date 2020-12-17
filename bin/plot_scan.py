#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

def plotting(args, state, Plotting_var):

    fig , ax = plt.subplots()

    # Axis ticks
    ax.xaxis.set_tick_params(top=False, direction='out', width=1)
    ax.xaxis.set_tick_params(bottom=True, direction='in', width=1)
    ax.yaxis.set_tick_params(right=False, direction='in', width=1)
    ax.yaxis.set_tick_params(bottom=True, direction='in', width=1)

    plt.rc('font', family='sans-serif')
    plt.tick_params(labelsize=Plotting_var.font_size*0.9)

    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))

    # Axis labels
    ax.set_ylabel(Plotting_var.yaxis_label, fontsize=Plotting_var.font_size)

    if args.bond:
        ax.set_xlabel("Distance (Angstron)", fontsize=Plotting_var.font_size)
    elif args.angle:
        ax.set_xlabel("Angle (Degree)", fontsize=Plotting_var.font_size)
    elif args.dihedral:
        ax.set_xlabel("Dihedral angle (Degree)", fontsize=Plotting_var.font_size)
    elif args.general:
        ax.set_xlabel("Steps", fontsize=Plotting_var.font_size)

    # Image size
    fig.set_size_inches(Plotting_var.size_x, Plotting_var.size_y)

    # Coloring
    num_cols = len(state.columns)
    if num_cols < 9:
        colors = plt.cm.Dark2.colors
    elif num_cols < 10:
        colors = plt.cm.Set1.colors
    elif num_cols < 11:
        colors = plt.cm.tab10.colors
    else:
        colors = plt.cm.tab20.colors
    
    # Plotting
    if (args.plot_linestyle is not None) and (len(args.plot_linestyle) % 3) != 0:
        print("The parameter --plot_linestyle has to be given in triples, the first element is the irrep (a, b) or multiplicity ({1}, {2}), the second is the marker style, the third is the line style.")
        exit()

    i = 0

    for S in state.columns:
        marker = 'o'
        style = '-'
        if args.plot_linestyle is not None:
            for key in range(int(len(args.plot_linestyle)/3)):
                if args.plot_linestyle[3*key] in S:
                    marker = args.plot_linestyle[3*key + 1] 
                    style = args.plot_linestyle[3*key + 2]

        ax.plot(state[S], \
                label=str(S), \
                marker=marker, \
                markersize=3, \
                color=colors[i], \
                linewidth=2, \
                linestyle = style)
        i += 1
        if i > 19:
            break

    # Legend
    if num_cols < 7:
        ncols = num_cols
    elif num_cols < 13:
        ncols = int(num_cols/2) + (num_cols % 2 > 0)
    elif num_cols < 19:
        ncols = int(num_cols/3) + (num_cols % 3 > 0)
    else:
        ncols = 6

    box = ax.get_position()
    ax.set_position([box.x0 + box.width*0.05, \
                     box.y0 + box.height*(0.05+Plotting_var.legend_height+int(num_cols/ncols)/30), \
                     box.width, \
                     box.height*(1-Plotting_var.legend_height)])
   
    ax.legend(loc='upper center', \
              bbox_to_anchor=(0.45, -Plotting_var.legend_height -int(num_cols/ncols)/30), \
              fancybox=True, \
              shadow=True, \
              ncol=ncols, \
              fontsize=Plotting_var.font_size*0.9)

    # Saving
    if args.save:
        diction = plt.gcf().canvas.get_supported_filetypes()

        if(args.save in diction):
            fig.savefig(args.output+Plotting_var.saving_name+args.save)
        else:
            print("Graphic plot type not supported, the available formats are:")
            for key in diction:
                print (key, " => ", diction[key])

    if args.noshow:
        return

    plt.show()

    return

def energy(args, energy):
    class Plotting_var:
        yaxis_label = r"$\Delta$E (eV)"
        size_x = 6.0
        size_y = 6.0
        font_size = 11.0
        legend_height = 0.1
        saving_name = '.energy.'

    plotting(args, energy, Plotting_var)

def oscillator_str(args, oscillator_strenght):
    class Plotting_var:
        yaxis_label = "Oscillator Strenght"
        size_x = 6.0
        size_y = 4.0
        font_size = 11.0
        legend_height = 0.2
        saving_name = '.oscillator_str.'

    plotting(args, oscillator_strenght, Plotting_var)

def charge_transfer(args, charge_tranfer):
    class Plotting_var:
        yaxis_label = "Charge Transfer"
        size_x = 6.0
        size_y = 4.0
        font_size = 11.0
        legend_height = 0.2
        saving_name = '.charge_transfer.'

    plotting(args, charge_tranfer, Plotting_var)

###################      TESTING     ###################
if __name__=='__main__':
    import argparse
    import pandas as pd

    parser = argparse.ArgumentParser(description='''\
    Script to extrat the energies from a scan calculation. \n \
    When passing the atom index, the first atom given is connected to the second, which is connected to the third, and so on. \n \
    Starting to count with the first atom in the coordinates list being 1.''')
    parser.add_argument('-s','--save',help="Determine the output file format",type=str,nargs='?',const='pdf')
    parser.add_argument('-o','--output',help="Base name of the output file",type=str,nargs='?',default='scan')
    parser.add_argument('--noshow',help="Stop from plotting the graph in the screen",nargs='?',type=bool,const=True,default=False)
    parser.add_argument('-pls','--plot_linestyle', help='Change line style of the different symmetry or multiplicity curves', type=str, nargs='+')
    parser.add_argument('-e','--energy', help='Read and plot the enrgy file', type=str)
    parser.add_argument('-os','--oscillator', help='Read and plot the oscillator strength file', type=str)
    parser.add_argument('-ct','--charge_transfer', help='Read and plot the charge transfer file', type=str)
    parameter = parser.add_mutually_exclusive_group(required=True)
    parameter.add_argument('-b','--bond',help="Atom index to calculate the bond distanced.",type=bool,nargs='?',const=True,default=False)
    parameter.add_argument('-a','--angle',help="Atom index to calculate the angle.", type=bool,nargs='?',const=True,default=False)
    parameter.add_argument('-d','--dihedral',help="Atom index to calculate the dihedral angle.",type=bool,nargs='?',const=True,default=False)
    parameter.add_argument('-g','--general', help="When the modification isn't symple", action='store_true')
    args=parser.parse_args()

    if args.energy:
        energy_df = pd.read_csv(args.energy, index_col=0)
        print(energy_df)
        energy(args, energy_df)

    if args.oscillator:
        oscillator_df = pd.read_csv(args.oscillator, index_col=0)
        print(oscillator_df)
        oscillator_str(args, oscillator_df)

    if args.charge_transfer:
        charge_transfer_df = pd.read_csv(args.charge_transfer, index_col=0)
        print(charge_transfer_df)
        charge_transfer(args, charge_transfer_df)