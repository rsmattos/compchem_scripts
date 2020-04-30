#! /usr/bin/env python

import numpy as np
import math
from sys import argv
from os import system
import argparse

"""
Modification of the script g09_spectrum.py, that started this whole project, available at
https://github.com/mdommett/compchem-scripts/blob/master/g09_spectrum.py

This program plots a UV/Vis absorption spectrum from a Gaussian 09 output file or Orca. Guassian broadening is used
to generate peaks about the excitation energy. see http://gaussian.com/uvvisplot/ for implementation details.
Multiple output files can be plotted at the same time.
Usage:
spectrum.py <outputfile>
A number of options can be used to control the output:
-gnu : plot using gnuplot
-prog : specify the program that generated the output being read
-mpl: plot using matplotlib (recommended)
-sticks: plot excitation energies as a stick g09_spectrum
-sd : set the standard deviation (in eV). Default is 0.4
-save : save resultant specturm as pdf file
-dat : save the spectrum data as a text file
"""
parser = argparse.ArgumentParser()
parser.add_argument("input",help="Log file of Gaussian 09 TD job", type=str, nargs='*')
parser.add_argument("-prog",help="Specify the quantum chem program used",type=str)
parser.add_argument("-mpl",help="Plot a spectrum using matplotlib",action="store_true")
parser.add_argument("-sticks",help="Plot the stick spectrum",action="store_true")
parser.add_argument("-sd",help="Standard deviation (in eV)",default=0.4,type=float)
parser.add_argument("-rng",help="Min and max values for the spectrum (in nm)",nargs=2,type=int)
parser.add_argument("-save",help="Save spectrum with matplotlib",nargs='?',const="spectrum", type=str)
parser.add_argument("-dat",help="Save data as a text file",nargs='?',const="spectrum",type=str)
parser.add_argument("-csv",help="Save data as a csv file",nargs='?',const="spectrum",type=str)
parser.add_argument("-dist", help="Choose the type of function to calculate the spectra", type=str)
parser.add_argument("-shift", help="Adds red or blue shift to the spectra", default=0.0, type=float)
args=parser.parse_args()

######################### READING OUTPUT FILES ########################
def find_program(line):
    for i in range(5):
        if "Gaussian" in line[i]:
            return "gaussian"

        if "* O   R   C   A *" in line[i]:
            return "orca"

    print("Could not identify the program, try specifying with -prog.")
    quit()

def read_g09(line):
    energies=[]
    os_strengths=[]
    for i in range(len(line)):
        if " Excited State " in line[i]:
            energies.append(float(line[i].split()[6]))
            os_strengths.append(float(line[i].split()[8][2:]))
    return energies,os_strengths

def read_orca(line):
    energies=[]
    os_strengths=[]
    for i in range(len(line)):
        if " nroots " in line[i]:
            nroots=int(line[i].split()[3])
        if " ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" in line[i]:
            for j in range(i+5,i+5+nroots):
                if float(line[j].split()[3]) > 10E-4:
                    energies.append(float(line[j].split()[2]))
                    os_strengths.append(float(line[j].split()[3]))
            break
    return energies,os_strengths

def read_files(f):
    # read output files
    infile=open(f,"r")
    lines=infile.readlines()

    # find out the program that generated the output
    if not args.prog:
        args.prog = find_program(lines)

    # read the energies and oscilator strenghts from output file
    if args.prog=="orca":
        energies, os_strengths=read_orca(lines)
    elif args.prog=="gaussian":
        energies, os_strengths=read_g09(lines)
    else:
        print("Program not supported.")
        quit()

    infile.close()
    return energies, os_strengths
#######################################################################

######################### CALCULATE SPECTRA ###########################
def gaussian_dist(f,lamb_max,lamb):
    # 1240 passes ards.sd from nm-1 to eV
    return f*np.exp(-(( (1/lamb - 1/lamb_max)*1240 + args.shift)/args.sd)**2)

def calculate_spectra(energies, os_strengths):
    # set x axis
    if args.rng:
        x=np.linspace(max(args.rng),min(args.rng),1000)

    else:
        x=np.linspace(max(energies)+200,min(energies)-200,1000)

    # calculate y axis
    sum = np.zeros(len(x))
    for i in range(len(energies)):
        sum += gaussian_dist(os_strengths[i],energies[i], x)

    sum *= 0.619/args.sd

    return x, sum
#######################################################################

############################## OUTPUTTING #############################
def out_data(xaxis,yaxis):
    with open(args.dat+str(n)+".dat","w") as d:
        d.write("# index       Wavelenght  Oscilator_strenght\n")
        for i in range(len(xaxis)):
            d.write("{0:>7} {1:>16}   {2:>16}\n".format(i+1,xaxis[i],yaxis[i]))
        d.close()

def out_csv(xaxis,yaxis):
    with open(args.csv+str(n)+".csv","w") as d:
        d.write("index,Wavelenght,Oscilator_strenght\n")
        for i in range(len(xaxis)):
            d.write("{0},{1},{2}\n".format(i+1,xaxis[i],yaxis[i]))
        d.close()

def setup_plot():
    import matplotlib.pyplot as plt

    plt.style.use("seaborn-deep")

    fig , ax = plt.subplots()

    plt.rc('font', family='sans-serif')
    plt.tick_params(labelsize=10)

    # Axis labels
    ax.set_xlabel("Energy (nm)", fontsize=10)
    ax.set_ylabel("Oscilator strenght", fontsize=10)

    # Image size
    fig.set_size_inches(5.0, 5.0)

    return plt, fig, ax

def mpl_plot(fig, ax, x, sum):
    ax.plot(x,sum,label=f[:-4]+str(n))
    
    if args.sticks:
        for i in range(len(energies)):
            ax.vlines(energies[i],0,os_strengths[i])
    if args.rng:
        ax.xlim(min(args.rng),max(args.rng))
    ax.legend()

    return
#######################################################################

if __name__=='__main__':
    if args.mpl or args.save:
        plt, fig, ax = setup_plot()

    for n,f in enumerate(args.input):

        energies, os_strengths = read_files(f)

        x, sum = calculate_spectra(energies, os_strengths)

        # output formats
        if args.dat:
            out_data(x,sum)

        if args.csv:
            out_csv(x,sum)

        if args.mpl or args.save:
            mpl_plot(fig, ax, x, sum)
    
    if args.save:
        plt.savefig(args.save+".pdf")
    if args.mpl:
        plt.show()
