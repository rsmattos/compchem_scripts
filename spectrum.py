#! /usr/bin/env python

import numpy as np
import math
from sys import argv
from os import system
import argparse
#########
#########
"""
Modification of the script g09_spectrum.py, available at
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
-raw : save the spectrum data as a text file
"""
parser = argparse.ArgumentParser()
parser.add_argument("input",help="Log file of Gaussian 09 TD job", type=str, nargs='*')
parser.add_argument("-gnu",help="Plot a spectrum using gnuplot",action="store_true")
parser.add_argument("-prog",help="Specify the quantum chem program used",default="orca",type=str)
parser.add_argument("-mpl",help="Plot a spectrum using matplotlib",action="store_true")
parser.add_argument("-sticks",help="Plot the stick spectrum",action="store_true")
parser.add_argument("-sd",help="Standard deviation (in eV)",default=0.4,type=float)
parser.add_argument("-r",help="Min and max values for the spectrum (in nm)",nargs=2,type=int)
parser.add_argument("-save",help="Save spectrum", type=str)
parser.add_argument("-raw",help="Save raw data as text file",default="data",type=str)
args=parser.parse_args()

def read_g09(file):
    energies=[]
    os_strengths=[]
    line=file.readlines()
    for i in range(len(line)):
        if " Excited State " in line[i]:
            energies.append(float(line[i].split()[6]))
            os_strengths.append(float(line[i].split()[8][2:]))
    return energies,os_strengths

def read_orca(file):
    energies=[]
    os_strengths=[]
    line=file.readlines()
    for i in range(len(line)):
        if " nroots " in line[i]:
            nroots=int(line[i].split()[3])
        if " ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" in line[i]:
            for j in range(i+5,i+5+nroots):
                energies.append(float(line[j].split()[2]))
                os_strengths.append(float(line[j].split()[3]))
    return energies,os_strengths

def abs_max(f,lam,ref):
    a=1.3062974e8
    b=f/(1e7/3099.6)
    c=np.exp(-(((1/ref-1/lam)/(1/(1240/args.sd)))**2))
    return a*b*c

def raw_data(xaxis,yaxis):
    with open(args.raw+".dat","w") as d:
        for i in range(len(xaxis)):
            d.write("{0} {1} {2}\n".format(i+1,xaxis[i],yaxis[i]))
        d.close()

def gnu_plot(xaxis,yaxis):
    with open("data","w") as d:
        for i in range(len(xaxis)):
            d.write("{0} {1} {2}\n".format(i+1,xaxis[i],yaxis[i]))
        d.close()

    with open("plot","w") as p:
        p.write("set xlabel \"Energy (nm)\"\nset ylabel \"Absorption Coeff.(E)\"\n")
        p.write("plot 'data' using 2:3 title '' with lines lt 1 lw 2,\\\n")
        p.close()

    system("gnuplot -persist plot")
    system("rm -f plot data")

    return

def mpl_plot(xaxis,yaxis):
    plt.scatter(xaxis,yaxis,s=2,c="r")
    plt.plot(xaxis,yaxis,color="k")
    plt.xlabel("Energy (nm)")
    plt.ylabel("$\epsilon$ (L mol$^{-1}$ cm$^{-1}$)")
    return

if __name__=='__main__':
    for n,f in enumerate(args.input):
        infile=open(f,"r")
        if(args.prog=="orca"):
            energies,os_strengths=read_orca(infile)
        elif(args.prog=="gaussian"):
            energies,os_strengths=read_g09(infile)
        else:
        	print("Program not supported.")
        infile.close()

        if args.r:
            x=np.linspace(max(args.r),min(args.r),1000)

        else:
            x=np.linspace(max(energies)+200,min(energies)-200,1000)

        sum=[]
        for ref in x:
            tot=0
            for i in range(len(energies)):
                tot+=abs_max(os_strengths[i],energies[i],ref)
            sum.append(tot)

        stick_intensities=[abs_max(os_strengths[i],energies[i],energies[i]) for i in range(len(energies))]

        if args.raw:
            raw_data(x,sum)

        if args.gnu:
            gnu_plot(x,sum)

        else:
            import matplotlib.pyplot as plt
            colours=["red","blue","green","orange","black","cyan","magenta"]
            plt.scatter(x,sum,s=2,c=colours[n])
            plt.plot(x,sum,color=colours[n],label=f[:-4])
            plt.xlabel("Energy (nm)")
            plt.ylabel("$\epsilon$ (L mol$^{-1}$ cm$^{-1}$)")
            if args.sticks:
                for i in range(len(energies)):
                    plt.plot((energies[i],energies[i]),(0,stick_intensities[i]),colours[n])



    if args.r:
        plt.xlim(min(args.r),max(args.r))
    plt.legend()
    if args.save:
        plt.savefig(args.save+".pdf")
    plt.show()
