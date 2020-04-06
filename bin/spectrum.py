#! /usr/bin/env python

import numpy as np
import math
from sys import argv
from os import system
import argparse
import matplotlib.pyplot as plt


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
parser.add_argument("-gnu",help="Plot a spectrum using gnuplot",action="store_true")
parser.add_argument("-prog",help="Specify the quantum chem program used",type=str)
parser.add_argument("-mpl",help="Plot a spectrum using matplotlib",action="store_true")
parser.add_argument("-sticks",help="Plot the stick spectrum",action="store_true")
parser.add_argument("-sd",help="Standard deviation (in eV)",default=0.4,type=float)
parser.add_argument("-rng",help="Min and max values for the spectrum (in nm)",nargs=2,type=int)
parser.add_argument("-save",help="Save spectrum with matplotlib",nargs='?',const="spectrum", type=str)
parser.add_argument("-dat",help="Save data as a text file",nargs='?',const="spectrum",type=str)
parser.add_argument("-csv",help="Save data as a csv file",nargs='?',const="spectrum",type=str)
args=parser.parse_args()

def find_program(line):
    for i in range(5):
        if "Gaussian" in line[i]:
            args.prog="gaussian"
            return
        if "* O   R   C   A *" in line[i]:
            args.prog="orca"
            return
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
                energies.append(float(line[j].split()[2]))
                os_strengths.append(float(line[j].split()[3]))
    return energies,os_strengths

#def abs_max(f,lamb_max,lamb):
#    a=1.3062974e8
#    b=f/(1e7/3099.6)
#    c=np.exp(-(((1/lamb-1/lamb_max)/(1/(1240/args.sd)))**2))
#    return a*b*c

def abs_max(f,lamb_max,lamb):
    # 1240 passes ards.sd from eV to nm-1
    return f*np.exp(-((1/lamb-1/lamb_max)*1240/args.sd)**2)

def out_data(xaxis,yaxis):
    with open(args.dat+".dat","w") as d:
        d.write("# index       Wavelenght  Oscilator_strenght\n")
        for i in range(len(xaxis)):
            d.write("{0:>7} {1:>16}   {2:>16}\n".format(i+1,xaxis[i],yaxis[i]))
        d.close()

def out_csv(xaxis,yaxis):
    with open(args.csv+".csv","w") as d:
        d.write("index,Wavelenght,Oscilator_strenght\n")
        for i in range(len(xaxis)):
            d.write("{0},{1},{2}\n".format(i+1,xaxis[i],yaxis[i]))
        d.close()

def gnu_plot(xaxis,yaxis):
    with open("data","w") as d:
        for i in range(len(xaxis)):
            d.write("{0} {1} {2}\n".format(i+1,xaxis[i],yaxis[i]))
        d.close()

    if args.sticks:
        with open("excit","w") as e:
            for i in range(len(energies)):
                e.write("{0} {1} {2}\n".format(i+1,energies[i],os_strengths[i]))
            e.close()

    with open("plot","w") as p:
        p.write("set xlabel \"Energy (nm)\"\nset ylabel \"Oscilator strenght\"\n")
        p.write("plot 'data' using 2:3 title '' with lines lt 1 lw 2,\\\n")
        if args.sticks:
            p.write("     'excit' using 2:3 title '' with impulse lt 2 lw 2,\\\n")
        p.close()

    system("gnuplot -persist plot")
    system("rm -f plot data")
    if args.sticks:
        system("rm -f excit")
    return

def mpl_plot(xaxis,yaxis):
    colours=["red","blue","green","orange","black","cyan","magenta"]
    plt.scatter(x,sum,s=2,c=colours[n])
    plt.plot(x,sum,color=colours[n],label=f[:-4])
    plt.xlabel("Energy (nm)")
    plt.ylabel("Oscilator strenght")
    if args.sticks:
        for i in range(len(energies)):
            plt.vlines(energies[i],0,os_strengths[i],color=colours[n])
    if args.rng:
        plt.xlim(min(args.rng),max(args.rng))
    plt.legend()

    return

if __name__=='__main__':
    for n,f in enumerate(args.input):
        infile=open(f,"r")
        line=infile.readlines()

        # find out the program that generated the output
        if not args.prog:
            find_program(line)

        # read the energies and oscilator strenghts from output file
        if args.prog=="orca":
            energies,os_strengths=read_orca(line)
        elif args.prog=="gaussian":
            energies,os_strengths=read_g09(line)
        else:
            print("Program not supported.")
            quit()

        infile.close()

        # set x axis
        if args.rng:
            x=np.linspace(max(args.rng),min(args.rng),1000)

        else:
            x=np.linspace(max(energies)+200,min(energies)-200,1000)

        # calculate y axis
        sum=[]
        for lamb in x:
            tot=0
            for i in range(len(energies)):
                tot+=abs_max(os_strengths[i],energies[i],lamb)
            sum.append(tot)

        # output formats

        if args.dat:
            out_data(x,sum)

        if args.csv:
            out_csv(x,sum)

        if args.gnu:
            gnu_plot(x,sum)

        if args.mpl or args.save:
            mpl_plot(x,sum)
    
    if args.save:
        plt.savefig(args.save+".pdf")
    if args.mpl:
        plt.show()
