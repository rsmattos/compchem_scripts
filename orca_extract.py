#! /usr/bin/env python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input",help="Log file of Gaussian 09 TD job", type=str, nargs='*')
parser.add_argument("-dat",help="Save data as text file",default="data",type=str)
parser.add_argument("-tex",help="Save data as latex formated table",default="latex_table",type=str)
args=parser.parse_args()

def read_abs_spectra(line):
    energies_cm=[]
    energies_nm=[]
    os_strengths=[]
    for i in range(len(line)):
        if " nroots " in line[i]:
            nroots=int(line[i].split()[3])
        if " ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" in line[i]:
            for j in range(i+5,i+5+nroots):
                energies_cm.append(float(line[j].split()[1]))
                energies_nm.append(float(line[j].split()[2]))
                os_strengths.append(float(line[j].split()[3]))
    return energies_cm,energies_nm,os_strengths

def out_dat():
    with open(args.dat+".dat","w") as d:
        d.write("# State   Wavelenght(nm) Oscilator_strength\n")
        for i in range(len(energies_nm)):
            d.write("{0:>6} {1:>16} {2:>16}\n".format(i+1,energies_nm[i],os_strengths[i]))
        d.close()

#def out_tex_table:

if __name__=='__main__':
    for n,f in enumerate(args.input):
        infile=open(f,"r")
        line=infile.readlines()

        # read vertical excitation energies
        energies_cm,energies_nm,os_strengths=read_abs_spectra(line)

        infile.close()

        if args.dat:
            out_dat()

#        if args.tex:
#            out_tex_table()