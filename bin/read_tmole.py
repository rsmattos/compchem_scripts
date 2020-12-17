#! /usr/bin/env python3

##############################################################################
#  Reads the geometry and energies from the utrbomole output, and outputs as a 
#  Dataframe
#  Input: enters the path and name of the files
#  Output: Dataframe with energies, column header has the states, vertical
#          index has the geometrica parameters
##############################################################################

import numpy as np
import pandas as pd
from collections import defaultdict
from geometric_param import *

def escf(args, paths):
    
    energy = defaultdict(dict)
    oscillator_str = defaultdict(dict)

    charge_transfer = pd.DataFrame()

    if args.theodore:
        tden_summ_header = pd.read_csv(paths[0][:-len(paths[0].split('/')[-1])]+'tden_summ.txt', sep='\s+', nrows=0).columns.tolist() 

    for i in range(len(paths)):
        p=[]

        file=open(paths[i],"r")
        lines=file.readlines()

        for line in range(len(lines)):
            # reads the parameter being evaluated
            if "| Atomic coordinate, charge and isotop information |" in lines[line]:
                if args.bond:
                    for atom in args.bond:
                        p.append(np.array([float(lines[line+3+atom].split()[0]),
                                           float(lines[line+3+atom].split()[1]),
                                           float(lines[line+3+atom].split()[2])]))
                    parameter=round(float(distance(p)))

                elif args.angle:
                    for atom in args.angle:
                        p.append(np.array([float(lines[line+3+atom].split()[0]),
                                           float(lines[line+3+atom].split()[1]),
                                           float(lines[line+3+atom].split()[2])]))
                    parameter=round(float(angle(p)))

                elif args.dihedral:
                    for atom in args.dihedral:
                        p.append(np.array([float(lines[line+3+atom].split()[0]),
                                           float(lines[line+3+atom].split()[1]),
                                           float(lines[line+3+atom].split()[2])]))
                    parameter=round(float(dihedral(p)))

                elif args.general:
                  	parameter=i+1

            # reads the ground state energy
            elif "              Ground state" in lines[line]:
                energy['GS'][parameter] = float(lines[line+3].split()[2])*27.2114

            elif "              I R R E P" in lines[line]:
                irrep = lines[line].split()[5]
                state = 1

            # reads the excited states energies
            elif ( "Excitation energy / eV:" ) in lines[line]:
                if lines[line-7].split()[1] == "singlet":
                    mult = "1"
                elif lines[line-7].split()[1] == "triplet":
                    mult = "3"

                if state > args.states:
                    continue
                energy[str(state)+'$^{'+mult+'}$'+irrep.upper()][parameter] = float(lines[line].split()[4])
                oscillator_str[str(state)+'$^{'+mult+'}$'+irrep.upper()][parameter] =  float(lines[line+9].split()[2])
                state += 1

        file.close()

        if args.theodore:
            tmp = pd.read_csv(paths[i][:-len(paths[i].split('/')[-1])]+'tden_summ.txt', sep='\s+', header=None, \
                              names=tden_summ_header, index_col='state', skiprows=2)
            charge_transfer[parameter] = tmp['CT']

    return pd.DataFrame(energy), pd.DataFrame(oscillator_str), charge_transfer.transpose()

def adc2(args, paths):
    
    energy = defaultdict(dict)

    charge_transfer = pd.DataFrame()
    oscillator_str = pd.DataFrame()

    if args.theodore:
        tden_summ_header = pd.read_csv(paths[0][:-len(paths[0].split('/')[-1])]+'tden_summ.txt', sep='\s+', nrows=0).columns.tolist() 

    for i in range(len(paths)):
        p=[]

        file=open(paths[i],"r")
        lines=file.readlines()

        for line in range(len(lines)):
            # reads the parameter being evaluated
            if "| Atomic coordinate, charge and isotop information |" in lines[line]:
                if args.bond:
                    for atom in args.bond:
                        p.append(np.array([float(lines[line+3+atom].split()[0]),
                                           float(lines[line+3+atom].split()[1]),
                                           float(lines[line+3+atom].split()[2])]))
                    parameter=round(float(distance(p)))

                elif args.angle:
                    for atom in args.angle:
                        p.append(np.array([float(lines[line+3+atom].split()[0]),
                                           float(lines[line+3+atom].split()[1]),
                                           float(lines[line+3+atom].split()[2])]))
                    parameter=round(float(angle(p)))

                elif args.dihedral:
                    for atom in args.dihedral:
                        p.append(np.array([float(lines[line+3+atom].split()[0]),
                                           float(lines[line+3+atom].split()[1]),
                                           float(lines[line+3+atom].split()[2])]))
                    parameter=round(float(dihedral(p)))

                elif args.general:
                  	parameter=i+1

            # reads the ground state energy
            elif "            Total Energy    :" in lines[line]:
                energy['GS'][parameter] = float(lines[line].split()[3])*27.2114

            # reads the excited states energies
            elif ( "Excited state reached by transition:" ) in lines[line]:
                state = int(lines[line+1].split()[4])
                irrep = lines[line+1].split()[5]
                mult = lines[line+1].split()[6]
                if state > args.states:
                    continue
                energy[str(state)+'$^{'+mult+'}$'+irrep.upper()][parameter] = float(lines[line+2].split()[5])

        file.close()

        if args.theodore:
            tmp = pd.read_csv(paths[i][:-len(paths[i].split('/')[-1])]+'tden_summ.txt', sep='\s+', header=None, \
                              names=tden_summ_header, index_col='state', skiprows=2)
            oscillator_str[parameter] = tmp['f']
            charge_transfer[parameter] = tmp['CT']

    return pd.DataFrame(energy), oscillator_str.transpose(), charge_transfer.transpose()

def energy(args, paths):
    if not paths[0].endswith('escf.out') and \
       not paths[0].endswith('adc2.out') and \
       not paths[0].endswith('ricc2.out'):
       print('Please specify the name of the file with the energies.')
       print('Acceptable files are escf.out, adc2.out or ricc2.out')
       exit()
    
    if paths[0].endswith('escf.out'):
        return escf(args, paths)

    if paths[0].endswith('adc2.out'):
        return adc2(args, paths)