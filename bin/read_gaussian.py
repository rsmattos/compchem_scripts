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

def energy(args, paths):
    
    energy = defaultdict(dict)
    oscillator_str = defaultdict(dict)

    if args.theodore:
        charge_transfer = defaultdict(dict)
        tden_summ_header = pd.read_csv(paths[0][:-len(paths[0].split('/')[-1])]+'tden_summ.txt', sep='\s+', nrows=0).columns.tolist() 

    for i in range(len(paths)):
        p=[]
        state = dict()

        file=open(paths[i],"r")
        lines=file.readlines()

        for line in range(len(lines)):
            # reads the parameter being evaluated
            if "Symbolic Z-matrix:" in lines[line]:
                if args.bond:
                    for atom in args.bond:
                        p.append(np.array([float(lines[line+1+atom].split()[1]),
                                           float(lines[line+1+atom].split()[2]),
                                           float(lines[line+1+atom].split()[3])]))
                    parameter=round(float(distance(p)))

                elif args.angle:
                    for atom in args.angle:
                        p.append(np.array([float(lines[line+1+atom].split()[1]),
                                           float(lines[line+1+atom].split()[2]),
                                           float(lines[line+1+atom].split()[3])]))
                    parameter=round(float(angle(p)))

                elif args.dihedral:
                    for atom in args.dihedral:
                        p.append(np.array([float(lines[line+1+atom].split()[1]),
                                           float(lines[line+1+atom].split()[2]),
                                           float(lines[line+1+atom].split()[3])]))
                    parameter=round(float(dihedral(p)))

                elif args.general:
                  	parameter=i+1


            # reads the ground state energy
            elif "SCF Done: " in lines[line]:
                energy['GS'][parameter] = float(lines[line].split()[4])*27.2114

            # reads the excited states energies
            elif " Excited State " in lines[line]:
                tmp = lines[line].replace("-", " ")
                irrep = tmp.split()[4]
                if irrep in state:
                    state[irrep] += 1
                else:
                    state[irrep] = 1

                if tmp.split()[3] == "Singlet":
                    mult = "1"
                elif tmp.split()[3] == "Triplet":
                    mult = "3"
                if state[irrep] > args.states:
                    continue
                energy[str(state[irrep])+'$^{'+mult+'}$'+irrep][parameter] = float(tmp.split()[5])
                oscillator_str[str(state[irrep])+'$^{'+mult+'}$'+irrep][parameter] =  float(tmp.split()[9][2:])

            elif "This state for optimization and/or second-order correction." in lines[line]:
                corrected_state = state[irrep]
                corrected_mult = mult
                corrected_irrep = irrep

            elif " Corrected transition energy" in lines[line]:
                energy[str(corrected_state)+'$^{'+corrected_mult+'}$'+corrected_irrep][parameter] = float(lines[line].split()[4])

        file.close()

        if args.theodore:
            tmp = pd.read_csv(paths[i][:-len(paths[i].split('/')[-1])]+'tden_summ.txt', sep='\s+', header=None, \
                              names=tden_summ_header, skiprows=2)
            state = dict()
            for index, row in tmp.iterrows():
                if 'Sing' in row['state']:
                    irrep = row['state'].split('Sing')[1]
                    if irrep in state:
                        state[irrep] += 1
                    else:
                        state[irrep] = 1

                    charge_transfer[str(state[irrep])+'$^{1}$'+irrep][parameter] = row['CT']

                elif 'Trip' in row['state']:
                    irrep = row['state'].split('Trip')[1]
                    if irrep in state:
                        state[irrep] += 1
                    else:
                        state[irrep] = 1

                    charge_transfer[str(state[irrep])+'$^{3}$'+irrep][parameter] = row['CT']

    return pd.DataFrame(energy), pd.DataFrame(oscillator_str), pd.DataFrame(charge_transfer)
