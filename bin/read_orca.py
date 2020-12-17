#! /usr/bin/env python3

##############################################################################
#  Reads the geometry and energies from the orca output, and outputs as a 
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

    charge_transfer = pd.DataFrame()
    oscillator_str = pd.DataFrame()

    if args.theodore:
        tden_summ_header = pd.read_csv(paths[0][:-len(paths[0].split('/')[-1])]+'tden_summ.txt', sep='\s+', nrows=0).columns.tolist() 

    for i in range(len(paths)):
        p=[]

        file=open(paths[i],"r")
        lines=file.readlines()

        state = 1

        for line in range(len(lines)):
            # reads the parameter being evaluated
            if "CARTESIAN COORDINATES (ANGSTROEM)" in lines[line]:
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
            elif "Total Energy " in lines[line]:
                energy['GS'][parameter] = float(lines[line].split()[5])

            # reads the excited states energies
            elif ( "STATE"+str(state) ) in lines[line].replace(" ",""):
                energy[state][parameter] = float(lines[line].split()[5])
                state += 1

        file.close()

        if args.theodore:
            tmp = pd.read_csv(paths[i][:-len(paths[i].split('/')[-1])]+'tden_summ.txt', sep='\s+', header=None, \
                              names=tden_summ_header, index_col='state', skiprows=2)
            oscillator_str[parameter] = tmp['f']
            charge_transfer[parameter] = tmp['CT']

    return pd.DataFrame(energy), oscillator_str.transpose(), charge_transfer.transpose()
