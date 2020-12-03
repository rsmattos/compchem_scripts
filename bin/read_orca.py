#! /usr/bin/env python3

import numpy as np
import pandas as pd
from collections import defaultdict
from geometric_param import *

def energy(args, paths, output):
    energy = defaultdict(dict)

    tden_summ = []

    if args.descriptors:
        tden_summ_header = pd.read_csv(paths[0]+'/'+'tden_summ.txt', sep='\s+', nrows=0).columns.tolist() 
        tden_summ = pd.DataFrame(columns=tden_summ_header)

    for i in range(len(paths)):
        p=[]

        file=open(paths[i]+output,"r")
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
                energy[0][parameter] = float(lines[line].split()[5])

            # reads the excited states energies
            elif ( "STATE"+str(state) ) in lines[line].replace(" ",""):
                energy[state][parameter] = float(lines[line].split()[5])
                state += 1

        file.close()

        if args.descriptors:
            tmp = pd.read_csv(paths[i]+'/'+'tden_summ.txt', sep='\s+', header=0, names=tden_summ_header, skiprows=1, nrows=1)
            tmp['index'] = parameter
            tden_summ = tden_summ.append(tmp, ignore_index=True)

    return pd.DataFrame(energy)

    ###################      MAIN     ###################
if __name__=='__main__':
    paths=""
    output="unbv.out"

    #energy = energy(paths,output)

    print(energy(paths,output))
