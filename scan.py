#! /usr/bin/env python3

import numpy as np
from collections import defaultdict
import argparse

# TODO
# include arguments for reading bond distance, angle and dihedral variation, inputing atom index
# include reading options, automatic iterates through folder, read a file with paths, read path input with default value, read a single output
# include option to print output as dat
# include general input for different variations of coordinates
# include simple gnuplot option
# improve matplotlib plotting

parser = argparse.ArgumentParser(description="Script to extrat the energies from a scan calculation.")
parser.add_argument('input',help="File with a list of outputs to be read", type=str, nargs='*')
variation = parser.add_mutually_exclusive_group(required=True)
variation.add_argument('-b','--bond_distance',help="Atom index to calculate the bond distanced. Starting to count with the first atom in the coordinates list being 1.", type=int, nargs=2)
variation.add_argument('-d','--dihedral',help="Atom index to calculate the dihedral angle,being that the first atom is bonded to the second and so on. Starting to count with the first atom in the coordinates list being 1.", type=int, nargs=4)
args=parser.parse_args()

# functino to calculate the dihedral angle from cartesian coordinates, taken from
# https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
def dihedral(p):
    """Praxeolitic formula 1 sqrt, 1 cross product"""
    p0 = p[0]
    p1 = p[1]
    p2 = p[2]
    p3 = p[3]

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))

def read_output_files(outs):
    file_list=[]
    line=outs.readlines()
    for i in range(len(line)):
        file_list.append(line[i].rstrip())
    return file_list

def read_energies(outputs):
    energy = defaultdict(dict)
    for i in range(len(outputs)):
        p=[]
        file=open(outputs[i],"r")
        line=file.readlines()
        for j in range(len(line)):
            if "CARTESIAN COORDINATES (ANGSTROEM)" in line[j]:
                if args.bond_distance:
                    for k in args.bond_distance:
                        p.append(np.array([float(line[j+1+k].split()[1]),
                                           float(line[j+1+k].split()[2]),
                                           float(line[j+1+k].split()[3])]))
                    variable=float(np.linalg.norm(p[0]-p[1]))

                elif args.dihedral:
                    for k in args.dihedral:
                        p.append(np.array([float(line[j+1+k].split()[1]),
                                           float(line[j+1+k].split()[2]),
                                           float(line[j+1+k].split()[3])]))
                    variable=int(dihedral(p))

            elif "Total Energy " in line[j]:
                energy[0][variable] = float(line[j].split()[5])
            elif "STATE " in line[j]:
                energy[int(line[j].split()[1].split(':')[0])][variable] = float(line[j].split()[5])

        file.close()
    return energy

def calc_energies_dic(state):

    minimal=state[0][min(state[0], key=state[0].get)]

    for i in state[0]:
        state[0][i]=state[0][i]-minimal

    for i in state:
        if not i==0:
            for j in state[i]:
                state[i][j]=state[i][j]+state[0][j]
    return state

if __name__=='__main__':
    for n,f in enumerate(args.input):

        outs=open(f,"r")
        outputs=read_output_files(outs)
        outs.close()

        state=read_energies(outputs)

        state=calc_energies_dic(state)

        import matplotlib.pyplot as plt

        for i in state:
            plt.plot(*zip(*sorted(state[i].items())))
        plt.show()
