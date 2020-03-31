#! /usr/bin/env python3

import numpy as np
from collections import defaultdict
import argparse
import os
import fnmatch

# TODO
# include option to print output as dat
# include general input for different variations of coordinates
# include simple gnuplot option
# improve matplotlib plotting

parser = argparse.ArgumentParser(description='''\
Script to extrat the energies from a scan calculation. \n \
When passing the atom index, the first atom given is connected to the second, which is connected to the third, and so on. \n \
Starting to count with the first atom in the coordinates list being 1.''')
parser.add_argument('input',help="Output files or directories containing them.",type=str,nargs='*',default='.')
parser.add_argument('-p','--paths_file',help="File with a list of outputs to be read.",type=str)
parser.add_argument('-e','--extension',help="Determine the type of extension to look for",type=str,default='.out')
parser.add_argument('-s','--save',help="Determine the output file format",type=str,nargs='?',const='pdf')
parser.add_argument('-o','--output',help="Base name of the output file",type=str,nargs='?',default='scan')
parser.add_argument('--noshow',help="Stop from plotting the graph in the screen",nargs='?',type=bool,const=True,default=False)
variation = parser.add_mutually_exclusive_group(required=True)
variation.add_argument('-b','--bond',help="Atom index to calculate the bond distanced.",type=int,nargs=2,metavar='ATOM')
variation.add_argument('-a','--angle',help="Atom index to calculate the angle.", type=int,nargs=3,metavar='ATOM')
variation.add_argument('-d','--dihedral',help="Atom index to calculate the dihedral angle.",type=int,nargs=4,metavar='ATOM')
args=parser.parse_args()

###################      CALCULATE PARAMETER VARIATION     ###################
def distance(p):
    return np.linalg.norm(p[0]-p[1])

def angle(p):
    x = p[0] - p[1]
    y = p[2] - p[1]

    # Making unitary vectors from arbitrary ones
    xu = x/np.linalg.norm(x)
    yu = y/np.linalg.norm(y)

    return np.degrees(np.arccos(np.dot(xu, yu)))

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

###################      READING FILES     ###################
# read file with a list of paths to the calculation outputs
def read_paths_file():
    file_list=[]

    outs=open(args.paths_file,"r")
    line=outs.readlines()

    print("The output files being used are:")
    for i in range(len(line)):
        print(line[i].rstrip())
        file_list.append(line[i].rstrip())

    outs.close()

    return file_list

# searches in the folder or passed arguments and creates a list
# of paths to the calculation outputs
def read_output_files():
    file_list=[]
    print("The output files being used are:")

    for inp in args.input:
        if os.path.isfile(inp):
            file_list.append(inp)

        elif os.path.isdir(inp):
            # if the path given is a directory, searches inside for output files
            for path,dir,file in os.walk(inp):
                for file_name in file:
                    if fnmatch.fnmatch(file_name, '*'+args.extension):
                        print(path+'/'+file_name)
                        file_list.append(path+'/'+file_name)

    return file_list

# read the calculation outputs in search for the parameter
# variation and the energies
def read_energies(outputs):
    energy = defaultdict(dict)

    for i in range(len(outputs)):
        p=[]

        file=open(outputs[i],"r")
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
                    parameter=float(distance(p))

                elif args.angle:
                    for atom in args.angle:
                        p.append(np.array([float(lines[line+1+atom].split()[1]),
                                           float(lines[line+1+atom].split()[2]),
                                           float(lines[line+1+atom].split()[3])]))
                    parameter=float(angle(p))

                elif args.dihedral:
                    for atom in args.dihedral:
                        p.append(np.array([float(lines[line+1+atom].split()[1]),
                                           float(lines[line+1+atom].split()[2]),
                                           float(lines[line+1+atom].split()[3])]))
                    parameter=float(dihedral(p))

            # reads the ground state energy
            elif "Total Energy " in lines[line]:
                energy[0][parameter] = float(lines[line].split()[5])

            # reads the excited states energies
            elif ( "STATE"+str(state) ) in lines[line].replace(" ",""):
                energy[state][parameter] = float(lines[line].split()[5])
                state += 1

        file.close()
    return energy

###################      MANIPULATING DATA     ###################
# the excited states energies are in relation to the ground state in
# its own geometry, for the plot is necessary to find the lowest gound
# state energy, set it to zero and calculate all energies in relation
# to the lowest ground state
def calc_energies_dic(state):

    minimal=state[0][min(state[0], key=state[0].get)]

    # dislocates all geometris ground states
    for i in state[0]:
        state[0][i]=state[0][i]-minimal

    # dislocates all excited state energies
    for i in state:
        if not i==0:
            for j in state[i]:
                state[i][j]=state[i][j]+state[0][j]
    return state

###################      OUTPUTTING     ###################
def plot_matplot(state):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FormatStrFormatter

    plt.style.use("seaborn-deep")

    fig , ax = plt.subplots()

    # Axis ticks
    ax.xaxis.set_tick_params(top=False, direction='out', width=1)
    ax.xaxis.set_tick_params(bottom=True, direction='in', width=1)
    ax.yaxis.set_tick_params(right=False, direction='in', width=1)
    ax.yaxis.set_tick_params(bottom=True, direction='in', width=1)

    plt.rc('font', family='sans-serif')
    plt.tick_params(labelsize=10)

    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))

    # Axis labels
    ax.set_ylabel(r"$\Delta$E (eV)", fontsize=10)

    if args.bond:
        ax.set_xlabel("Distance (Angstron)", fontsize=10)
    elif args.angle:
        ax.set_xlabel("Angle (Degree)", fontsize=10)
    elif args.dihedral:
        ax.set_xlabel("Dihedral angle (Degree)", fontsize=10)

    # Image size
    fig.set_size_inches(5.0, 5.0)

    # Plotting
    for i in state:
        ax.plot(*zip(*sorted(state[i].items())),label="S"+str(i),marker='o',markersize=3)

    # Legend
    box = ax.get_position()
    ax.set_position([box.x0 + box.width*0.05, box.y0 + box.height * 0.1, box.width, box.height * 0.95])

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), fancybox=True, shadow=True, ncol=6, fontsize='small')

    # Saving
    if args.save:
        diction = plt.gcf().canvas.get_supported_filetypes()

        if(args.save in diction):
            fig.savefig(args.output+'.'+args.save)
        else:
            print("Graphic plot type not supported, the available formats are:")
            for key in diction:
                print (key, " => ", diction[key])

    if args.noshow:
        return

    plt.show()


###################      MAIN     ###################
if __name__=='__main__':
    if args.paths_file:
        outputs=read_paths_file()
    else:            
        outputs=read_output_files()

    if len(outputs)==0:
        print("No output files found.")
        quit()

    state=read_energies(outputs)
    print()
    
    state=calc_energies_dic(state)

    plot_matplot(state)