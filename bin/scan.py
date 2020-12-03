#! /usr/bin/env python3

import argparse

import read_directory
import read_orca
import formatting
import plotting

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
parser.add_argument('--noenergy', help="Doesn't print the energy scan", action='store_false')
parser.add_argument('--descriptors', help="Print the TheoDORE descriptors", action='store_true')
parser.add_argument('-e','--extension',help="Determine the type of extension to look for",type=str,default='.out')
parser.add_argument('-s','--save',help="Determine the output file format",type=str,nargs='?',const='pdf')
parser.add_argument('-o','--output',help="Base name of the output file",type=str,nargs='?',default='scan')
parser.add_argument('--noshow',help="Stop from plotting the graph in the screen",nargs='?',type=bool,const=True,default=False)
variation = parser.add_mutually_exclusive_group(required=True)
variation.add_argument('-b','--bond',help="Atom index to calculate the bond distanced.",type=int,nargs=2,metavar='ATOM')
variation.add_argument('-a','--angle',help="Atom index to calculate the angle.", type=int,nargs=3,metavar='ATOM')
variation.add_argument('-d','--dihedral',help="Atom index to calculate the dihedral angle.",type=int,nargs=4,metavar='ATOM')
variation.add_argument('-g','--general', help="When the modification isn't symple", action='store_true')
args=parser.parse_args()

###################      MAIN     ###################
if __name__=='__main__':
    paths, file =read_directory.find_output_files(args)

    if len(paths)==0:
        print("No output files found.")
        quit()

    energy = read_orca.energy(args, paths, file)

    if args.noenergy:
        energy = formatting.energies( args, energy )

        energy.to_csv('energy.csv')
        
        plotting.energy(args, energy)

    if args.descriptors:
        tden_summ.to_csv('descriptors.csv')

        plotting.descriptors(args, tden_summ )
