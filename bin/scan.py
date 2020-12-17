#! /usr/bin/env python3

import pandas as pd
import argparse

import read_directory
import read_orca
import read_tmole
import read_gaussian
import formatting
import plot_scan

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
parser.add_argument('--noenergy', help="Doesn't print the energy scan", action='store_true')
parser.add_argument('-t','--theodore', help="Print the TheoDORE descriptors", action='store_true')
parser.add_argument('-e','--extension',help="Determine the type of extension to look for",type=str,default='.out')
parser.add_argument('-s','--save',help="Determine the output file format",type=str,nargs='?',const='pdf')
parser.add_argument('--states',help="Number of states of each irrep to be plotted",type=int,default=1000)
parser.add_argument('-o','--output',help="Base name of the output file",type=str,nargs='?',default='scan')
parser.add_argument('--noshow',help="Stop from plotting the graph in the screen",nargs='?',type=bool,const=True,default=False)
parser.add_argument('-p','--program', help='Program that generated the results', type=str)
parser.add_argument('-pls','--plot_linestyle', help='Change line style of the different symmetry or multiplicity curves', type=str, nargs='+')
parameter = parser.add_mutually_exclusive_group(required=True)
parameter.add_argument('-b','--bond',help="Atom index to calculate the bond distanced.",type=int,nargs=2,metavar='ATOM')
parameter.add_argument('-a','--angle',help="Atom index to calculate the angle.", type=int,nargs=3,metavar='ATOM')
parameter.add_argument('-d','--dihedral',help="Atom index to calculate the dihedral angle.",type=int,nargs=4,metavar='ATOM')
parameter.add_argument('-g','--general', help="When the modification isn't symple", action='store_true')
args=parser.parse_args()

if args.program is None and (args.input[0].endswith('escf.out') or \
                             args.input[0].endswith('adc2.out') or \
                             args.input[0].endswith('ricc2.out')):
    args.program = 'turbomole'

elif args.program is None and args.input[0].endswith('.log'):
    args.program = 'gaussian'

elif args.program is None:
    args.program = 'orca'

###################      MAIN     ###################
if __name__=='__main__':

    paths =read_directory.find_outputs(args)

    if len(paths)==0:
        print("No output files found.")
        quit()

    if args.program == 'orca':
        print('Working with ORCA')
        energy, oscillator_str, charge_transfer = read_orca.energy(args, paths)
    
    elif args.program == 'turbomole':
        print('Working with TURBOMOLE')
        energy, oscillator_str, charge_transfer = read_tmole.energy(args, paths)

    elif args.program == 'gaussian':
        print('Working with GAUSSIAN')
        energy, oscillator_str, charge_transfer = read_gaussian.energy(args, paths)

    else:
        print("Only supported outputs from ORCA and TURBOMOLE, please select one of them with:")
        print("(-p|--program) (orca|turbomole)")
        exit()

    energy, oscillator_str = formatting.energies( args, energy, oscillator_str )

    if args.noenergy is False:
        energy.to_csv('energy.csv')
        plot_scan.energy(args, energy)

    if not oscillator_str.empty:
        oscillator_str.to_csv('oscillator_str.csv')
        plot_scan.oscillator_str(args, oscillator_str)

    if not charge_transfer.empty:
        charge_transfer.to_csv('charge_transfer.csv')
        plot_scan.charge_transfer(args, charge_transfer)
