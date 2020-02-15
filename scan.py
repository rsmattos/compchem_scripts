#! /usr/bin/env python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input",help="File with a list of outputs to be read", type=str, nargs='*')
args=parser.parse_args()

def read_output_files(outs):
    file_list=[]
    line=outs.readlines()
    for i in range(len(line)):
        file_list.append(line[i].rstrip())
    return file_list

def read_energies(outputs):
    geom=[]
    for i in range(len(outputs)):
        energies=[]
        file=open(outputs[i],"r")
        line=file.readlines()
        for j in range(len(line)):
            if "Total Energy " in line[j]:
                energies.append(float(line[j].split()[5]))
            elif "STATE " in line[j]:
                energies.append(float(line[j].split()[5]))
        geom.append(energies)
        file.close()

    state=[]
    for i in range(len(geom[0])):
        tmp=[]
        for j in range(len(geom)):
            tmp.append(geom[j][i])
        state.append(tmp)
    return state

def calc_energies(state):
    tot_en=[]
    for i in range(len(state[0])):
        tot_en.append(state[0][i])

    minimal=min(tot_en)
    for i in range(len(state[0])):
        state[0][i]=tot_en[i]-minimal
        for j in range(len(state)-1):
            state[j+1][i]=state[j+1][i]+state[j+1][i]

    return state

if __name__=='__main__':
    for n,f in enumerate(args.input):
        outs=open(f,"r")
        outputs=read_output_files(outs)
        outs.close()

        state=read_energies(outputs)

        state=calc_energies(state)

        import matplotlib.pyplot as plt

        x=list(range(len(state[0])))

        for i in range(len(state)):
            plt.plot(x,state[i],label=i)
        plt.show()
