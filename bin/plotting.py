#! /usr/bin/env python3

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

def energy(args, state):

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
    elif args.general:
        ax.set_xlabel("Steps", fontsize=10)

    # Image size
    fig.set_size_inches(5.0, 5.0)

    # Plotting
    for S in state.columns:
        ax.plot(state[S],label="S"+str(S),marker='o',markersize=3)

    # Legend
    box = ax.get_position()
    ax.set_position([box.x0 + box.width*0.05, box.y0 + box.height * 0.1, box.width, box.height * 0.95])

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), fancybox=True, shadow=True, ncol=6, fontsize='small')

    # Saving
    if args.save:
        diction = plt.gcf().canvas.get_supported_filetypes()

        if(args.save in diction):
            fig.savefig(args.output+'.energy.'+args.save)
        else:
            print("Graphic plot type not supported, the available formats are:")
            for key in diction:
                print (key, " => ", diction[key])

    if args.noshow:
        return

    plt.show()

    return

def descriptors(args, tden_summ):

    plt.style.use("seaborn-deep")

    fig , axs = plt.subplots(2, 1)

    # Axis ticks
    axs[0].xaxis.set_tick_params(top=False, direction='out', width=1)
    axs[0].xaxis.set_tick_params(bottom=True, direction='in', width=1)
    axs[0].yaxis.set_tick_params(right=False, direction='in', width=1)
    axs[0].yaxis.set_tick_params(bottom=True, direction='in', width=1)

    plt.rc('font', family='sans-serif')
    plt.tick_params(labelsize=10)

    axs[0].yaxis.set_major_formatter(FormatStrFormatter("%.2f"))

    # Axis labels
#    axs[0].set_ylabel(r"$\Delta$E (eV)", fontsize=10)

    if args.bond:
        axs[0].set_xlabel("Distance (Angstron)", fontsize=10)
    elif args.angle:
        axs[0].set_xlabel("Angle (Degree)", fontsize=10)
    elif args.dihedral:
        axs[0].set_xlabel("Dihedral angle (Degree)", fontsize=10)
    elif args.general:
        axs[0].set_xlabel("Steps", fontsize=10)

    # Image size
    fig.set_size_inches(5.0, 5.0)

    # Plotting
    axs[1].plot(tden_summ['CT'],label='CT',marker='o',markersize=3)
    axs[1].plot(tden_summ['f'],label='f',marker='o',markersize=3)
    axs[0].plot(tden_summ['POS'],label='POS',marker='o',markersize=3)
    axs[0].plot(tden_summ['PR'],label='PR',marker='o',markersize=3)
    axs[0].plot(tden_summ['PRNTO'],label='$PR_{NTO}$',marker='o',markersize=3)

    # Legend
    box = axs[0].get_position()
    axs[0].set_position([box.x0 + box.width*0.05, box.y0 + box.height * 0.1, box.width, box.height * 0.95])

    axs[0].legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), fancybox=True, shadow=True, ncol=3, fontsize='small')

    box = axs[1].get_position()
    axs[1].set_position([box.x0 + box.width*0.05, box.y0 + box.height * 0.1, box.width, box.height * 0.95])

    axs[1].legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), fancybox=True, shadow=True, ncol=2, fontsize='small')

    # Saving
    if args.save:
        diction = plt.gcf().canvas.get_supported_filetypes()

        if(args.save in diction):
            fig.savefig(args.output+'.descriptors.'+args.save)
        else:
            print("Graphic plot type not supported, the available formats are:")
            for key in diction:
                print (key, " => ", diction[key])

    if args.noshow:
        return

    plt.show()

    return
