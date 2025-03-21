"""
summary_plots.py
10 January 2025

ENERGY, RMSD, RG, RMSF

RAFE, COORDS (Pending)

Secondary structure information (plotting residue-or-secondary structure wise deviations/dynamics in the trajectory)

Energy plots the energy fluctuations (potential, kinetic and total)
RAFE plots the relative average force errors (particle wise or total system)
COORDS plots the pearson correlation coefficient between various simulated trajectories
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.ndimage import gaussian_filter1d
import pandas as pd
import numpy as np
import sys
import os

simtime = 10

def get_pdb_list(simtime):
    """
    Gets the list of available pdb paths and corresponding legends according to the simulation time
    Change according to simulations conducted, file paths, and preferred sigma values for Gaussian smoothening
    """
    if simtime == 1:
        pdblist = ['../previous_simfiles/1xda_sp_merged_chain', '../previous_simfiles/1xda_dp_merged_chain', '../gromacs_sp', '../gromacs_dp', '../Calligo_sp', '../Calligo_dp_merged']
        legends = ['Gromacs_SP_R1', 'Gromacs_DP_R1', 'Gromacs_SP_R2', 'Gromacs_DP_R2', 'Calligo_SP','Calligo_DP']
        sigma = 4
    elif simtime == 10:
        pdblist = ['../previous_simfiles/1xda_sp_merged_chain', '../previous_simfiles/1xda_dp_merged_chain', '../gromacs_sp', '../gromacs_dp']
        legends = ['Gromacs_SP_R1', 'Gromacs_DP_R1', 'Gromacs_SP_R2', 'Gromacs_DP_R2']
        sigma = 20
    elif simtime == 100:
        pdblist = ['../gromacs_sp', '../gromacs_dp']
        legends = ['Gromacs_SP', 'Gromacs_DP']
        sigma = 200

    return pdblist, legends, sigma

# Plotting parameters and y-axis labels 
plt.rcParams['font.size'] = 16
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plot_yaxes={'rmsd':'RMSD (nm)', 'rmsf': 'RMSF (nm)', 'rg':'Radius of gyration (nm)', 'potential':'Potential energy (kJ/mol)', 'total_energy': 'Total energy (kJ/mol)', 'kinetic':'Kinetic energy (kJ/mol)','density': 'Density (kg/m^3)', 'temperature': 'Temperature (K)', 'pressure': 'Pressure (bar)', 'distance':"Distance (nm)"}

# Helper functions

def get_cols(prop):
    if prop in ['rmsd', 'potential', 'total_energy', 'pressure', 'density', 'temperature', 'kinetic', 'rmsf', 'distance']:
        cols = ['xaxis', 'property']
    elif prop == "rg":
        cols= ['xaxis', 'property', 'propx', 'propy', 'propz']
    return cols

def skiprows(prop):
    if prop == "rmsd":
        skip = 18
    elif prop == "rmsf":
        skip = 17
    elif prop == "rg":
        skip = 27
    elif prop == "ss":
        skip = 33
    elif prop in ["potential", "kinetic", "total_energy", "pressure", "temperature", "density", "distance"]:
        skip = 24
    return skip

def label_resids(file, chain):
    """
    Returns a list of residue names for a chain
    """
    chaindict = {'_chainA': 0, '_chainB': 1}
    with open(file, 'r') as file:
        lines = file.readlines()
        row = lines[chaindict[chain]].strip()
        return list(row)

def get_ranges(chain):
    if chain == "_chainA":
        return [1,9,12,18]
    elif chain == "_chainB":
        return [1,20]
    
# Plotting functions

def plot_energies(plotfile, simtime):
    """
    Plots energies averaged over the entire system (simulation quality analysis)
    """
    plt.figure(figsize=(13,8))
    prop = ['kinetic', 'potential', 'total_energy']

    pdblist, legends = get_pdb_list(simtime)

    for p in prop:
        cols = get_cols(p)
        for pdbnum in range(len(pdblist)):
            filename = pdblist[pdbnum]+"/"+str(simtime)+"_ns/results_"+str(simtime)+"/"+p+".xvg"
            resultsdir ="plots/"
            data = pd.read_csv(filename, sep="\s+", names=cols, skiprows=skiprows(p))

            xaxis = data['xaxis']
            yaxis = data['property']

            plt.plot(xaxis, yaxis, linewidth=2)

            plt.xlabel('Time (ps)', fontsize=20)
            plt.ylabel('Energies (kJ/mol)', fontsize=20)
            plt.grid('both')

    #     legend_handles = [
    # mpatches.Patch(color='blue', label='Gromacs_SP'),
    # mpatches.Patch(color='green', label='Gromacs_DP'),
    # mpatches.Patch(color='black', label='Calligo_SP'),
    # mpatches.Patch(color='orange', label='Calligo_DP')
    # ]
        # plt.legend(handles=legend_handles)
    plt.savefig(resultsdir+plotfile)
    plt.close()

def plot_thermo(prop, plotfile, simtime):
    """
    Plots thermodynamic properties of the entire system
    """
    cols = get_cols(prop)
    plt.figure(figsize=(12,8))

    pdblist, legends = get_pdb_list(simtime)
        
    for pdbnum in range(len(pdblist)):
        filename = pdblist[pdbnum]+"/"+str(simtime)+"_ns/results_"+str(simtime)+"/"+prop+".xvg"
        resultsdir = "plots/"
        data = pd.read_csv(filename, sep="\s+", names=cols, skiprows=skiprows(prop))

        xaxis = data['xaxis']
        yaxis = data['property']

        plt.plot(xaxis, yaxis, linewidth=2,label=legends[pdbnum])
        plt.xlabel('Time (ps)', fontsize=20)
        plt.ylabel(plot_yaxes[prop], fontsize=20)
        plt.grid('both')
        plt.legend()
    plt.savefig(resultsdir+plotfile)
    plt.close()

def plot_structs(prop, simtime, smooth):
    """
    Plots structural properties of the entire system
    """
    cols = get_cols(prop)
    suffixes = ['_chainA', '_chainB']
    pdblist, legends, sigma = get_pdb_list(simtime)

    for suffix in suffixes:
        plt.figure(figsize=(14,8))
        for pdbnum in range(len(pdblist)):
            filename = pdblist[pdbnum]+"/"+str(simtime)+"_ns/results_"+str(simtime)+"/"+prop + suffix +".xvg"
            resultsdir = "plots/smooth/"
            data = pd.read_csv(filename, sep="\s+", names=cols, skiprows=skiprows(prop))

            xaxis = data['xaxis']
            yaxis = data['property']
            y_smooth = gaussian_filter1d(yaxis, sigma)

            if prop == 'rmsf':
                helix_ranges = get_ranges(suffix)
                print(helix_ranges)
                plt.plot(xaxis, yaxis,linewidth=2, label=legends[pdbnum], marker="o")
                plt.xlabel('Residues', fontsize=20)
                plt.xlim([min(xaxis), max(xaxis)])
                plt.xticks(np.arange(0, max(xaxis), 2))
                plt.ylim([0,0.4])
                resid_codes = label_resids("sequence.txt", suffix)
                plt.xticks(ticks=xaxis, labels=[f'{i}\n{name}' for i, name in zip(xaxis, resid_codes)])
                for i in range(0, len(helix_ranges), 2):
                    if len(helix_ranges) != 0:
                        plt.fill_between(xaxis, yaxis, where=(xaxis >= helix_ranges[i]) & (xaxis <= helix_ranges[i+1]), color='red', alpha=0.1)
                if suffix == "_chainB":
                    plt.fill_between(xaxis, yaxis, where=(xaxis >= 21) & (xaxis <= 23), color='green', alpha=0.1)  
                    plt.fill_between(xaxis, yaxis, where=(xaxis >= 24) & (xaxis <= 26), color='blue', alpha=0.1)  
                plt.fill_between(xaxis, yaxis, color='grey', alpha=0.1)
       

            else:
                if smooth:
                    plt.plot(xaxis, y_smooth, linewidth=2, label=legends[pdbnum])
                else:
                    plt.plot(xaxis, yaxis, linewidth=2, label=legends[pdbnum])
                plt.xlabel('Time (ps)',fontsize=20)
                print("Mean "+prop+" ("+pdblist[pdbnum]+"): "+str(np.mean(yaxis)))
                print("Std_dev "+prop+" ("+pdblist[pdbnum]+"): "+str(np.std(yaxis)))
                # print("RMSD:"+prop+" ("+pdblist[pdbnum]+"): "+str(np.sqrt(((yaxis-np.ones(len(yaxis))*np.mean(yaxis)) ** 2).mean())))
            
            if prop == 'rg':
                plt.ylim([min(yaxis) - 0.05, max(yaxis) + 0.05])

            plt.ylabel(plot_yaxes[prop], fontsize=20)

            plt.grid('both')
            plt.legend()

        plotfile=prop+suffix+".png"

        plt.savefig(resultsdir+plotfile)
        plt.close()

def normalise_and_compare(prop,simtime):
    """
    Mean centers Rg values and plots them for different systems
    """
    cols = get_cols(prop)
    suffixes = ['_chainA', '_chainB']
    pdblist, legends = get_pdb_list(simtime)

    for suffix in suffixes:
        plt.figure(figsize=(13,8))

        for pdbnum in range(len(pdblist)):
            filename =pdblist[pdbnum]+"/"+str(simtime)+"_ns/results_"+str(simtime)+"/"+ prop + suffix +".xvg"
            resultsdir = "plots/"
            data = pd.read_csv(filename, sep="\s+", names=cols, skiprows=skiprows(prop))

            xaxis = data['xaxis']
            yaxis = data['property']

            mean = np.mean(yaxis)
            mean_centered = [yaxis[i] - mean for i in range(len(yaxis))]
            # normalised = stats.zscore(np.array(yaxis))
            plt.plot(xaxis, mean_centered, linewidth=2, label=legends[pdbnum])
            plt.xlabel('Time (ps)',fontsize=20)
            plt.ylabel(plot_yaxes[prop], fontsize=20)
            plt.ylim([min(mean_centered) - 0.05, max(mean_centered) + 0.05])
            plt.grid('both')
            plt.legend()

        plotfile=prop+suffix+"_mean.png"

        plt.savefig(resultsdir+plotfile)
        plt.close()
           
def get_statistics(prop,simtime):
    """Get average, RMSD, for thermo and structural properties"""

    pdblist, legends = get_pdb_list(simtime)
    suffixes = ['_chainA', '_chainB']
    for p in prop:
        if prop in ['rmsd', 'rmsf', 'rg']:
            for suffix in suffixes:
                cols = get_cols(p)
                for pdbnum in range(len(pdblist)):
                    filename = pdblist[pdbnum]+"/"+str(simtime)+"_ns/results_"+str(simtime)+"/"+p+suffix+".xvg"
                    data = pd.read_csv(filename, sep="\s+", names=cols, skiprows=skiprows(p))
                    print(legends[pdbnum])
                    xaxis = data['xaxis']
                    yaxis = data['property']

                    mean = np.mean(yaxis)
                    rmsd = np.sqrt(np.mean((np.ones(len(yaxis))*mean - yaxis) ** 2))
                    print(p+suffix,":\t",mean,"\t",rmsd)
        else:
            cols = get_cols(p)
            for pdbnum in range(len(pdblist)):
                filename = pdblist[pdbnum]+"/"+str(simtime)+"_ns/results_"+str(simtime)+"/"+p+".xvg"
                data = pd.read_csv(filename, sep="\s+", names=cols, skiprows=skiprows(p))
                print(legends[pdbnum])
                xaxis = data['xaxis']
                yaxis = data['property']

                mean = np.mean(yaxis)
                rmsd = np.sqrt(np.mean((np.ones(len(yaxis))*mean - yaxis) ** 2))
                print(p,":\t",mean,"\t",rmsd)

def plot_ss(chain, simtime):

    plt.figure(figsize=(13,11))
    pdblist, legends = get_pdb_list(simtime)

    for pdbnum in range(len(pdblist)):
        chainpath=pdblist[pdbnum]+"/"+str(simtime)+"_ns/results_"+str(simtime)+"/ss_data_"+chain+".xvg"
        resultsdir = "plots/"
        chaindata=pd.read_csv(chainpath, names=['time', 'loops', 'breaks', 'bends', 'turns', 'pp_helices', 'pi_helices', '3_10_helices', 'beta_strands', 'beta_bridges', 'alpha_helices'], sep="\s+", skiprows=skiprows("ss"))
        # print(chaindata)

        col_means=np.round(np.mean(chaindata, axis=0))
        col_means[0] = 0
        print(legends[pdbnum])
        print(col_means)

        # mean_centered = chaindata - col_means
        # mean_centered = mean_centered.astype({"time":'int', 'loops':'int', 'breaks':'int', 'bends':'int', 'turns':'int', 'pp_helices':'int', 'pi_helices':'int', '3_10_helices':'int', 'beta_strands':'int', 'beta_bridges':'int', 'alpha_helices':'int'})
        # # print(mean_centered)
        # mean_centered.to_csv(pdb+"/results/mean_centered_"+chain+".txt", sep="\t")

        time = chaindata['time']
        Loops = chaindata['loops']                     
        Bends = chaindata['bends']                        
        Turns = chaindata['turns']                        
        Helices_3_10  = chaindata['3_10_helices']                
        Alpha_helices = chaindata['alpha_helices'] 

        listofss=[Loops, Bends, Turns, Helices_3_10, Alpha_helices]
        legends=['Loops', 'Bends', 'Turns', '3-10 H', 'Alpha H']
        
        plt.suptitle('Number of secondary structural elements in chain '+chain+ ' with time', fontsize=18)
        for ss in range(len(listofss)):
            plt.subplot(len(listofss), 1, ss+1)
            plt.plot(time, listofss[ss])
            plt.ylabel(legends[ss], fontsize=18)
            plt.grid('both')

        plt.xlabel('Time (ps)',fontsize=18)
       
        # plt.ylabel('Number of secondary structural elements')
    plt.savefig(resultsdir+"ssplot_"+chain+".png")
    plt.close() 

def traj_correlations(pdblist):
    """get pearson's correlation coefficient from trajectory coordinate information (need to copy traj data)"""

def plot_cluster_rmsd_distributions(pdblist, simtime):
    """
    Plots the distribution of RMS with respect to frames
    """
    # Open file
    suffixes = ['A', 'B', 'all']
    for suffix in suffixes:
        for pdbnum in range(len(pdblist)):
            filepath=pdblist[pdbnum]+"/"+str(simtime)+"_ns/results_"+str(simtime)+"/clusters/rmsd-dist_"+suffix[-1]+".xvg"

def plot_terminal_RMSDs(prop, simtime):
    """
    Plots terminal RMSDs and their gradients with respect to time
    """
    cols = get_cols(prop)
    suffixes = ['_ChainA_N', '_ChainB_C', '_ChainB_N', '_Rf_terminal']
    pdblist, legends, sigma = get_pdb_list(simtime)

    for suffix in suffixes:
        plt.figure(figsize=(10,8))
        for pdbnum in range(len(pdblist)):

            filename = pdblist[pdbnum]+"/"+str(simtime)+"_ns/results_"+str(simtime)+"/"+prop + suffix +".xvg"
            resultsdir = "plots/terminal_RMSDs/"
            data = pd.read_csv(filename, sep="\s+", names=cols, skiprows=skiprows(prop))

            xaxis = data['xaxis']
            yaxis = data['property']  

            # Find the gradient of RMSDs with respect to time
            y_smooth = gaussian_filter1d(yaxis, sigma)
            gradient = np.gradient(y_smooth, xaxis)

            # plt.plot(xaxis, gradient, linewidth=2, label=legends[pdbnum], alpha=0.7)
            # plt.xlabel('Time (ps)',fontsize=20)
            # # print("Mean "+prop+" ("+legends[pdbnum]+"): "+str(np.mean(yaxis))+";\t Std_dev: "+str(np.std(yaxis)))
            # plt.ylabel(plot_yaxes[prop], fontsize=20)
            # plt.ylim([-0.003, 0.003])
            # plt.grid('both')
            # plt.legend()

            # Plot the original and smoothed data
            plt.subplot(2, 1, 1)
            plt.plot(xaxis, y_smooth, linewidth=2, label=legends[pdbnum])
            plt.legend(fontsize=12)
            plt.ylabel(plot_yaxes[prop], fontsize=20)
            plt.grid('both')

            # Plot the derivative of the smoothed data
            plt.subplot(2, 1, 2)
            plt.plot(xaxis, gradient, label=legends[pdbnum])
            plt.legend(fontsize=12)
            plt.grid('both')
            plt.ylabel('Gradient_RMSD')
            
            plt.xlabel('Time (ps)',fontsize=20)
            plt.tight_layout()
            

        plotfile=prop+suffix+"_grad_"+str(simtime)+"ns.png"

        plt.savefig(resultsdir+plotfile)
        plt.close()

def plot_all_transitions(prop, simtime):
    """
    For each system, plots all transitions (N termini of both chains, C terminus of B chain, and Rf terminus)
    to check whether they are in tandem, to understand the rate of conformational transitions 
    """

    cols = get_cols(prop)
    suffixes = ['_ChainA_N', '_ChainB_C', '_ChainB_N', '_Rf_terminal']
    pdblist, legends, sigma = get_pdb_list(simtime)

    plt.figure(figsize=(14,8))
    for pdbnum in range(len(pdblist)):
        for suffix in suffixes:

            filename = pdblist[pdbnum]+"/"+str(simtime)+"_ns/results_"+str(simtime)+"/"+prop + suffix +".xvg"
            resultsdir = "plots/terminal_RMSDs/"
            data = pd.read_csv(filename, sep="\s+", names=cols, skiprows=skiprows(prop))

            xaxis = data['xaxis']
            yaxis = data['property']  

            # Find the gradient of RMSDs with respect to time
            y_smooth = gaussian_filter1d(yaxis, sigma)
            gradient = np.gradient(y_smooth, xaxis)

            # Plot the original and smoothed data
            plt.subplot(1, 2, pdbnum+1)
            plt.plot(xaxis, y_smooth, linewidth=2, label=suffix[1:], alpha=0.7)
            plt.ylabel(plot_yaxes[prop], fontsize=15)
            plt.grid('both')
            plt.title(legends[pdbnum], fontsize=15)
            plt.ylim([0.1,0.5])
            # # Plot the derivative of the smoothed data
            # plt.subplot(2, 1, 2)
            # plt.plot(xaxis, gradient, label=legends[pdbnum])
            plt.legend(fontsize=15)
            # plt.grid('both')
            # plt.ylabel('Gradient_RMSD')
            
            plt.xlabel('Time (ps)',fontsize=15)
            plt.tight_layout()
        

    plotfile=prop+str(simtime)+"ns.png"
    # plt.legend(ncol=2, loc='upper center', bbox_to_anchor=(0.5, 0.95), fontsize=14)
    plt.savefig(resultsdir+plotfile)
    plt.close()

def plot_all_distances(prop, simtime):
    """
    For each system, plots all important distances (N-C termini, B5-B13, B26-B12, B28-B8)
    to check whether they are in tandem, to understand the rate of conformational transitions 
    """

    cols = get_cols(prop)
    prefixes = ['NC_', 'B5B13_', 'B26B12_', 'B28B8_']
    pdblist, legends, sigma = get_pdb_list(simtime)

    plt.figure(figsize=(14,8))
    for pdbnum in range(len(pdblist)):
        for prefix in prefixes:

            filename = pdblist[pdbnum]+"/"+str(simtime)+"_ns/results_"+str(simtime)+"/"+prefix+prop+".xvg"
            resultsdir = "plots/distances/"
            data = pd.read_csv(filename, sep="\s+", names=cols, skiprows=skiprows(prop))

            xaxis = data['xaxis']
            yaxis = data['property']  

            # Find the gradient of RMSDs with respect to time
            y_smooth = gaussian_filter1d(yaxis, sigma)
            gradient = np.gradient(y_smooth, xaxis)

            # Plot the original and smoothed data
            plt.subplot(1, 2, pdbnum+1)
            plt.plot(xaxis, y_smooth, linewidth=2, label=prefix[:-1], alpha=0.7)
            plt.ylabel(plot_yaxes[prop], fontsize=15)
            plt.grid('both')
            plt.title(legends[pdbnum], fontsize=15)
            plt.ylim([0,2.5])
            # # Plot the derivative of the smoothed data
            # plt.subplot(2, 1, 2)
            # plt.plot(xaxis, gradient, label=legends[pdbnum])
            plt.legend(fontsize=15)
            # plt.grid('both')
            # plt.ylabel('Gradient_RMSD')
            
            plt.xlabel('Time (ps)',fontsize=15)
            plt.tight_layout()
        

    plotfile=prop+str(simtime)+"ns.png"
    # plt.legend(ncol=2, loc='upper center', bbox_to_anchor=(0.5, 0.95), fontsize=14)
    plt.savefig(resultsdir+plotfile)
    plt.close()

def plot_heatmaps():
    """
    Plots heatmaps of RMSDs of structural and thermodynamic properties for different simulations 
    """
    import seaborn as sns
    # Sample data
    data = {
        'A': np.random.rand(10),
        'B': np.random.rand(10),
        'C': np.random.rand(10)
    }
    df = pd.DataFrame(data)

    # Plot heat maps for each column
    for column in df.columns:
        plt.figure(figsize=(5, 1))
        sns.heatmap(df[[column]].T, annot=True, cmap='coolwarm', cbar=False)
        plt.title(f'Heatmap for {column}')
        plt.show()

    # Main

def mainplotcode(simtime):
    """
    Calls all plotcodes and scripts to get statistics of all properties of the simulated systems
    """
    # Define all properties of interest
    energyprops = ['potential', 'kinetic', 'total_energy']
    structprops = ['rmsd', 'rg', 'rmsf']
    thermoprops = ['density', 'pressure', 'temperature']

    # for p in thermoprops:
    #     plotfile=p+".png"
    #     plot_thermo(p, plotfile, simtime)

    for p in structprops:
        plot_structs(p, simtime)

    # plotfile="energies.png"
    # plot_energies(plotfile, simtime)

    # normalise_and_compare('rg', simtime)

    # allprops = energyprops +  thermoprops
    # get_statistics(allprops, simtime)

    # plot_ss("A", simtime)
    # plot_ss("B", simtime)

    # plot_terminal_RMSDs("rmsd", simtime)
    # plot_all_transitions("rmsd", simtime)
    # plot_all_distances('distance', simtime)


mainplotcode(simtime)





