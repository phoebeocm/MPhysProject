#from plotting import combo_plot
import numpy as np
from time import sleep
from tqdm import tqdm
import process_dump_file as process


#protein_list, polymer_list, Rg_list, timestep = combo_plot("gyration_time","polymer_bound","protein_cluster")
"""

def mean(value_list, polymer_list, Rg_list, timestep_list):

    ## find a way to find the equilibrium mean_values
    particle_val = 0
    equilibrium_point = 0
    cutoff = 0
    # find the equilibrium point
    type = input("What average are you trying to measure? (ans: Rg, polymer clusters, particle clusters)")
    if type == "polymer clusters":
        cutoff = 7
        value_list = polymer_list
    elif type == "particle clusters":
        cutoff = 7
        value_list = protein_list
    elif type == "Rg":
        cutoff = 1.0
        value_list = Rg_list
    for i in range(len(value_list)):
        if max([abs(value_list[i]-value_list[i+j]) for j in range(20)])  < cutoff: # or some other tolerance
            equilibrium_point = i
            break
    # use equilibrium point to calculate the mean
    standard_dev = []
    for i in range(equilibrium_point, len(value_list)):
        particle_val += value_list[i]
        standard_dev.append(value_list[i])
    standard = np.std(standard_dev)
    mean = particle_val/(len(value_list)- (equilibrium_point+1))
    print("equilibrium is " + str(equilibrium))
    return mean, standard, equilibrium
"""
#mean(protein_list, polymer_list, Rg_list, timestep)


def easy_mean(pol_attraction, prot_attraction):
    clusterfilename = "protein_cluster_pp%s_pc%s" %(prot_attraction, pol_attraction)
    cluster = open(clusterfilename,"r")
    particle_val = 0
    equilibrium_point = 0
    cutoff = 1.8
    standard_dev = []
    # read first line of the file to skip titles
    line = cluster.readline()
    # create plotting list
    timestep = []
    protein_list = []
    # find number of lines in our gyration-time file
    line_number = process.lines_in_file(clusterfilename)
    # loop over the rest of the lines in the g-t file to populate our lists
    for i in range(line_number-1):
        line = cluster.readline()
        line = line.split()
        timestep.append(int(line[0]))
        protein_list.append(int(line[1]))
    for i in range(len(protein_list)):
        if max([abs(protein_list[i]-protein_list[i+j]) for j in range(5)])  < cutoff: # or some other tolerance
            equilibrium_point = i
            break
    for i in range(equilibrium_point, len(protein_list)):
        particle_val += protein_list[i]
        standard_dev.append(protein_list[i])
    standard = np.std(standard_dev)
    mean = particle_val/(len(protein_list) - (equilibrium_point+1))
    print(mean)
    print(standard)
    print(prot_attraction)
    total_file = open("averages_file_polymer_attraction_%s"%(pol_attraction),"a")
    total_file.write("%s %.5f %.5f\n"%(prot_attraction, mean, standard))
    total_file.close()

easy_mean(1.0, 1.0)
