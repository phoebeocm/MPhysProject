from plotting import combo_plot
import numpy as np
from time import sleep
from tqdm import tqdm


protein_list, polymer_list, Rg_list, timestep = combo_plot("gyration_time","polymer_bound","protein_cluster")


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
    print("Mean is " + str(mean) + ", equilibrium point is " + str(equilibrium_point))

mean(protein_list, polymer_list, Rg_list, timestep)
