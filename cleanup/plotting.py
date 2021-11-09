import matplotlib.pyplot as plt
import process_dump_file as process
import numpy as np

##########
# Plotting gyration radius against times
#########

def gyration_time_plot(infile):
    """
    Takes in gyration-time info and plots as a line graph
    @param {string} file name
    """
    # unpack file
    dumpfilename = infile
    inf = open(dumpfilename, 'r')
    # read first line of the file to skip titles
    line = inf.readline()
    # create plotting list
    timestep = []
    Rg_list = []
    # find number of lines in our gyration-time file
    line_number = process.lines_in_file(infile)
    # loop over the rest of the lines in the g-t file to populate our lists
    for i in range(line_number-1):
        line = inf.readline()
        line = line.split()
        timestep.append(int(line[0]))
        Rg_list.append(np.float64(line[1]))

    #plot graph
    a = max(timestep)

    plt.plot(timestep,Rg_list,"-r")
    plt.xlabel("Time")
    plt.xticks(np.arange(min(timestep), max(timestep)+1, 50000))
    plt.ylabel("Radius of Gyration")
    plt.show()
    inf.close()

def polymer_plot(infile):
    """
    takes in the polymer-time info and plots
    """
    # unpack file
    dumpfilename = infile
    inf = open(dumpfilename, 'r')
    # read first line of the file to skip titles
    line = inf.readline()
    # create plotting list
    # find number of lines in our gyration-time file
    timestep = []
    polymer_list = []
    # find the number of lines in the file
    line_number = process.lines_in_file(infile)
    # loop over the rest of the lines in the g-t file to populate our lists
    for i in range(line_number-1):
        line = inf.readline()
        line = line.split()
        timestep.append(int(line[0]))
        polymer_list.append(int(line[1]))
    #plot graph

    plt.plot(timestep,polymer_list,"-k")
    plt.xlabel("Time")
    plt.xticks(np.arange(min(timestep), max(timestep)+1, 50000))
    plt.ylabel("Number of proteins bound to the polymer")
    plt.show()
    inf.close()

def cluster_plot(infile):
    """
    takes in the protein-time info and plots
    """
    # unpack file
    dumpfilename = infile
    inf = open(dumpfilename, 'r')
    # read first line of the file to skip titles
    line = inf.readline()
    # create plotting list
    timestep = []
    protein_list = []
    # find number of lines in our gyration-time file
    line_number = process.lines_in_file(infile)
    # loop over the rest of the lines in the g-t file to populate our lists
    for i in range(line_number-1):
        line = inf.readline()
        line = line.split()
        timestep.append(int(line[0]))
        protein_list.append(int(line[1]))
    #plot graph
    plt.plot(timestep,protein_list,"-r")
    plt.xlabel("Time")
    plt.xticks(np.arange(min(timestep), max(timestep)+1, 50000))
    plt.ylabel("Number of proteins in a cluster")
    plt.show()
    inf.close()


def combo_plot(Rgfile, polymerfile, proteinfile):
    """
    """
    inf = open(Rgfile,'r')
    # read first line of the file to skip titles
    line = inf.readline()
    # create plotting list
    timestep = []
    Rg_list = []
    # find number of lines in our gyration-time file
    line_number = process.lines_in_file(Rgfile)

    # loop over the rest of the lines in the g-t file to populate our lists
    for i in range(line_number-1):
        line = inf.readline()
        line = line.split()
        timestep.append(int(line[0]))
        Rg_list.append(np.float64(line[1]))

    inf.close()

    inf = open(polymerfile,'r')
    # read first line of the file to skip titles
    line = inf.readline()
    # create plotting list
    timestep = []
    polymer_list = []
    # find number of lines in our gyration-time file
    line_number = process.lines_in_file(polymerfile)

    # loop over the rest of the lines in the g-t file to populate our lists
    for i in range(line_number-1):
        line = inf.readline()
        line = line.split()
        timestep.append(int(line[0]))
        polymer_list.append(int(line[1]))

    inf.close()

    # unpack file

    inf = open(proteinfile,'r')
    # read first line of the file to skip titles
    line = inf.readline()
    # create plotting list
    timestep = []
    protein_list = []
    # find number of lines in our gyration-time file
    line_number = process.lines_in_file(proteinfile)
    # loop over the rest of the lines in the g-t file to populate our lists
    for i in range(line_number-1):
        line = inf.readline()
        line = line.split()
        timestep.append(int(line[0]))
        protein_list.append(int(line[1]))
    inf.close()


    #print(len(polymer_list))
    #print(len(Rg_list))
    #print(len(timestep))


    #plot graph
    fig, (ax1,ax2,ax3) = plt.subplots(1,3, figsize = (10,10))
    ax1.plot(timestep,Rg_list,"-r")
    ax2.plot(timestep,polymer_list,"-m")
    ax3.plot(timestep,protein_list,"-k")
    ax1.set_xlabel("Time")
    ax2.set_xlabel("Time")
    ax3.set_xlabel("Time")
    ax1.set_xticks(np.arange(min(timestep), max(timestep)+1, 50000))
    ax2.set_xticks(np.arange(min(timestep), max(timestep)+1, 50000))
    ax3.set_xticks(np.arange(min(timestep), max(timestep)+1, 50000))
    ax1.set_ylabel("Radius of gyration")
    ax2.set_ylabel("Number of proteins in attached to a polymer")
    ax3.set_ylabel("Number of proteins in a cluster")
    plt.show()
    plt.tight_layout()
    return protein_list, polymer_list, Rg_list, timestep


combo_plot("gyration_time","polymer_bound","protein_cluster")
#gyration_time_plot("gyration_time")
#polymer_plot("polymer_bound")
#cluster_plot("protein_cluster")
