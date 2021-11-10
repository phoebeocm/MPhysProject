from process_dump_file import Atom
from process_dump_file import *
import numpy as np
from averages import *

############################################################################
#### Start of the main program

def main():

    #dumpfilename = 'dump.DNA+proteins'   # this is hardcoded here, could instead read as command line argument
    pol_attraction = input("What's your polymer-protein attraction? ")
    prot_attraction = input("What's your protein-protein attraction? ")

    dumpfilename = 'dump.pp_%s_pc_%s' % (prot_attraction, pol_attraction)

    Natoms = 700   # this is hardcoded here, could instead read as command line argument

    Nlines = lines_in_file(dumpfilename)  # get length of file
    Nframes = int(Nlines / (Natoms+9))         # there are 9 header lines in each frame

    # open the input file
    inf = open(dumpfilename, 'r')

    # open the output files and print a header
    #ouf_rg = open('r_g_1.dat', 'w')
    #ouf_rg.write("# frame number, radius of gyration\n")
    # open file to write out gyration and time data
    #gt = open("gyration_time",'w')
    #gt.write("time, radius of gyration\n")
    # open file to write bound polymer data
    polymer = open("polymer_bound_pp%s_pc%s" %(prot_attraction,pol_attraction),"w")
    polymer.write("time, number of proteins bound to the polymer\n")
    # open file to write protein cluster data
    cluster = open("protein_cluster_pp%s_pc%s" %(prot_attraction, pol_attraction),"w")
    cluster.write("time, number of proteins in clusters\n")
    # open file to write the pp rdf
    rdf_protein = open("rdf_protein","w")
    rdf_protein.write("bin, g(r)\n")
    # open file to write the pc rdf
    rdf_polymer = open("rdf_polymer", "w")
    rdf_polymer.write("bin, g(r)\n")
    # open file to write cc rdf
    rdf_cc = open("rdf_cc", "w")
    rdf_cc.write("bin, g(r)\n")

    #info for plotting rdf
    bins = 30
    g_of_r_pp = np.zeros(bins)
    g_of_r_pc = np.zeros(bins)
    g_of_r_cc = np.zeros(bins)
    volumes = np.zeros(bins)

    # go through the file frame by frame
    for frame in range(Nframes):

        #read the frame, unwrapping periodic coordinates
        atoms, L, timestep = readframe_unwrap(inf,Natoms)

       # unwarp period boundary coordinates -- needed for radius of gyration
        for i in range(len(atoms)):
            atoms[i].unwrap(L)

        # calculate radius of gyration
        #Rg = radius_of_gyration(atoms,L)

        # calculating the time?
        time = int(timestep[0])

        # calculate the number of proteins bound to the polymer in this frame
        polymer_number = bound_polymer(atoms)

        # calculate the number of proteins in clusters in this frame
        protein_number = protein_clusters(atoms)

        # calculate the radial distribution function

        proteins, polymers, g_of_r_pp, g_of_r_pc, g_of_r_cc, volumes, box_volume, radii = rdf(atoms, bins, g_of_r_pp, g_of_r_pc, g_of_r_cc,  volumes)

        # output the frame index, radius of gyration to a file
        #ouf_rg.write("%i %.5f\n"%(frame+1,Rg) )
        # output timestep, radius of gyration to a file
        #gt.write("%i %.5f\n"%(time,Rg))
        # output timestep, number of proteins bound to a polymer to file
        polymer.write("%i %i\n"%(time,polymer_number))
        # output timestep, number of proteins in a cluster to file
        cluster.write("%i %i\n"%(time,protein_number))


    #ouf_rg.close()
    #gt.close()
    polymer.close()
    cluster.close()

    ############################################################################
    # THIS SECTION IS MEANT TO SUBTRACT CONTRIBUTIONS PRIOR TO EQUILIBRIUM POINT
    equilibrium = find_equilibrium(pol_attraction,prot_attraction)
    g_of_r_pp_2 = np.zeros(bins)
    g_of_r_pc_2 = np.zeros(bins)
    g_of_r_cc_2 = np.zeros(bins)
    volumes_2 = np.zeros(bins)

    for i in range(equilibrium):
        atoms, L, timestep = readframe_unwrap(inf,Natoms)

       # unwarp period boundary coordinates -- needed for radius of gyration
        for i in range(len(atoms)):
            atoms[i].unwrap(L)

        proteins, polymers, g_of_r_pp_2, g_of_r_pc_2, g_of_r_cc_2, volumes_2, box_volume, radii = rdf(atoms, bins, g_of_r_pp_2, g_of_r_pc_2, g_of_r_cc_2,  volumes_2)

    g_of_r_cc = g_of_r_cc - g_of_r_cc_2
    g_of_r_pp = g_of_r_pp - g_of_r_pp_2
    g_of_r_pc = g_of_r_pc - g_of_r_pc_2


    inf.close()

     ###########################################################################
    for i, value in enumerate(g_of_r_pp):

        #g_of_r[i] = (value/volumes[i])
        g_of_r_pp[i] = (value/volumes[i])*(box_volume/len(proteins))
        rdf_protein.write("%i %.5f\n"%(radii[i],g_of_r_pp[i]))

    for i, value in enumerate(g_of_r_pc):
        g_of_r_pc[i] = (value/volumes[i])*(box_volume/len(polymers)) #should I divide by prots or pols?
        rdf_polymer.write("%i %.5f\n"%(radii[i],g_of_r_pc[i]))

    for i, value in enumerate(g_of_r_cc):
        g_of_r_cc[i] = (value/volumes[i])*(box_volume/len(polymers)) # divide by what?
        rdf_cc.write("%i %.5f\n"%(radii[i],g_of_r_cc[i]))

    # close the files

    rdf_protein.close()
    rdf_polymer.close()
    rdf_cc.close()

    # Can now reach into the bound values we created above to calculate the average
    #polymer = open("polymer_bound_pp%s_pc%s" %(prot_attraction,pol_attraction),"r")

    #easy_mean(pol_attraction,prot_attraction)

    #polymer.close()
    #cluster.close()

# Finished!
if __name__ == "__main__":
    main()
