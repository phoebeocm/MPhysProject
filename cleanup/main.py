from process_dump_file import Atom
from process_dump_file import *
import numpy as np


############################################################################
#### Start of the main program
#print('provide input file name')
#dumpfilename = input()
def main():
    dumpfilename = 'dump.DNA+proteins'   # this is hardcoded here, could instead read as command line argument
    Natoms = 700   # this is hardcoded here, could instead read as command line argument

    outfile_Rg = 'r_g_1.dat'

    Nlines = lines_in_file(dumpfilename)  # get length of file
    Nframes = int(Nlines / (Natoms+9))         # there are 9 header lines in each frame

    # open the input file
    inf = open(dumpfilename, 'r')

    # open the output files and print a header
    ouf_rg = open(outfile_Rg, 'w')
    ouf_rg.write("# frame number, radius of gyration\n")
    # open file to write out gyration and time data
    gt = open("gyration_time",'w')
    gt.write("time, radius of gyration\n")
    # open file to write bound polymer data
    polymer = open("polymer_bound","w")
    polymer.write("time, number of proteins bound to the polymer\n")
    # open file to write protein cluster data
    cluster = open("protein_cluster","w")
    cluster.write("time, number of proteins in clusters\n")

    rdf_protein = open("rdf_protein","w")
    rdf_protein.write("bin, g(r)")

    # go through the file frame by frame
    for frame in range(Nframes):
        #read the frame, unwrapping periodic coordinates
        atoms, L, timestep = readframe_unwrap(inf,Natoms)

       # unwarp period boundary coordinates -- needed for radius of gyration
        for i in range(len(atoms)):
            atoms[i].unwrap(L)

        # calculate radius of gyration
        Rg = radius_of_gyration(atoms,L)

        # calculating the time?
        time = int(timestep[0])

        # calculate the number of proteins bound to the polymer
        polymer_number = bound_polymer(atoms)

        # calculate the number of proteins in clusters
        protein_number = protein_clusters(atoms)

        # calculate the radial distribution function
        bins = 20
        g_of_r = np.zeros(bins)
        volumes = np.zeros(bins)
        box_size = 20.0
        box_volume = 20.0**3
        dr = box_size/bins
        radii = np.linspace(0.0, bins * dr, bins)
        proteins = rdf(atoms,10,g_of_r,volumes,radii)
        # output some results
        ouf_rg.write("%i %.5f\n"%(frame+1,Rg) )
        gt.write("%i %.5f\n"%(time,Rg))
        polymer.write("%i %i\n"%(time,polymer_number))
        cluster.write("%i %i\n"%(time,protein_number))

    for i, value in enumerate(g_of_r):
        print("value is " + str(value))
        print("volume is " + str(volumes[i]))
        #g_of_r[i] = (value/volumes[i])
        #g_of_r[i] = (value/volumes[i])*(box_volume/len(proteins))
        #rdf_protein.write("%i %.5f\n"%(radii[i],g_of_r[i]))
    # close the files

    inf.close()
    ouf_rg.close()
    gt.close()
    polymer.close()
    cluster.close()
    rdf_protein.close()

# Finished!
if __name__ == "__main__":
    main()
