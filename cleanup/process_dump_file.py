#!/usr/bin/env python3

# Program to read in a dump file frame by frame
# and do some calculations.

# Example here is radius of gyration
# Here Rg of all atoms in the file are calculated

# Assumes the dump file has the following format
# id type x y z ix iy iz
# were x  etc. are atom coordinates
#      ix etc. are image flags for taking periodic boundaries into account
# need to adjust these to range from -L/2 to L/2


import numpy as np
import operator
import matplotlib.pyplot as plt
#from averages import *


class Atom:
    """ A Class for storing atom information """

    def __init__(self):
        """ Initialise the class """
        self.id = 0                              # id of the atom
        self.type = 0                            # type of the atom
        self.x = np.array([0.0,0.0,0.0],dtype=np.float64)     # position of the atom
        self.image = np.array([0,0,0],dtype=np.int32)         # image flags for atoms
        self.unwrap_flag = False

    def sep(self,B):
        """ Find separation of this Atom and Atom B """
        return np.sqrt( (self.x[0]-B.x[0])**2 +
                        (self.x[1]-B.x[1])**2 + (self.x[2]-B.x[2])**2 )

    def minus(self,B):
        """ Subtract B.x vector from self.x vector """
        AminusB = np.array([0.0,0.0,0.0],dtype=np.float64)
        for i in range(3):
            AminusB[i] = self.x[i] - B.x[i]
        return AminusB

    def xdot(self,B):
        """ Find dot product of position x of this Atom and Atom B """
        AdotB = np.array([0.0,0.0,0.0],dtype=np.float64)
        # TO DO : find AdotB
        return AdotB

    def unwrap(self,L):
        """ Unwraps the coordinates for periodic box, and overwrites x """
        if not self.unwrap_flag:   # first check it has not already been done
            for j in range(3):
                self.x[j] = self.x[j] + self.image[j]*L[j] # unwrap
            unwrap_flag = True



def readframe_unwrap(infile,N):
    """ Read a single frame of atoms from a dump file
        Rescale coordinates to be in rnage -L/2 to L/2
        DOES NOT Unwrap corrdinates for periodic box """

    atoms = []
    L = []
    timestep = []

    # read in the 9 header lines of the dump file

    for i in range(9):
        line = infile.readline()
        line = line.split()
        if i == 1:
            # get the timestep
            timestep.append(line[0])
        if i==5 or i==6 or i==7:
            # get the box size
            L.append( np.float64(line[1]) - np.float64(line[0]) )

    # now read the atoms
    for i in range(N): #N
        line = infile.readline()
        #need to bake in some way of avoiding strings
        line = line.split()
        newatom = Atom()#lower down, this is being looped over the whole list file, including
        #strings

        newatom.id = int(line[0])
        newatom.type = int(line[1])
        for j in range(3):
            newatom.x[j] = np.float64(line[j+2]) # scale
        atoms.append(newatom)

#    make sure atoms are sorted by id
    atoms.sort(key=operator.attrgetter('id'))

    return atoms,L,timestep


def lines_in_file(filename):
    """ Get the number of lines in the file """

    with open(filename) as f:
        for i, l in enumerate(f):
            pass

    return i + 1

def protein_list(atoms):
    proteins = []
    for i in range(len(atoms)): #
        if atoms[i].type == int(2):
            proteins.append(atoms[i])
    return proteins

def polymer_list(atoms):
    polymer = []
    for i in range(len(atoms)): #
        if atoms[i].type == int(1):
            polymer.append(atoms[i])
    return polymer


def radius_of_gyration(atoms,L):
    #Calculate the radius of gytation -- Rg^2 = (1/N) sum ( r_k - r_mean)^2
    #remember to unwrap periodic boundaries "

    # get mean position
    r_mean = np.zeros(3,dtype=np.float64)
    for i in range(len(atoms)):
        r_mean[0] += atoms[i].x[0]
        r_mean[1] += atoms[i].x[1]
        r_mean[2] += atoms[i].x[2]
    r_mean = r_mean/len(atoms)

    # get Rg2
    Rg2 = 0.0
    for i in range(len(atoms)):
        Rg2 += np.sum( np.square( atoms[i].x - r_mean ) )
    Rg2 = Rg2/len(atoms)

    return np.sqrt( Rg2 )


def bound_polymer(atoms):
    """
    returns number of proteins bound to polymer
    """
    count = 0
    proteins = protein_list(atoms)
    polymer = polymer_list(atoms)
    for i in proteins:
        for j in polymer:
            if Atom.sep(i,j) < 1.8:
                count += 1
                break
    return count


def protein_clusters(atoms):
    """
    returns the number of proteins in clusters
    """
    count = 0
    proteins = protein_list(atoms)
    polymer = polymer_list(atoms)
    for i in range(len(proteins)):
        for j in range(i+1,len(proteins)):
            if Atom.sep(proteins[i],proteins[j]) < 1.8:
                count += 1
                break
    return count

def volume(radius):
    return 4/3*np.pi*radius**3

def rdf(atoms, bins, g_of_r_pp, g_of_r_pc, g_of_r_cc, volumes):

    box_size = 10.0
    box_volume = box_size**3
    dr = 10.0/bins # or 1.8/bins
    radii = np.linspace(0.0, bins * dr, bins)

    ## find the protein and polymer lists
    proteins = protein_list(atoms)
    polymers = polymer_list(atoms)

    ## calculate protein-protein
    for i in range(len(proteins)):
        # find shell volumes
        for j in range(bins):
            r1 = j*dr
            r2 = r1 + dr
            v1 = volume(r1)
            v2 = volume(r2)
            shell_volume = v2 - v1
            volumes[j] += shell_volume

        for j in range(i+1,len(proteins)):
            sep = Atom.sep(proteins[i],proteins[j])
            index = int(sep / dr)
            if 0 < index < bins:
                g_of_r_pp[index] += 2

        for j in range(len(polymers)):
            sep = Atom.sep(proteins[i],polymers[j])
            index = int(sep / dr)
            if 0 < index < bins:
                g_of_r_pc[index] += 2
    # persistence length?
    for i, polymer in enumerate(polymers):
        for polymer2 in polymers[i:]:
            sep = Atom.sep(polymer,polymer2)
            if sep > 4.0:
                index = int(sep/dr)
                if 0 < index < bins:
                    g_of_r_cc[index] += 2

    return proteins, polymers, g_of_r_pp, g_of_r_pc, g_of_r_cc, volumes, box_volume, radii



# Finished!
