#!/usr/bin/env python
#

#For all exercises, instead of reading a local file, use PDB id as input and download the structure from PDB. 

from Bio.PDB import PDBList

def get_structure_file(pdb_id):
    pdbl = PDBList()
    filename = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='pdb')
    return filename

#Determine the list of pairs of residues whose CA atoms are closer than a given distance
#Parameters: PDB file name, distance.

""" Simple program to search contacts """

from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

def pairs_residues(PDB_file, distance):

    MAXDIST = distance

    parser = PDBParser(PERMISSIVE=1)

    # load structure from PDB file

    filename = get_structure_file(PDB_file)
    st = parser.get_structure('struc', filename)
    select = []

    #Select only CA atoms

    for at in st.get_atoms():
        if at.id == 'CA':
            select.append(at)
            print(f"ATOM: {at.get_parent().get_resname()}, {at.get_parent().id[1]}, {at.id}")

    # Preparing search
    nbsearch = NeighborSearch(select)

    print("NBSEARCH:")

    #Searching for contacts under HBLNK

    ncontact = 1

    for at1, at2 in nbsearch.search_all(MAXDIST):
        print(f"Contact: {ncontact}")
        print(f"at1: {at1}, {at1.get_serial_number()}, {at1.get_parent().get_resname()}")
        print(f"at2: {at2}, {at2.get_serial_number()}, {at2.get_parent().get_resname()}")
        print()
        ncontact += 1



#!/usr/bin/env python
#

#Generate a list of all atoms for a given residue number 
#Parameters: PDB file name, Residue number (Including Chain if applicable)

""" Simple program to print ARG residues iteration over atoms """

from Bio.PDB.PDBParser import PDBParser

def atoms_residue_number(PDB_file, residue_number):
    parser = PDBParser()

    filename = get_structure_file(PDB_file)
    st = parser.get_structure('struc', filename)

    selected = []

    for at in st.get_atoms():
        if at.get_parent().id[1]==residue_number:
            selected.append(at)

    print("Coordinates:")
    for atom in selected:
        print(f"{atom.get_parent().get_resname()}, {atom.get_parent().id}, {atom.get_name()}, {atom.get_coord()}")



#Determine all possible hydrogen bonds (Polar atoms at less than 3.5 Å). 
#Parameters: PDB file name. Optional: cut-off distance (defaults to 3.5)

""" Simple program to search contacts """

from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

def hydrogen_bonds(PDB_file, distance=3.5):

    MAXDIST = distance

    parser = PDBParser(PERMISSIVE=1)

    # load structure from PDB file

    filename = get_structure_file(PDB_file)
    st = parser.get_structure('struc', filename)

    select = []

    for at in st.get_atoms():
        if at.id in ('O','N','S'):
            select.append(at)
            print(f"ATOM: {at.get_parent().get_resname()}, {at.get_parent().id[1]}, {at.id}")

    # Preparing search
    nbsearch = NeighborSearch(select)

    print("NBSEARCH:")

    #Searching for contacts under HBLNK

    ncontact = 1

    for at1, at2 in nbsearch.search_all(MAXDIST):
        print(f"Contact: {ncontact}")
        print(f"at1: {at1}, {at1.get_serial_number()}, {at1.get_parent().get_resname()}")
        print(f"at2: {at2}, {at2.get_serial_number()}, {at2.get_parent().get_resname()}")
        print()
        ncontact += 1


#Generate a list of all CA atoms of given residue type with coordinates
#Parameters: PDB file name, residue type. 
#Optional: accept residue codes in one- or three-letter formats automatically

""" Simple program to search contacts """

from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

def residue_type_coord(PDB_file, residue_type):

    parser = PDBParser(PERMISSIVE=1)

    # load structure from PDB file

    filename = get_structure_file(PDB_file)
    st = parser.get_structure('struc', filename)

    select = []

    #Select only CA atoms

    for at in st.get_atoms():
        if at.id == 'CA' and at.get_parent().get_resname()==residue_type:
            select.append(at)
            print(f"ATOM: {at.get_parent().get_resname()}, {at.get_parent().id[1]}, {at.id}, {at.get_coord()}")


#Generate a list of backbone connectivity (i.e. which residues are linked by ordinary peptide bonds).
#Parameters: PDB file name. Optional: Cut-off distance for peptide bonds (defaults to 2.5)

""" Simple program to search contacts """

from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

def backbone_connectivity(PDB_file, distance=2.5):

    MAXDIST = distance

    parser = PDBParser(PERMISSIVE=1)

    # load structure from PDB file

    filename = get_structure_file(PDB_file)
    st = parser.get_structure('struc', filename)

    select = []

    #Select only CA atoms

    for at in st.get_atoms():
        if at.id in ('C','N'):
            select.append(at)
            print(f"ATOM: {at.get_parent().get_resname()}, {at.get_parent().id[1]}, {at.id}")

    # Preparing search
    nbsearch = NeighborSearch(select)

    print("NBSEARCH:")

    #Searching for contacts under HBLNK

    ncontact = 1

    for at1, at2 in nbsearch.search_all(MAXDIST):
        if at1.get_parent() != at2.get_parent():
            print(f"Contact: {ncontact}")
            print(f"at1: {at1}, {at1.get_serial_number()}, {at1.get_parent().get_resname()}")
            print(f"at2: {at2}, {at2.get_serial_number()}, {at2.get_parent().get_resname()}")
            print()
            ncontact += 1

#Id 4, but for disulphide bonds.

""" Simple program to search contacts """

from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser

def disulphide_bonds(PDB_file, distance=1.9):

    MAXDIST = distance

    parser = PDBParser(PERMISSIVE=1)

    # load structure from PDB file

    filename = get_structure_file(PDB_file)
    st = parser.get_structure('struc', filename)

    select = []

    for at in st.get_atoms():
        if at.id == 'SG' and at.get_parent().get_resname()=='CYS':
            select.append(at)
            print(f"ATOM: {at.get_parent().get_resname()}, {at.get_parent().id[1]}, {at.id}")

    if len(select) < 2:
        print("No possible disulphide bonds found -> Error")
        return

    # Preparing search
    nbsearch = NeighborSearch(select)

    print("NBSEARCH:")

    #Searching for contacts under HBLNK

    ncontact = 1

    for at1, at2 in nbsearch.search_all(MAXDIST):
        print(f"Contact: {ncontact}")
        print(f"at1: {at1}, {at1.get_serial_number()}, {at1.get_parent().get_resname()}")
        print(f"at2: {at2}, {at2.get_serial_number()}, {at2.get_parent().get_resname()}")
        print()
        ncontact += 1


#Print distances between all atom pairs of two given residues
#Parameters: PDB file name, Residue 1, Residue 2 

""" Simple program to print distances between atoms """

from Bio.PDB.PDBParser import PDBParser
import numpy as np

def dist_by_residue(PDB_file, residue1, residue2):
    parser = PDBParser()

    filename = get_structure_file(PDB_file)
    st = parser.get_structure('struc', filename)

    #Selection of Residue 10 of Chain A
    res10 = st[0]["A"][residue1]
    #Selection of Residue 20 of Chain A
    res20 = st[0]["A"][residue2]

    print(f"Residue", {residue1}, "is", res10.get_resname())
    print(f"Residue", {residue2}, "is", res20.get_resname())

    print("\nAtom1 Atom2 dist1 dist2\n-------------------------")
    for at10 in res10.get_atoms():      # Replace get_atoms with get_atom if you get an Error!
        for at20 in res20.get_atoms():
            dist = at20 - at10     # Direct procedure with (-) to compute distances
            vector = at20.coord - at10.coord  # Or using numpy coordinates
            distance = np.sqrt(np.sum(vector ** 2))
            print(at10, at20, dist, distance)



##TEST 

#1. 
pairs_residues('1UBQ', 4)

#2.
atoms_residue_number('1UBQ', 8)

#3.
hydrogen_bonds('1UBQ')

#4.
residue_type_coord('1UBQ', 'ARG')

#5.
backbone_connectivity('1UBQ')

#6.
disulphide_bonds('1UBQ')

#7.
dist_by_residue('1UBQ', 10, 20)