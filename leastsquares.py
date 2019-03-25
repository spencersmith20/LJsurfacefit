from numpy import linspace, power, random, inf, asarray
import pandas as pd
from math import sqrt
from ase import Atoms
from os import chdir

#calculated / look-up values
sotolon_vac = -86.96520401259681; graphene_vac = -305.49689755541885
carbon_sig = 3.55; carbon_eps = 0.066 / 627.509 #kcal/mol -> Hartree

#get the distance between two atom positions
def distance_between_coords(p1, p2):
    return sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)

#evaluate lennard jones potential
def lennard_jones(r, eps, sig):
    return 4 * eps * ((sig/r)**12 - (sig/r)**6)

#calculate all atom-atom energies outside of the cutoff. self-interaction in molecule is excluded
def get_energy(epsilon, sigma, atoms, m):
    E = [];
    for i in range(len(atoms)):
        p1 = atoms[i]; e = []

        #calculate p1's interactions with p2 thru pn, then increment p1 to p2, dealing with p3 thru pn
        for j in range(i+1,len(atoms)):
            p2 = atoms[j]; r = distance_between_coords(p1, p2)
            if r > 10: #cut-off
                continue;

            if i < m: # deal with sotolon-carbon interaction ... will need to deal with an array of parameters
                e.append(lennard_jones(r, epsilon, sigma))
            else: # deal with carbon-carbon interaction
                e.append(lennard_jones(r, carbon_eps, carbon_sig))

            # residual to be minimized ... v[i] is the electronic energy from this run
            E.append(sum(e))

    #return summed energy of configuration
    return(sum(E))

#load files
xyz_file = open('filenames.txt'); xyz = [line.strip('\n') for line in xyz_file.readlines()]
xyz_file.close()

#initial guesses
sot_eps = 1e-7; sot_sig = 1

# mixing rules for LJ parameters
mix_eps = sot_eps * carbon_eps; mix_sig = (sot_sig + carbon_sig) / 2

#read energies
energy_file = open('energies.txt')
v = [float(line.strip('\n').split('\t')[2]) for line in energy_file.readlines()]
energy_file.close()

#for each configuration of the molecule, calculate total E and find residual
energies = []; m = 1
for dir in xyz:

    #move into particular folder
    chdir(dir)

    #get appropriate coordinates
    df = pd.read_csv(dir + '.xyz', header=2, names=['type', 'x', 'y', 'z'], sep='\s+')
    surf_z = df.iloc[0,3]; n = len(df.iloc[:,0])

    #separates molecule and surface by checking z coordinate
    for i in range(n):
        if df.iloc[i, 3] != surf_z:
            surf_df = df.iloc[:i,:]; mol_df = df.iloc[i:,:]; break

    #create Atoms objects for surface and molecule
    surf = Atoms(symbols=surf_df.iloc[:,0].values, positions=surf_df.iloc[:,1:].values)
    mol = Atoms(symbols=mol_df.iloc[:,0].values, positions=mol_df.iloc[:,1:].values)

    #list for iteration of energies. m is number of bodies interaction (currently dealing with sphere approximation)
    atoms = mol.get_center_of_mass() + surf.positions

    #get the sum squared error for this set of parameters
    energies.append(get_energy(mix_eps, mix_sig, atoms, m))

    #go into bulk directory
    chdir('..')

#get residuals for each configuration
residuals = []
for i in range(len(v)):
    residuals.append((energies[i] - v[i])**2)
