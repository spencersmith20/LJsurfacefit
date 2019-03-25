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

#for each configuration of the molecule, calculate total E and find residual
def get_residuals(dir_list, V, epsilon, sigma, m):

    energies = []; residuals = []
    for dir in dir_list:

        #move into particular folder
        chdir(dir)

        #get appropriate coordinates
        df = pd.read_csv(dir + '.xyz', header=2, names=['type', 'x', 'y', 'z'], sep='\s+')
        coords = df[['x', 'y', 'z']]; types = df['type']
        surf_z = coords['z'].iloc[0]; n = len(coords['z'])

        #separates molecule and surface by checking z coordinate
        for i in range(n):
            if df['z'].iloc[i] != surf_z:
                surf_coords = coords.iloc[:i,:]; mol_coords = coords.iloc[i:,:]
                surf_types = types.iloc[:i]; mol_types = types.iloc[i:]; break

        #create Atoms objects for surface and molecule
        surf = Atoms(positions=surf_coords.values, symbols=surf_types.values)
        mol = Atoms(positions=mol_coords.values, symbols=mol_types.values)

        #list for iteration of energies. m is number of bodies interaction (currently dealing with sphere approximation)
        atoms = mol.get_center_of_mass() + surf.positions

        #get the sum squared error for this set of parameters
        energies.append(get_energy(epsilon, sigma, atoms, m))

        #go into bulk directory
        chdir('..')

    #get residuals for each configuration
    for i, energ in enumerate(energies):
        residuals.append((energ - v[i])**2)
    return [energies, residuals]

#calculate all atom-atom energies outside of the cutoff. self-interaction in molecule is excluded
def get_energy(epsilon, sigma, atoms, m):
    E = [];
    for i, p1 in enumerate(atoms):
        e = []

        #calculate p1's interactions with p2 thru pn, then increment p1 to p2, dealing with p3 thru pn
        for j, p2 in enumerate(atoms[(i+1):]):

            r = distance_between_coords(p1, p2)
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

## begin main ##

#load files
files = open('filenames.txt'); dirs = [line.strip('\n') for line in files.readlines()]
files.close()

#initial guesses
sot_eps = 1e-7; sot_sig = 1

# mixing rules for LJ parameters
mix_eps = sot_eps * carbon_eps; mix_sig = (sot_sig + carbon_sig) / 2

#read energies
energy_file = open('energies.txt')
v = [float(line.strip('\n').split('\t')[2]) for line in energy_file.readlines()]
energy_file.close()

m = 1; [energies, residuals] = get_residuals(dirs, v, mix_eps, mix_sig, m)
