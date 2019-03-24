#written by spencer smith | wilmer lab | october 2018
from ase import Atoms; from ase.io import read, write
import numpy as np; from os import mkdir, chdir; from random import uniform
from shutil import copy

#read filenames
print('enter name of molecule .xyz file:'); f1=input()
print('enter name of surface .xyz file:'); f2=input()

#read in surface and molecule, get COMs for alignment
mol = read(f1); surface = read(f2)

#shift coms
shift = surface.get_center_of_mass() - mol.get_center_of_mass()

#shift each coordinate
translated = np.empty(mol.positions.shape)
for i in range(len(mol.positions)):
    translated[i] = np.add(mol.positions[i], shift)
mol.set_positions(translated)

#create directory and get system name
dir_name = f1.split('.')[0] + '_' + f2.split('.')[0]
mkdir(dir_name)

#num is desired number of stochastic files generated from combination
print('\nfiles to be generated:'); num = int(input())

#copy files into directory
for file in ['setslurm_non.sh', 'energy-geoopt.inp', 'cp2k-template.slurm', 'scrape.sh']:
    copy(file, dir_name)

#change directories
chdir(dir_name); print('\n.xyz files in directory: ' + dir_name)

#outfile for writing filename generated (used by bash), writing z dists
filenames = open('filenames.txt', 'w'); distances = open('distances.txt', 'w')

#create configurations
ang_max = 360; r_list = [uniform(1, 6) for i in range(num)]
[ang_phi, ang_psi, ang_theta] = [[uniform(0, ang_max) for i in range(num)] for j in range(3)]

for j in range(1, num+1):

    #generate random values for phi, and theta. beginning of loop
    phi = ang_phi[j-1]; psi = ang_psi[j-1]; theta = ang_theta[j-1]; r = r_list[j-1]

    #transform the molecule according to random angles
    mol.euler_rotate(center='COM', phi=phi, psi=psi, theta=theta)

    #include distance from center of mass to z minimum once rotated
    r = r + (mol.get_center_of_mass()[2] - mol.positions[:, 2].min())

    #format angles and radial distance for .xyz header
    for thing in [phi, theta, psi, r]:
        thing = format(thing, '.3f')

    #write to outputfiles
    filenames.write(dir_name + str(j) + '\n')
    distances.write(str(j) + '\t' + str(r) +'\n')

    #save original info for next iteration of loop
    out_mol = mol.copy(); output = surface.copy()

    #translate molecule in z direction
    translated = np.empty(mol.positions.shape)
    shift = np.array([0, 0, r])

    for i in range(len(mol.positions)):
        translated[i] = np.add(mol.positions[i], shift)
    out_mol.set_positions(translated)

    #attach molecule to surface, output to xyz file in "sysname{$j}.xyz"
    output.extend(out_mol)
    write(filename=(dir_name + str(j) + '.xyz'), images=output, format='xyz')

filenames.close(); distances.close()
