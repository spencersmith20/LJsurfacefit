#written by spencer smith for wilmerlab nov 2018 - jan 2019

#imports
from ase import Atoms; from ase.io import read, write
import numpy as np; from math import cos, sin, pow
from os import mkdir, chdir; from random import uniform
from shutil import copy

#file input
print('enter your molecule .xyz file'); f1 = input()
print('enter your surface .cif file'); f2 = input()

#surface displacement for molecule c.o.m.
mol = read(f1); print('enter desired buffer size'); buffer = 7
dists = mol.positions.max(0) - mol.positions.min(0)

#create systemname
dir_name = f1.split('.')[0] + '_' + f2.split('.')[0]

#get cif file info
surf = read(f2); omega = surf.get_volume()
[a, b, c, alpha, beta, gamma] = surf.get_cell_lengths_and_angles()
fractional = surf.get_scaled_positions()

#convert from degrees to radians
alpha = np.deg2rad(alpha); beta = np.deg2rad(beta); gamma = np.deg2rad(gamma)

#this is the desired box size in x,y,z directions
box = dists + 2 * buffer + mol.get_center_of_mass()

#create mapping matrix -- convert fractional into cartesian unit cell
M = np.array([[a, b*cos(gamma), c*cos(beta)],
        [0, b*sin(gamma), c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma)],
        [0, 0, omega/(a*b*sin(gamma))]])

#generate cartesian coordinates for unit cell
unit_coords = []; unit_types = []
for l in fractional:
    l.reshape(3,1); unit_coords.append(np.dot(M, l)); unit_types.append('C')

#start at molecular center of mass so they're in line
[x0, y0, z0] = unit_coords[0]; [x, y, z] = [x0, y0, z0]; i = 0
coords = []; types = []; vec = [unit_coords[1][0] - x, unit_coords[1][1] - y];

#each loop run completes a column of linked unit cells
while(y >= (y0 - box[1])):
    while (x <= (x0 + box[1]) and x >= x0):

        #draw the items in the unit cell (this could be turned into a loop)
        coords.append([x, y, z]); types.append(unit_types[0])
        [x, y] = map(sum, zip([x, y], vec))

        coords.append([x, y, z]); types.append(unit_types[1])
        [x, y] = map(sum, zip([x, y], [-k for k in vec]))
        x = x + a

    #move to the next column of unit cells
    y = y - b*sin(gamma) + vec[1]; x = x0

    #reverse vector to make mirror image for other half of hexagon
    vec = [-k for k in vec]

y_min = min(coords[1])

#fix periodicity in y-axis
for coord in coords:
    if (coord[1] == y_min):
        coords.remove(coord)

#create an Atoms object of the surface unit cell
surface = Atoms(positions=coords, symbols=types)

#find box dimensions
x_box = max(coords[0]) - min(coords[0]); z_box = box[2]
y_box = max(coords[1]) - min(coords[1]) + b*sin(gamma) - vec[1]

#write to file to be used by simulation
box_file = open('box.txt', 'w')
box_file.write(str(x_box) + '\n' + str(y_box) + '\n' + str(z_box))
box_file.close()

#line up COMs, "zero" the system
shift = surface.get_center_of_mass() - mol.get_center_of_mass()

#shift each coordinate (ASE translate functionality is bad, IMO)
translated = np.empty(mol.positions.shape)
for i in range(len(mol.positions)):
    translated[i] = np.add(mol.positions[i], shift)
mol.set_positions(translated)

#generate configurations
print('files to be generated:'); num = int(input())

#create directory, save name of directory for individual .xyz files
mkdir(dir_name)

#copy the bash script, input file, slurm submission file
for file in ['setslurm_periodic.sh', 'periodic.inp', 'cp2k-per-template.slurm', 'scrape.sh', 'box.txt']:
    copy(file, dir_name)

#change directories
chdir(dir_name); print('\n.xyz files in directory: ' + dir_name)

#outfile for writing filename generated (used by bash) & for logging radii
filenames = open('filenames.txt', 'w'); distances = open('distances.txt', 'w')

#list of random floats for radial distances and angles
ang_max = 360; r_list = [uniform(1,6) for i in range(num)]
[ang_phi, ang_psi, ang_theta] = [[uniform(0, ang_max) for i in range(num)] for j in range(3)]

for j in range(1, num+1):

    #get values for this iteration
    phi = ang_phi[j-1]; psi = ang_psi[j-1]; theta = ang_theta[j-1]; r = r_list[j-1]

    #transform the molecule according to random angles
    mol.euler_rotate(center='COM', phi=phi, psi=psi, theta=theta)

    #include distance from center of mass to z minimum once rotated
    r = r + (mol.get_center_of_mass()[2] - mol.positions[:, 2].min())

    #create final xyz file
    filename3 = dir_name + str(j) + '.xyz'

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

    #attach molecule to surface, output to xyz file in "sysname_j.xyz"
    output.extend(out_mol)
    write(filename=(dir_name + str(j) + '.xyz'), images=output, format='xyz')

#close output files
for file in [filenames, distances]:
    file.close()
