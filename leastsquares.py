import numpy as np; import sympy as sp; from numpy.linalg import inv
from pandas import read_csv; from math import sqrt; from ase import Atoms
from os import chdir; import matplotlib.pyplot as plt

#calculated / look-up values
sotolon_vac = -86.96520401259681; graphene_vac = -305.49689755541885; methane_vac = -8.07599867375029
carbon_sig = 3.55; carbon_eps = 0.066 / 627.509 #kcal/mol -> Hartree
sotolon_pervac = -86.96759558341823; graphene_pervac = -7.39805875741705

#get the distance between two atom positions
def distance_between_coords(p1, p2):
    return sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)

#for each configuration of the molecule, calculate total E then find residual
def get_obj_funcs(dir_list, v, b, f):

    #define initial residual to be added for each .xyz file
    residual = []; surf_E = 0; energ = []

    for j, dir in enumerate(dir_list):

        #move into particular folder
        chdir(dir)

        #get & separate coordinates
        df = read_csv(dir + '.xyz', header=2, names=['type', 'x', 'y', 'z'], sep='\s+')
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

        #only need to calculate surface energy once, b/c identical from trial-to-trial
        if j == 0:
            surf_E = get_surface_energy(f, surf.positions)

        #calculate energy for each config.
        energ.append(get_energy(f, [mol.get_center_of_mass()], surf.positions) + surf_E)
        rp1 = (1/2)*(v[j] - energ[j])**2
        residual.append(rp1); #print(rp1)

        #go into bulk directory
        chdir('..')

    return residual

#calculate the vacuum energy of the surface with given LJ params -- only do this once!
def get_surface_energy(f, surf_list):
    E = 0
    for i, p1 in enumerate(surf_list):

        #calculate p1's interactions with p2 thru pn, then increment p1 to p2
        for j, p2 in enumerate(surf_list[(i+1):]):
            dist = distance_between_coords(p1, p2)

            if dist > 6: #cut-off
                continue;

            #calculate surface energy if less than cut-off
            e = f.subs([(r, dist), (sig, carbon_sig), (eps, carbon_eps)])

            E = E + e
    return E

#calculate all atom-atom energies outside of the cutoff. self-interaction in molecule is excluded
def get_energy(f, mol_list, surf_list):
    E = 0
    for i, p1 in enumerate(mol_list):

        #calculate p1's interactions with p2 thru pn, then increment p1 to p2, dealing with p3 thru pn
        for j, p2 in enumerate(surf_list):
            dist = distance_between_coords(p1, p2)

            if dist > 6: #cut-off
                continue;

            #calculate the energy if less than cut-off
            E = E + f.subs(r, dist)

    return E

#load files
files = open('filenames.txt'); dirs = [line.strip('\n') for line in files.readlines()]
files.close()

#declare symbols (variables ... x0 epsilon, x1 sigma), max attempts
max_attempts = 300
eps, sig, r = sp.symbols("eps sig r")

#declare function of variables and radius
f = 4 * eps * ((sig/r)**12 - (sig/r)**6); b = 1

#read energies and radii
energy_file = open('energies.txt'); v = []; rad = []; lines = energy_file.readlines()
r_list = [float(l.strip('\n').split('\t')[1]) for l in lines]
e_list = [(float(l.strip('\n').split('\t')[2]) - methane_vac - graphene_vac) / 27.2114 for l in lines]

st_dev = np.std(np.array(e_list).astype(np.float64))
avg = np.average(np.array(e_list).astype(np.float64))

plt.plot(r_list, e_list, 'go'); plt.show()

#filter out non-LJ regions
energy_file.close(); directories = []; rad = []; v = []
for i, V in enumerate(e_list):
    if (avg - st_dev*.6 <= V <= avg + st_dev*.4):
        rad.append(r_list[i]); v.append(V); directories.append(dirs[i])
plt.plot(rad, v, 'bo'); plt.show()

#size of Jacobian matrix -- m is # objective functions, n is # params
m = len(v); n = 2

#find objective function to be minimized -- this is extremely time limiting
obj = sp.Matrix(1, m, sp.factor(get_obj_funcs(directories, v, b, f))).T

#fill in Jacobian matrix
J = sp.Matrix(m, n, [sp.factor(sp.diff(o, s)) for o in obj for s in [eps, sig]]); Jt = J.T

#define initial guesses for (epsilon, sotolon) parameter pair -- stored as tuple
params = [(1e-5, 5.3)]

#calculate moore-penrose pseudo-inverse
a = Jt * J

for k in range(max_attempts):

    x = params[k]
    ap = a.subs([(eps, x[0]), (sig, x[1])]); a_inv = ap**-1

    #calculate Moore-Penrose pseudoinverse, make substitutions, calcualte delta
    Jtp, obj_sub = [ma.subs([(eps, x[0]), (sig, x[1])]) for ma in [Jt, obj]]
    MPPI = a_inv * Jtp; delta = MPPI * obj_sub; new_params = []

    #keep from going negative, could implement a line search here
    for j, X in enumerate(x):
        if (X - delta[j]*.7) > 0:
            new_params.append(x[j] - delta[j]*.7)
        else:
            new_params.append(x[j] + delta[j]*.7)
    params.append(new_params)

    newobj_sub = obj.subs([(eps, new_params[0]), (sig, new_params[1])])
    print('new residual calcualted')
    print(x); print(new_params)

    #conv_criteria = obj_sub - newobj_sub; conv_criteria = [conv_criteria[i] / r2 for i, r2 in enumerate(obj_sub)]
    conv_criteria = [abs((X - new_params[i])/X) for i, X in enumerate(x)]
    print(conv_criteria)

    if all(np.absolute(conv) < 0.05 for conv in conv_criteria):
        break

# use mixing rules to extract sotolon sigma and epsilon at the end
