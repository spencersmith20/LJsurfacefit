from numpy import linspace, power, random, inf, asarray, absolute
from pandas import read_csv; from math import sqrt
from ase import Atoms; from os import chdir
from sympy import symbols, diff, factor, Matrix; from numpy.linalg import inv

#calculated / look-up values
sotolon_vac = -86.96520401259681; graphene_vac = -305.49689755541885
carbon_sig = 3.55; carbon_eps = 0.066 / 627.509 #kcal/mol -> Hartree
sotolon_pervac = -86.96759558341823; graphene_pervac = -7.39805875741705

#get the distance between two atom positions
def distance_between_coords(p1, p2):
    return sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)

#for each configuration of the molecule, calculate total E and find residual
def get_obj_func(dir_list, v, m, f):

    #define initial residual to be added for each .xyz file
    residual = 0; print(str(len(dir_list)))

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

        #list for iteration of energies. m is number of bodies interaction (currently dealing with sphere approximation)
        atoms = mol.get_center_of_mass() + surf.positions

        #get the sum squared error for this set of parameters
        residual = residual + (1/2)*(v[j] - get_energy(f, atoms, m))**2

        print('j = ' + str(j))

        #go into bulk directory
        chdir('..')

    return residual

#calculate all atom-atom energies outside of the cutoff. self-interaction in molecule is excluded
def get_energy(f, atoms, m):

    E = 0
    for i, p1 in enumerate(atoms):

        #calculate p1's interactions with p2 thru pn, then increment p1 to p2, dealing with p3 thru pn
        for j, p2 in enumerate(atoms[(i+1):]):
            dist = distance_between_coords(p1, p2)

            if dist > 6: #cut-off
                continue;

            if i < m: # deal with sotolon-carbon interaction ... will need to deal with an array of parameters
                E = E + f.subs(r, dist)
            else: # deal with carbon-carbon interaction
                E = E + f.subs([(r, dist), (eps, carbon_eps), (sig, carbon_sig)])
    return E

## begin main ##

#load files
files = open('filenames.txt'); dirs = [line.strip('\n') for line in files.readlines()]
files.close()

#declare symbols (variables ... x0 epsilon, x1 sigma), max attempts
max_attempts = 300
eps, sig, r = symbols('eps sig r')

#declare function of variables and radius
f = -4 * eps * ((sig/r)**12 - (sig/r)**6)

#read energies
energy_file = open('energies.txt')
v = [float(line.strip('\n').split('\t')[2]) for line in energy_file.readlines()]
energy_file.close(); m = 1

#find objective function to be minimized -- this takes forever
obj = get_obj_func(dirs, v, m, f)
obj = factor(obj)

print('we out here we past the objective function yeet yeet yet')

#create in jacobian matrix
Jt = Matrix([f.diff(eps), f.diff(sig)])
J = Jt.T

#define initial guesses for (epsilon, sotolon) parameter pair -- stored as tuple
params = [(1e-7, 2.5)]

for k in range(max_attempts):

    #extract parameters
    print('uh beginning of loop'); x = params[k]

    #substitute parameters into jacobian, transpose, and residual
    #J_sub, Jt_sub, obj_sub = [m.subs([(eps, x[0]), (sig, x[1])]) for m in [J, Jt, obj]]

    #calculate moore-penrose pseudo-inverse
    #MPPI = ((Jt_sub * J_sub)**-1) * Jt_sub; print(MPPI)

    MPPI = ((Jt * J)**-1) * Jt; print(MPPI)
    delta = MPPI * obj; print(delta)
    delta.sub([(eps, x[0]), (sig, x[1])])

    #calcualte delta for params, increment & append
    new_params = (x[0] + delta[0], x[1] + delta[1]); params.append(new_params)

    print('we made it'); print(params)

    new_resid = obj.subs([(eps, new_params[0]), (sig, new_params[1])])
    print('new residual calcualted')

    #check for parameter convergence
    conv_criteria = absolute((resid - new_resid)/resid)
    print(conv_criteria)

    if conv_criteria < 0.0001:  #convergence achieved
        print('epsilon: ' + str(new_params[0]) + '\tsigma: ' + str(new_params[1]))
        break
    print('uh end of loop')

# use mixing rules to extract sotolon sigma and epsilon at the end
