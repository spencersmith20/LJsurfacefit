#written by spencer smith | wilmerlab, feb 2019
import matplotlib.pyplot as plt
from numpy import linspace, power, random, inf, asarray, ones
from pandas import read_csv
from lmfit import Model, Parameters
from decimal import Decimal

files = []
print('enter path to each result file (press enter when done)')

while(True):

    #entering a newline moves to the rest of the script
    filename = input()
    if (filename is ''):
        break

     #add file to list of inputs to be scanned
    files.append(filename)

files=['energies_NP.txt']

#concatenate file data
cat_file = open('big.dat', 'w')
for file in files:
    f = open(file)

    for line in f:
        cat_file.write(line)
    f.close()
cat_file.close()

#calculated / look-up values
sotolon_vac = -86.96520401259681; graphene_vac = -305.49689755541885
methane_vac = -8.07599867375029
carbon_sig = 3.55; carbon_eps = 0.066 / 627.509 #kcal/mol -> hartree

#initial guesses
sot_eps = 1e-7; sot_sig = 1

#get data from each file
df = read_csv('big.dat', header=None, names=['num', 'distance', 'energy'], sep='\t')
print(df.to_string(index=False))

#subtract out vacuum energies
r_data = df['distance']; v_data = df['energy']
v_data = [(dat - (sotolon_vac + graphene_vac)) for dat in v_data]

#puts smaller weight on values < 2
w = ones((len(r_data)))
for i in range(len(r_data)):
    if r_data[i] < 4.4:
        w[i] = 0

#define functions
def lennard_jones(r, eps, sig):
    return 4 * eps * ((sig/r)**12 - (sig/r)**6)

#define model
p = Parameters()
v_min = min(v_data); r_corr = r_data[v_data.index(v_min)]
p.add('eps', value=.2, min=0); p.add('sig', value=r_corr, min=0)
model = Model(lennard_jones)

#get parameters
result = model.fit(v_data, p, r=r_data, weights=w)
params= result.best_values

#display results
print(result.fit_report())

#get epsilon into scientific notation
ep = '%.2E' % Decimal(params['eps'])
print('fit epsilon is ' + str(ep) + ', fit sigma is ' + str(params['sig']))

sot_eps = params['eps'] / carbon_eps; sot_sig = 2*params['sig'] - carbon_sig
s_ep = '%.2E' % Decimal(sot_eps)

print('sotolon epsilon is ' + str(s_ep) + ', sotolon sigma is ' + str(sot_sig))

#get curve fit
fit_r = linspace(1.2, 6, 100)
fit_v = lennard_jones(fit_r, params['eps'], params['sig'])

#plot line and data
plt.plot(r_data, v_data, 'bo', label='electronic calc')
#plt.plot(fit_r, fit_v, 'k-', label='LJ fit')

#plt.ylim(-.05, .05)
plt.show()
