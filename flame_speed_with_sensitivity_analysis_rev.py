from __future__ import print_function
from __future__ import division

import cantera as ct
import numpy as np
import matplotlib.pylab as plt
import pandas as pd

print("Running Cantera Version: " + str(ct.__version__))

def compute_thickness(flame_obj):
    dTdx = np.gradient(flame_obj.T, flame_obj.flame.grid)
    max_dTdx = np.max(dTdx)
    thickness = (flame_obj.T[-1] - flame_obj.T[0])/max_dTdx
    return thickness 

def calcHRR(T, P, Yk):
    HRR = np.zeros((T.shape))
    Yspc = np.zeros(nSpcs)
    for i in range(flame.T.shape[0]):
        Yspc[:] = Yk[:,i]
        gas.TPY = T[i], P, Yspc
        reactor = ct.IdealGasReactor(gas)
        HRR[i] = -np.dot(reactor.kinetics.net_production_rates,\
        reactor.thermo.partial_molar_enthalpies)
    return HRR
# HRR = calcHRR(flame.T,flame.P,flame.Y)

# ========================================================================================
#! Define the reactant conditions, gas mixture and kinetic mechanism associated with the gas
# ========================================================================================
#Inlet Temperature in Kelvin and Inlet Pressure in Pascals
#In this case we are setting the inlet T and P to room temperature conditions
To = 500
Po = 101325

#Define the gas-mixutre and kinetics
#In this case, we are choosing a GRI3.0 gas
gas = ct.Solution('gri30.xml')

spcsName = gas.species_names
nSpcs = gas.n_species

# Create a stoichiometric CH4/Air premixed mixture 
gas.set_equivalence_ratio(1.0, 'CH4', {'O2':1.0, 'N2':3.76})
gas.TP = To, Po


# ========================================================================================
##! Define flame simulation conditions
# ========================================================================================
# Domain width in metres
width = 0.014

# Create the flame object
flame = ct.FreeFlame(gas, width=width)

# Define tolerances for the solver
flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)

# Define logging level
loglevel = 1



# solve
flame.solve(loglevel=loglevel, auto=True)
Su0 = flame.u[0]
print("Flame Speed is: {:.2f} cm/s".format(Su0*100))

# compute the heat release rate
HRR = calcHRR(flame.T,flame.P,flame.Y)

# Determine the flame thickness
flame_thickness = compute_thickness(flame)
print('Flame thickness = {0:4e} m'.format(flame_thickness))

#Note that the variable Su0 will also be used downsteam in the sensitivity analysis

# ========================================================================================
###! Plot figures
# ========================================================================================
# Import plotting modules and define plotting preference
# import matplotlib.pylab as plt

plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.figsize'] = (8,6)

# Get the best of both ggplot and seaborn
plt.style.use('ggplot')
plt.style.use('seaborn-deep')

plt.rcParams['figure.autolayout'] = True

# #### Temperature Plot
# plt.figure()
# plt.plot(flame.grid*100, flame.T, '-o')
# plt.xlabel('Distance (cm)')
# plt.ylabel('Temperature (K)')
# plt.savefig('TemperatureProfile', dpi=300)
# plt.show()


fig = plt.figure()
ax = fig.add_subplot()

lns1 = ax.plot(flame.grid*100, flame.T, '-', label = 'T')
ax2 = ax.twinx()
lns2 = ax2.plot(flame.grid*100, HRR, '-', label = 'HRR')

# add legends together
lns = lns1+lns2
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc=0)

ax.grid()
ax.set_xlabel("Distance (cm)")
ax.set_ylabel(r" T (K)")
ax2.set_ylabel(r"HRR [$J/m^3/s$]")
# ax2.set_ylim(0, 35)
# ax.set_ylim(-20,100)
plt.show()

# plt.savefig('1Dflame_T_HRR', dpi=300,
plt.savefig('1Dflame_T_HRR.pdf',bbox_inches='tight')
plt.show()

# ========================================================================================
####! Major species' plot
# ========================================================================================
"""
# To plot species, we first have to identify the index of the species in the array
# For this, cut & paste the following lines and run in a new cell to get the index
for i, specie in enumerate(gas.species()):
    print(str(i) + '. ' + str(specie))
"""

# Extract concentration data
X_CH4 = flame.X[13]
X_CO2 = flame.X[15]
X_H2O = flame.X[5]

plt.figure()

plt.plot(flame.grid*100, X_CH4, '-o', label=r'$CH_{4}$')
plt.plot(flame.grid*100, X_CO2, '-s', label=r'$CO_{2}$')
plt.plot(flame.grid*100, X_H2O, '-<', label=r'$H_{2}O$')

plt.legend(loc=2)
plt.xlabel('Distance (cm)')
plt.ylabel('MoleFractions')
# plt.savefig('MajorSpecies', dpi=300)
plt.savefig('MajorSpecies.pdf')

plt.show()


## Sensitivity Analysis
# See which reactions effect the flame speed the most


# ========================================================================================
#!Import a data frame module. This simplifies the code
# import pandas as pd
# ========================================================================================
# Create a dataframe to store sensitivity-analysis data
sensitivities = pd.DataFrame(data=[], index=gas.reaction_equations(range(gas.n_reactions)))
### Compute sensitivities

# Set the value of the perturbation
dk = 1e-2

# Create an empty column to store the sensitivities data
sensitivities["baseCase"] = ""


for m in range(gas.n_reactions):
    gas.set_multiplier(1.0) # reset all multipliers                                                                     
    gas.set_multiplier(1+dk, m) # perturb reaction m   
    
    # Always force loglevel=0 for this
    # Make sure the grid is not refined, otherwise it won't strictly 
    # be a small perturbation analysis
    flame.solve(loglevel=0, refine_grid=False)
    
    # The new flame speed
    Su = flame.u[0]
    
    sensitivities["baseCase"][m] = (Su-Su0)/(Su0*dk)

# This step is essential, otherwise the mechanism will have been altered
gas.set_multiplier(1.0)


###! Make plots

# Reaction mechanisms can contains thousands of elementary steps. Choose a threshold
# to see only the top few
threshold = 0.03

firstColumn = sensitivities.columns[0]

# For plotting, collect only those steps that are above the threshold
# Otherwise, the y-axis gets crowded and illegible
sensitivitiesSubset = sensitivities[sensitivities[firstColumn].abs() > threshold]
indicesMeetingThreshold = sensitivitiesSubset[firstColumn].abs().sort_values(ascending=False).index
sensitivitiesSubset.loc[indicesMeetingThreshold].plot.barh(title="Sensitivities for GRI 3.0",
                                                          legend=None)
plt.gca().invert_yaxis()

plt.rcParams.update({'axes.labelsize': 20})
plt.xlabel(r'Sensitivity: $\frac{\partial\:\ln{S_{u}}}{\partial\:\ln{k}}$');

# Uncomment the following to save the plot. A higher than usual resolution (dpi) helps
# plt.savefig('sensitivityPlot', dpi=300)
plt.savefig('sensitivityPlot.pdf')
plt.show()

