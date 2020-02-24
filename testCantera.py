#!/usr/bin/python3
# simple test to ensure that cantera importeds and runs properly
import cantera as ct

# to work with numpy array
import numpy as np
# for plotting
import matplotlib.pyplot as plt
import matplotlib

# check install path
print(ct.__path__)
# check cantera version
print(ct.__version__)

# print out the cantera modules
dir(cantera) 


# set the mechanism name
gas = ct.Solution('gri30.xml')

# get helps
help(gas.species_index)
help(gas.__class__.T)

# print out the gas classess/functions
dir(gas)


# ------------------------------------------------------------------
print(gas.n_reactions) # number of reactions
print(gas.n_species) # number of species
print(gas.species_names) # returns a list of species names
print(gas.species_names[0]) # returns the name of selected species


# ====================================================================
#! Global reaction, stoichiometric CH4/air mixture
# CH4 + 2(O2 + 3.76 N2) --> CO2 + 2H2O + 2*3.76 N2
# ====================================================================
# intial condition:
Tin = 1000.0 # K
Pin = 40.0 # atm
composition = 'CH4:1.0, O2:2,N2:7.52' 

# set up the reactor
gas.TPX = Tin, Pin*ct.one_atm, composition


print(gas.X) # returns an array of moles fractions
print(gas.Y) # returns an array of mass fractions
print(gas.X[0])  # returns the mole fraction of selected species

print(gas['CH4'].X)
print(gas['O2'].X)
print(gas['N2'].X)
print(gas['CH4'].partial_molar_cp)


# ====================================================================
#! select reactor: UV fixed specific internal energy and specific volume
# ====================================================================
gas.TPX = Tin, Pin*ct.one_atm, composition
gas.equilibrate('UV')
Teq = gas.T # get equilibrium (adiabatic) temperature
print(Teq)


# ====================================================================
#! select reactor: HP fixed specific enthalpy and pressure 
# ====================================================================
gas.TPX = Tin, Pin*ct.one_atm, composition
gas.equilibrate('HP')
print(gas.T)
# ====================================================================


# ====================================================================
# find the corresponding adibatic temperature with the UV reactor
# by using a loop over a  range of equivalence ratio, phi
# ====================================================================
npoints = 20
phi = np.linspace(0.1, 2.0, npoints)
T_ad = np.zeros(npoints)
for i in range(npoints):
    gas.TP = 300, ct.one_atm
    gas.set_equivalence_ratio(phi[i], 'CH4', 'O2:1, N2:3.76')
    gas.equilibrate('UV')
    T_ad[i] = gas.T
    print(phi[i], T_ad[i])

print(T_ad)

# save the result to a file
np.savetxt('save_phi_Tad.dat', np.transpose([phi, T_ad]), delimiter=' ', fmt='%.2f %12.5e')

# plotting
plt.plot(phi,T_ad)
plt.xlabel('phi')
plt.ylabel('T_ad (K)')
# show the plot
plt.show()
# ====================================================================

