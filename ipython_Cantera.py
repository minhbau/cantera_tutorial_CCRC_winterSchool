# ====================================================================
#! download the tutorial 
# ====================================================================
cd $HOME
git clone https://github.com/minhbau/cantera_tutorial_CCRC_winterSchool.git
cd cantera_tutorial_CCRC_winterSchool

# activate cantera enviroment
source /opt/cantera-2.4/bin/setup_cantera

# run an example test
python3 example.py

# ====================================================================
#! activate interactive ipython enviroment
ipython3
# run a python script inside ipython enviroment
run ./example.py

# get helps in IPython shell
gas.species_index?

# exit ipython enviroment
exit()
# ====================================================================


# ====================================================================
# load cantera
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


# ====================================================================
# set the mechanism name
gas = ct.Solution('gri30.xml')
# the mechanism, gri30.xml
# GRI-Mech 3.0 mechanism consists of 325 reactions that involve 53 species 
# to model natural gas combustion.

# get helps
help(gas.species_index)
help(gas.__class__.T)

# print out the gas classess/functions
dir(gas)


# ====================================================================
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

# set up a state
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
# intial condition:
Tin = 1000.0 # K
Pin = 40.0 # atm
composition = 'CH4:1.0, O2:2,N2:7.52' 

gas.TPX = Tin, Pin*ct.one_atm, composition
gas.equilibrate('UV')
Teq = gas.T # get equilibrium (adiabatic) temperature
print(Teq)


# ====================================================================
#! select reactor: HP fixed specific enthalpy and pressure 
# ====================================================================
# intial condition:
Tin = 1000.0 # K
Pin = 40.0 # atm
composition = 'CH4:1.0, O2:2,N2:7.52' 

gas.TPX = Tin, Pin*ct.one_atm, composition
gas.equilibrate('HP')
print(gas.T)
# ====================================================================


# ====================================================================
#! find the indices of just those reactions which convert CO into CO2:
# ====================================================================
gas = ct.Solution('gri30.xml')
II = [i for i,r in enumerate(gas.reactions())
    if 'CO' in r.reactants and 'CO2' in r.products]
for i in II:
    print(i, gas.reaction(i).equation)
gas.net_rates_of_progress[II]

# ===========================================================================


# ===========================================================================
#! Equation of state (EoS):
# ===========================================================================
# P =  m/V *RT/MW = rho*T*R/MW = rho*T*R_specific
# Where P in Pa=J/m^3, m in kg;  V m^3, T (K),
# ideal, or universal gas constant R = 8.314 J/(K·mol), 

# intial condition:
Tin = 1000.0 # K
Pin = 40.0 # atm
composition = 'CH4:1.0, O2:2,N2:7.52' 

# set up a state
gas.TPX = Tin, Pin*ct.one_atm, composition
print(P/ct.one_atm,ct.gas_constant,gas.T,gas.density,gas.mean_molecular_weight)

# compute P by using EoS
P = ct.gas_constant*gas.density*gas.T/gas.mean_molecular_weight
print(P/ct.one_atm)
# ====================================================================

