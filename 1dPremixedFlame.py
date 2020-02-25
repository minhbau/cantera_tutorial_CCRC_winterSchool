"""
A freely-propagating, premixed flat flame.
"""

################################################################################
# Load modules
################################################################################
import cantera as ct
import numpy as np


################################################################################
# function to compute flame thickness
################################################################################
def compute_thickness(flame_obj):
    # print(flame_obj.T)
    # print(flame_obj.flame.grid)
    dTdx = np.gradient(flame_obj.T, flame_obj.flame.grid)
    max_dTdx = np.max(dTdx)
    thickness = (flame_obj.T[-1] - flame_obj.T[0])/max_dTdx
    return thickness 


################################################################################
#
# main code
#
################################################################################
# Simulation parameters
p = 1.0*ct.one_atm  # pressure [Pa]
Tin = 300.0  # unburned gas temperature [K]
reactants = 'CH4:1, O2:2, N2:7.52'  # premixed gas composition
eq_ratio = 1.0
width = 0.05  # domain length [m]
loglevel = 0  # amount of diagnostic output (0 to 8)


# IdealGasMix object used to compute mixture properties, set to the state of the
# upstream fuel-air mixture
gas = ct.Solution('gri30.xml')
gas.TPX = Tin, p, reactants

# Set equivalence ratio
gas.set_equivalence_ratio( phi=eq_ratio, fuel={'CH4':1}, oxidizer={'O2':1.0, 'N2':3.76} )


# Set up flame object
f = ct.FreeFlame(gas, width=0.05)
f.set_refine_criteria(ratio=3, slope=0.05, curve=0.03)
#f.show_solution()


# Turn on/off radiation
f.radiation_enabled = False


################################################################################
# Solve with mixture-averaged transport model
f.transport_model = 'Mix'
f.solve(loglevel=loglevel, auto=True)
# Save solution
f.save('sol_freeflame.xml', 'mix', 'solution with mixture-averaged transport')
#f.show_solution()
print('Mixture-averaged flamespeed = {0:7f} m/s'.format(f.u[0]))
# Write the velocity, temperature, density, and mass/mole fractions to a CSV file
f.write_csv('sol_freeflame_mix.csv', species='X', quiet=False)
# # Determine the flame thickness
# flame_thickness = compute_thickness(f)
# print('Flame thickness = {0:4e} m'.format(flame_thickness))


################################################################################
# Solve with multi-component transport properties
f.transport_model = 'Multi'
f.solve(loglevel) # don't use 'auto' on subsequent solves
#f.show_solution()
print('Multicomponent flamespeed = {0:7f} m/s'.format(f.u[0]))
# Save solution
f.save('sol_freeflame.xml', 'multi', 'solution with multicomponent transport')
# Write the velocity, temperature, density, and mass/mole fractions to a CSV file
f.write_csv('sol_freeflame_multi.csv', species='X', quiet=False)
# Determine the flame thickness
flame_thickness = compute_thickness(f)
print('Flame thickness = {0:4e} m'.format(flame_thickness))


################################################################################
# Solve with multi-component transport properties and thermal diffusion
f.transport_model = 'Multi'
f.soret_enabled = True
f.solve(loglevel) # don't use 'auto' on subsequent solves
#f.show_solution()
print('Multicomponent with Soret effect flamespeed = {0:7f} m/s'.format(f.u[0]))
# Save solution
f.save('sol_freeflame.xml', 'multi_soret', 'solution with multicomponent transport and Soret effect')
# Write the velocity, temperature, density, and mass/mole fractions to a CSV file
f.write_csv('sol_freeflame_multi_soret.csv', species='X', quiet=False)
# # Determine the flame thickness
# flame_thickness = compute_thickness(f)
# print('Flame thickness = {0:4e} m'.format(flame_thickness))

