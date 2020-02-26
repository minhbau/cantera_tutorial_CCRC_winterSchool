"""
A freely-propagating, premixed flat flame.
"""

################################################################################
# Load modules
################################################################################
import cantera as ct
import numpy as np
import matplotlib.pylab as plt

def calcHRR(T, P, Yk):
    HRR = np.zeros((T.shape))
    Yspc = np.zeros(nSpcs)
    for i in range(T.shape[0]):
        Yspc[:] = Yk[:,i]
        gas.TPY = T[i], P, Yspc
        reactor = ct.IdealGasReactor(gas)
        HRR[i] = -np.dot(reactor.kinetics.net_production_rates,\
        reactor.thermo.partial_molar_enthalpies)
    return HRR
# HRR = calcHRR(flame.T,flame.P,flame.Y)

################################################################################
# function to compute flame thickness
################################################################################
def compute_thickness(flame_obj):
    dTdx = np.gradient(flame_obj.T, flame_obj.grid)
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
nSpcs = gas.n_species
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

# compute the heat release rate
HRR = calcHRR(f.T,f.P,f.Y)

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

#### Temperature Plot
plt.figure()
plt.plot(f.grid*100, f.T, '-o')
plt.xlabel('Distance (cm)')
plt.ylabel('Temperature (K)')
plt.savefig('Tmulti_Profile.pdf')
plt.show()

plt.figure()
plt.plot(f.grid*100, HRR, '-o')
plt.xlabel('Distance (cm)')
plt.ylabel(r"HRR [$J/m^3/s$]")
plt.savefig('HRRmulti_Profile.pdf', dpi=300)
plt.show()
