# Load cantera modules
import cantera as ct

def soundSpeed(gas):
    import cantera as ct
    gamma = gas.cp/gas.cv
    Rgas = ct.gas_constant/gas.mean_molecular_weight
    sound_speed = (gamma * Rgas * gas.T)**0.5
    return sound_speed

# Main program
if __name__ == "__main__":
    # define the gas
    gas = ct.Solution('dme39.xml')
    # intial condition:
    Tin = 300.0 # K
    Pin = 40.0 # atm
    composition = 'CH4:1.0, O2:2,N2:7.52' 

    # set up a state
    gas.TPX = Tin, Pin*ct.one_atm, composition
    print('Universal gas constant = %.2f [J/kmol/K]' %(ct.gas_constant))
    print('Gas density = %12.5e [kg/m^3]' %(gas.density))
    print('Mean molecular weight = %.2f[kg/kmol].' %(gas.mean_molecular_weight))

    thermal_diffusivity = gas.thermal_conductivity/(gas.density*gas.cp)
    kinematic_viscosity = gas.viscosity/gas.density

    print('thermal_diffusivity = %12.5e' %(thermal_diffusivity))
    print('kinematic_viscosity = %12.5e [m^2/s]' %(kinematic_viscosity))

    soundSpeed_freshGas = soundSpeed(gas)
    print('soundSpeed_freshGas = %.2f [m/s].' %(soundSpeed_freshGas))

    # get the equilibrium condition conscant U & V 
    gas.equilibrate('UV')
    Teq = gas.T
    Peq = gas.P
    soundSpeed_eq = soundSpeed(gas)

    print('Equilibrium (adiabatic) temperature = %.2f [K]' %(Teq))
    print('Equilibrium pressure = %.2f [atm].' %(Peq/ct.one_atm))
    print('Equilibrium soundSpeed_eq = %.2f [m/s].' %(soundSpeed_eq))

