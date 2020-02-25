# ========================================================================================
#  Load modules
# ========================================================================================
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

# ========================================================================================
# Mixture composition from equivalence ratio
# ========================================================================================
def composition_CxHyOz_Air(gas, fuel_species, phi):
    iFuel = gas.species_index(fuel_species)
    iO2 = gas.species_index('O2')
    iN2 = gas.species_index('N2')
    C_atoms = gas.n_atoms(fuel_species, 'C')
    H_atoms = gas.n_atoms(fuel_species, 'H')
    O_atoms = gas.n_atoms(fuel_species, 'O')
    stoich_O2 = C_atoms + H_atoms/4.0 - O_atoms/2.0
    X = np.zeros(gas.n_species)
    X[iFuel] = phi
    X[iO2] = stoich_O2
    X[iN2] = 3.76*stoich_O2
    return X
# ========================================================================================
# Heat release rate [W/m^3]
# ========================================================================================
def heat_release_rate(reactor):
    wdot_molar = reactor.kinetics.net_production_rates # [kmol/m^3/s]
    h_molar = reactor.thermo.partial_molar_enthalpies  # [J/kmol]
    nsp = reactor.thermo.n_species
    hrr = 0.0
    hrr = - np.dot(wdot_molar,h_molar) # faster

    # for i in range(nsp): # loop slower
    #     hrr += -(wdot_molar[i]*h_molar[i])

    return hrr

# ========================================================================================
#! Reactor time-stepping
# ========================================================================================
# set the temperature range to optimize the time step
LowT   =  600 
InterT =  750 
HighT  =  1000
# set time-stepping parameters
time_max_veryLowT   =  1e+4    # if  T =< 600 K 
time_max_LowT       =  1e+2    # if 600 < T =< 750 K 
# time_max_InterT     =  1e+0    # if 750 < T =< 1000 K 
time_max_InterT     =  50.0*1e-3    # if 750 < T =< 1000 K 
time_max_HighT      =  10.0*1e-3 # if T > 1000 K   

# set suitable time_step for each T range 
# time_step = 1e-6
time_step_veryLowT  =  5e-5 # if T   < 600 K 
time_step_LowT      =  5e-6 # if 600 < T =< 750 K  
time_step_InterT    =  5e-7 # if 750 < T =< 1000 K 
time_step_HighT     =  1e-7 # if T > 1000 K   

def RunReactor(mechFile, fuel_species, temp, press, phi, time_step, time_max):
    gas = ct.Solution(mechFile)
    gas.TPX = temp, press, composition_CxHyOz_Air(gas, fuel_species, phi)

    gas2 = ct.Solution(mechFile)
    gas2.TPX = temp, press, composition_CxHyOz_Air(gas, fuel_species, phi)
    gas2.equilibrate('UV')

    # set the reactor and
    r = ct.IdealGasReactor(gas)
    sim = ct.ReactorNet([r])

    # initialize some parameters
    Teq = gas2.T
    Told = r.T
    dT = 0.0
    dTmax = 1.0
    HHRmax = 0.0
    iOH = gas.species_index('OH')
    max_Y_OH = 0.0
    HRRmax_final = 0.0
    time_HRRmax_final = 0.0
    time_max_Y_OH = 0.0
    T_HTHR = 1200 # (K) - T ~1200 K  at the onset of the main combustion
    T_ITHR = 900 # (K) - T ~950 K  at the onset of H2O2 decomposition

    time = 0.0
    timeig_dT500 = 0.0
    time_maxdT = 0.0
    time_max_HRR = 0.0
    ignited = False
    ncounts = 0

    # set output file name
    stem, ext = mechFile.split('.')
    out_file = 'saveIgnData/'+ stem + '_p' + str('{:.0f}'.format(P/ct.one_atm)) + '_phi' + str(phi) + '_T' + str(temp) + '.dat'
    zone = 'T' + str(temp) + 'K_phi' + str(phi) + '_p' + str(press/ct.one_atm)
    
    if write_out_file == 1:
        file = open(out_file,'w')
        # header info for the output file
        # file.write('#t(s) HRR(W/m3) Y_OH T(K) P(atm)\n ')
        file.write('variables = "time (ms)" "HRR (W/m3)" "Y_OH" "T (K)" "P (atm)"\n ')

        # write the  mass fraction of all the species
        for j in range(gas.n_species): 
            file.write('"%s"\n '%(gas.species_names[j]))
        file.write('zone t = "%s" \n'%(zone))

    while True:
        time += time_step
        ncounts += ncounts
        sim.advance(time)
        hrr = heat_release_rate(r)
        if r.T > T_ITHR and hrr > HHRmax :
            HHRmax = hrr
            time_max_HRR = time
        if r.T > T_HTHR and r.T < (T_HTHR+1):
            HRRmax_final = hrr
            time_HRRmax_final = time
        if (hrr > HRRmax_final) and (r.T > T_HTHR):
            HRRmax_final = hrr
            time_HRRmax_final = time
        if r.T > T_ITHR and r.T < (T_ITHR+1):
            max_Y_OH = gas.Y[iOH]
            time_max_Y_OH = time
        if (gas.Y[iOH] > max_Y_OH) and (r.T > T_ITHR):
            max_Y_OH = gas.Y[iOH]
            time_max_Y_OH = time*0.95
        # if gas.Y[iOH] > max_Y_OH:
            # max_Y_OH = gas.Y[iOH]
            # time_max_Y_OH = time
        dT = r.T-Told
        if dT > dTmax:
            dTmax = dT
            time_maxdT = time
        Told = r.T
        if ignited == False and r.T > (temp+450.0):
            timeig_dT500 = time
            ignited = True
        # output to file
        if ((ncounts%nStepSkip) ==0) and (write_out_file == 1):

            # output time, HRR, Y_OH, T, P
            file.write('%12.5e %10.3e %10.3e %10.3f %10.3f'%(sim.time*1e3, hrr, gas.Y[iOH], r.T, r.thermo.P/ct.one_atm))

            # output the species mass fraction
            for j in range(gas.n_species): file.write(' %10.3e'%(gas.Y[j]))
            file.write('\n')
        if (time > time_max) or (r.T >= (Teq-10.0)):
            actual_time_stop = time 
            break 
    if write_out_file == 1:
        file.close
    #print('Ignition time: %g'%timeig)
    return time_max_HRR*1e3, time_HRRmax_final*1e3, actual_time_stop*1e3, time_max_Y_OH*1e3, time_maxdT*1e3, timeig_dT500*1e3


# ========================================================================================
# Generate list of values
# ========================================================================================
def Generate_list(min_val, max_val, delta):
    n_pts = (max_val-min_val)/delta + 1
    tmp_list = []
    for i in range( int(n_pts) ):
        tmp_list.append( min_val + i*delta )
    if tmp_list[-1] < max_val:
        tmp_list.append( max_val )
    return tmp_list


# ========================================================================================
#!  MAIN Program
# ========================================================================================
if __name__ == "__main__":

    # define the gas and fuel
    mechFile = "dme39.xml"  
    fuel_species =  'CH3OCH3'

#-------------------------------------------------------- 
    #! choose to write write_out_file to track  (time, HRR, OH, T, P, species...)
    write_out_file = 1 # 1-write, 0-not write 

    # skips a number of steps/interval not write the data, big for low T, nStepSkip = ~ 1e4 
    nStepSkip = 1e3
#-------------------------------------------------------- 

    #! set range of pressure values in atm
    min_val = 40 
    max_val = 40 
    step = 1.0
    P_values = Generate_list(min_val, max_val, step)
    # convert to Pascals
    for i in range(len(P_values)):
        P_values[i] *= ct.one_atm
    print('\nInitial pressure values:',P_values)
    
    #! set range of equivalence ratio values
    max_val = 1.00
    min_val = 1.00
    step =    0.01
    phi_values = Generate_list(max_val, min_val, -step)
    # phi_values = Generate_list(min_val, max_val, step)
    print('\nEquivalence ratio values:',phi_values)

    #! set range of temperature values
    min_val = 700.0   
    max_val = 1200.0  
    step = 100.0
    T_values = Generate_list(min_val, max_val, step)
    print('\nInitial temperature values:',T_values)
    print('\n')

    n_Tvals = len(T_values)
    ignition_time = np.zeros(n_Tvals)

    # # main loop
    # fig1 = plt.figure(1)
    # # loop over pressure
    for P in P_values:
        # loop over equivalence ratio
        for phi in phi_values:
            name0D_ignition_results = 'saveIgnData/0Dignition_DME_' + 'phi' + str('{:.2f}'.format(phi)) + '_p' + str('{:.0f}'.format(P/ct.one_atm))+'atm'  
            ignitionResults = name0D_ignition_results + '.dat'
            ignitionResults_multiCriteria = name0D_ignition_results + '_multiCriteria' + '.dat'

            # file2 = open(ignitionResults_multiCriteria,'w')
            # # header info for the output file
            # file2.write('variables = "T(K)" "time_max_HRR (ms)" "time_HRRmax_final" "time_max_Y_OH" "time_maxdT" "timeig_dT500" "actual_time_stop" \n' )

            # header info for the output file
            file3 = open(ignitionResults,'w')
            file3.write('variables = "T(K)" "time_HRRmax_final"\n' )
            file3.write('zone t = "%s"\n'%(ignitionResults) )
            # loop over temperature
            for i in range(n_Tvals):
                # set optimal values for different T range
                if T_values[i] < LowT:  
                   time_max  = time_max_veryLowT
                   time_step = time_step_veryLowT
                if T_values[i] <= InterT  and T_values[i] >= LowT:
                   time_max  = time_max_LowT
                   time_step = time_step_LowT
                if T_values[i] <= HighT  and T_values[i] > InterT:
                   time_max  = time_max_InterT
                   time_step = time_step_InterT
                if T_values[i] > HighT:  
                   time_max  = time_max_HighT
                   time_step = time_step_HighT
                
                # run the reactor
                time_max_HRR, time_HRRmax_final,  actual_time_stop, time_max_Y_OH, time_maxdT, timeig_dT500 = \
                RunReactor(mechFile, fuel_species, T_values[i], P, phi, time_step, time_max)
                print('Ignition time for P =',P/ct.one_atm, ', phi =',phi, ', T =',T_values[i], '==>', time_max_HRR, time_HRRmax_final, time_max_Y_OH, time_maxdT, timeig_dT500,  actual_time_stop, 'ms')
                #time_ig = timeig_dT500
                time_ig = time_max_HRR
                ignition_time[i] = time_ig


                # output to file
                # file2.write('%10.2f %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e \n'%(T_values[i], time_max_HRR, \
                # time_HRRmax_final,  actual_time_stop, time_max_Y_OH, time_maxdT, timeig_dT500))
                file3.write('%10.2f %12.5e \n'%(T_values[i], time_HRRmax_final))

            # close the file
            # file2.close
            file3.close

            print('\nignition_time (ms):', ignition_time)
            # plots
            plt.figure(1)
            p1, = plt.plot(T_values,ignition_time,'-o')
            plt.xlabel('Temperature [K]',fontsize=14)
            plt.ylabel('t_ig [s]',fontsize=14)
            plt.grid(True)
    # show plots
    plt.show()