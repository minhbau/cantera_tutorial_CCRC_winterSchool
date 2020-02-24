#!/usr/bin/python3
import cantera as ct
import numpy as np # to work with numpy array
import matplotlib.pyplot as plt # for plotting
import matplotlib

# ====================================================================
# find the corresponding adiabatic temperature with the UV reactor
# by using a loop over a  range of equivalence ratio, phi
# ====================================================================


# ------------------------------------------------------------------------
#! Generate list functions
# ------------------------------------------------------------------------
def Generate_list(min_val, max_val, delta):
    n_pts = (max_val-min_val)/delta + 1
    tmp_list = []
    for i in range( int(n_pts) ):
        tmp_list.append( min_val + i*delta )
    if tmp_list[-1] < max_val:
        tmp_list.append( max_val )
    return np.asarray(tmp_list)

# Main program
if __name__ == "__main__":
    # set the mechanism name
    gas = ct.Solution('gri30.xml')

    phi_min = 0.1
    phi_max = 2.0
    delta = 0.1

    phi_arr = Generate_list(phi_min, phi_max, delta)
    npoints = len(phi_arr)

    T_ad = np.zeros(npoints)
    for i in range(npoints):
        gas.TP = 300, ct.one_atm
        gas.set_equivalence_ratio(phi_arr[i], 'CH4', 'O2:1, N2:3.76')
        gas.equilibrate('UV')
        T_ad[i] = gas.T
        print(phi_arr[i], T_ad[i])

    print(T_ad)

    # save the result to a file
    np.savetxt('save_phi_Tad.dat', np.transpose([phi_arr, T_ad]), delimiter=' ', fmt='%.2f %12.5e')

    # plotting
    plt.plot(phi_arr,T_ad, '-o')
    plt.xlabel('phi_arr')
    plt.ylabel('T_ad (K)')
    # show the plot
    plt.show()
    # ====================================================================

