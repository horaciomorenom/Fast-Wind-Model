#!/opt/anaconda3/bin/python
import os
from datetime import datetime
from idlpy import *
import numpy as np
import matplotlib.pyplot as plt
import hissw
import math

ssw = hissw.Environment(
    ssw_packages=['packages/chianti/'], ssw_paths=['chianti'])

#parameters#

c1 = 4e-17
c2 = 2.43
c3 = 3e6
c4 = 0.324
c5 = 0.54675
c6 = 2e-15
c7 = 33
c8 = 0.8

c1 = 4e-17
c2 = 3
c3 = 2e6
c4 = 0.4
c5 = 0.75  # 0.739
c6 = 2e-15
c7 = 30
c8 = 0.8

c1 = 2.871e-17
c2 = 5.8461513
c3 = 2e6
c4 = 0.36
c5 = 0.7425
c6 = 2.2e-15
c7 = 24.3
c8 = 0.8

def cutAndTransformToArray(l, ind):
    l = l[ind:]
    l = np.array(l)
    return l


def parametrize(parameters):

    overflow = False
    c = parameters

    IDL.run("cd, '/home/hmorenom/SSW_Files/OgModel'  \n")
    IDL.run("restore,'fast_wind.save'\n", stdout=True)

    height = IDL.HEIGHT
    velocity = IDL.VELOCITY
    bfield = IDL.BFIELD
    alfven = IDL.ALFVEN

    original_mass = IDL.mass
    original_temp = IDL.temp

    mass = []

    temp = []

    for r in height:
        try:
            m = (c[0]/(r**c[1])) + (c[5]/(r**c[6]))
        except OverflowError:
            overflow = True
        try:
            T = c[2]*(((c[3]*r)**c[7])/(r-c[4]))*math.exp(-(c[3]/(r-c[4]))**3)
        except OverflowError:
            overflow = True
        mass.append(m)
        temp.append(T)

    ind_20k = 0

    ind_vel = 0

    for index, elem in enumerate(temp):
        if elem > 2e4:
            ind_20k = index
            break

    for index, elem in enumerate(velocity):
        if elem > 1:
            ind_vel = index
            break

    index = max(ind_20k, ind_vel)

    height = cutAndTransformToArray(height, index)
    mass = cutAndTransformToArray(mass, index)
    temp = cutAndTransformToArray(temp, index)
    velocity = cutAndTransformToArray(velocity, index)
    bfield = cutAndTransformToArray(bfield, index)
    alfven = cutAndTransformToArray(alfven, index)

    original_temp = cutAndTransformToArray(original_temp, index)
    original_mass = cutAndTransformToArray(original_mass, index)
    
    IDL.mass = mass
    IDL.temp = temp
    IDL.height = height
    IDL.velocity = velocity
    IDL.bfield = bfield
    IDL.alfven = alfven
    
    #IDL.run("ratio = proton_dens(alog10(temp), /hydrogen, abund_file=\'/usr/local/ssw/packages/chianti/dbase/abundance/sun_photospheric_2011_caffau.abund\', ioneq_file = \'/usr/local/ssw/packages/chianti/dbase/ioneq/chianti.ioneq\')\n")

    IDL.run("save, units, height, mass, velocity, temp, bfield, alfven, filename= 'param_fast_wind.save', /verb ")
    #ratio = IDL.RATIO


    return np.array([height, original_mass, original_temp, mass, temp]), index
    #IDL.run("save, units, height, mass, velocity, temp, bfield, alfven, filename= 'param_fast_wind.save', /verb ")

    #logging.info('Saved param_fast_wind.save file.')


def get_density():
    
    IDL.run("cd, '/home/hmorenom/SSW_Files/OgModel/density'  \n")
    IDL.run("restore,'density_original.save'\n", stdout=True)

    original_density = IDL.DENSITY

    original_density = np.array(original_density)

    IDL.run("cd, '/home/hmorenom/SSW_Files/OgModel/density'  \n")
    IDL.run("restore,'density_parameters.save'\n", stdout=True)

    parameters_density = IDL.DENSITY
    parameters_density = np.array(parameters_density)

    IDL.run("cd, '/home/hmorenom/SSW_Files/OgModel/density'  \n")
    IDL.run("restore,'density_starting_params.save'\n", stdout=True)

    start_parameters_density = IDL.DENSITY
    start_parameters_density = np.array(parameters_density)

    return np.array([original_density, parameters_density, start_parameters_density])

def get_carbon_ion():


    IDL.run("cd, '/home/hmorenom/SSW_Files/FastWindData'  \n")
    IDL.run("restore,'fast_wind_measurements.save', /v\n", stdout=True)
    IDL.run("measured_ion = carbon")

    measured_ion = IDL.MEASURED_ION

    IDL.run("cd, '/home/hmorenom/SSW_Files/OgModel/temp'  \n")
    IDL.run("restore,'pred_c.save'\n", stdout=True)
    IDL.run("predicted_ion = ioneq_evol(*,-1)")

    predicted_ion = IDL.PREDICTED_ION

    return np.array([measured_ion, predicted_ion])

def get_oxygen_ion():


    IDL.run("cd, '/home/hmorenom/SSW_Files/FastWindData'  \n")
    IDL.run("restore,'fast_wind_measurements.save', /v\n", stdout=True)
    IDL.run("measured_ion = oxygen")

    measured_ion = IDL.MEASURED_ION

    IDL.run("cd, '/home/hmorenom/SSW_Files/OgModel/temp'  \n")
    IDL.run("restore,'pred_o.save'\n", stdout=True)
    IDL.run("predicted_ion = ioneq_evol(*,-1)")

    predicted_ion = IDL.PREDICTED_ION

    return np.array([measured_ion, predicted_ion])

def get_iron_ion():


    IDL.run("cd, '/home/hmorenom/SSW_Files/FastWindData'  \n")
    IDL.run("restore,'fast_wind_measurements.save', /v\n", stdout=True)
    IDL.run("measured_ion = iron")

    measured_ion = IDL.MEASURED_ION

    IDL.run("cd, '/home/hmorenom/SSW_Files/OgModel/temp'  \n")
    IDL.run("restore,'pred_fe.save'\n", stdout=True)
    IDL.run("predicted_ion = ioneq_evol(*,-1)")

    predicted_ion = IDL.PREDICTED_ION

    return np.array([measured_ion, predicted_ion])
    
def plot_temp(values_end,  values_start):

    fig_dens, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
    ax1.plot(values_end[0], values_end[2], color='green', label='Theoretical Temperature')
    ax1.plot(values_end[0], values_end[4], color='blue',
             label='Final Parametrized Temperature')
    ax1.plot(values_start[0], values_start[4], 'g--', color='red',
             label='Initial Parametrized Temperature')
    ax1.set_xlabel('Distance from the Sun [solar radii] (log scale)')
    ax1.set_ylabel('Temperature [K] (log scale)')
    ax1.set(title='Temperature vs. Distance from the Sun')
    ax1.grid()
    #ax2.plot(values[0], values[1], color='green', label='Original Mass')
    #ax2.plot(values[0], values[3], color='blue',
    #         label='First Parametrized Mass')
    #ax2.set_xlabel('Distance (log(solar radii))')
    #ax2.set_ylabel('Mass')
    #ax2.set(title='Mass')
    ax1.legend()
    ax1.set_yscale('log')
    ax1.set_xscale('log')

    #ax2.set_yscale('log')
    #ax2.set_xscale('log')


def plot_densities(height, dens):

    fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
    ax1.plot(height, dens[0], color='green', label='Theoretical Electron Density')
    ax1.plot(height, dens[1], color='blue',
             label='Final Parametrized Electron Density')
    ax1.plot(height, dens[2], 'g--', color='red',
             label='Starting Parametrized Electron Density')
    ax1.set_xlabel('Distance from the Sun[solar radii] (log scale) ')
    ax1.set_ylabel('Electron Density (log scale)')
    ax1.set(title='Electron Density vs. Distance from the Sun')
    ax1.grid()	

    #ax2.plot(values[0], values[1], color='green', label='Original Mass')
    #ax2.plot(values[0], values[3], color='blue',
    #         label='First Parametrized Mass')
    #ax2.set_xlabel('Distance (log(solar radii))')
    #ax2.set_ylabel('Mass')
    #ax2.set(title='Mass')

    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.legend()

    #ax2.set_yscale('log')
    #ax2.set_xscale('log')

def plot_ion(ions):
    
    fig_ion, (ax1) = plt.subplots(1,1, figsize=(26,15))
    X = np.arange(7)
    ax1.bar(X - 0.20, ions[0], color = 'b', width = 0.4, label='Measured Fractions')
    ax1.bar(X + 0.20, ions[1], color = 'g', width = 0.4, label='Predicted Fractions')
    
    ax1.set_xlabel('Ionization Stages [+e]')
    ax1.set_ylabel('Fraction of Total Number of C Atoms')	
    ax1.set(title='Fraction of Measured and Predicted Carbon Ions')
    ax1.legend()

def plot_ions(ions_c, ions_o, ions_fe):
    
    fig_ion, (ax1,ax2,ax3) = plt.subplots(1,3, figsize=(20,8))
    ions_c_cut = [ions_c[0,3:7], ions_c[1,3:7]]
    X = np.arange(3,7)
    ax1.bar(X - 0.20, ions_c_cut[0], color = 'b', width = 0.4, label='Measured Fractions')
    ax1.bar(X + 0.20, ions_c_cut[1], color = 'g', width = 0.4, label='Predicted Fractions')
    
    ax1.set_xlabel('Ionization Stages [+e]')
    ax1.set_ylabel('Fraction of Total Number of C Atoms')	
    ax1.set(title='Fraction of Measured and Predicted Carbon Ions')
    ax1.legend()
    ax1.set_xticks(np.arange(min(X), max(X)+1, 1.0))    

    Z = np.arange(4,9)
    ions_o_cut = [ions_o[0,4:9], ions_o[1,4:9]]
    ax2.bar(Z - 0.20, ions_o_cut[0], color = 'b', width = 0.4, label='Measured Fractions')
    ax2.bar(Z + 0.20, ions_o_cut[1], color = 'g', width = 0.4, label='Predicted Fractions')
    
    ax2.set_xlabel('Ionization Stages [+e]')
    ax2.set_ylabel('Fraction of Total Number of O Atoms')	
    ax2.set(title='Fraction of Measured and Predicted Oxygen Ions')
    ax2.legend()
    ax2.set_xticks(Z)

    ions_fe_cut = [ions_fe[0,7:17], ions_fe[1,7:17]]
    print(ions_fe_cut)
    Y =  np.arange(7,17)
    ax3.bar(Y - 0.20, ions_fe_cut[0], color = 'b', width = 0.4, label='Measured Fractions')
    ax3.bar(Y + 0.20, ions_fe_cut[1], color = 'g', width = 0.4, label='Predicted Fractions')
    
    ax3.set_xlabel('Ionization Stages [+e]')
    ax3.set_ylabel('Fraction of Total Number of Fe Atoms')	
    ax3.set(title='Fraction of Measured and Predicted Iron Ions')
    ax3.legend()
    ax3.set_xticks(Y)

c = [c1, c2, c3, c4, c5, c6, c7, c8]

c = [2.36196e-17, 3.63, 2000000.0, 0.4, 0.75, 2.42e-15, 21.87, 0.7128]

c =[2.1892047693750008e-17, 4.00684079296875, 2000000.0, 0.4, 0.7129687499999999, 2.4805000000000004e-15, 22.97716875, 0.7867978284374999]

starting_c = [4e-17, 4, 2e6, 0.5, 0.5, 2e-15, 30, 0.2]



values, index = parametrize(c)

#ssw.run('run_wind')

start_vals, start_index = parametrize(starting_c)


plot_temp(values, start_vals)

#densities = get_density()

#print(densities[0].shape)

#densities[0] = densities[0][index:]
 
#densities[1] = densities[1, index:]

#plot_densities(values[0], densities)

c_ions = get_carbon_ion()
o_ions = get_oxygen_ion()
fe_ions = get_iron_ion()



#print(fe_ions)

plot_ions(c_ions, o_ions, fe_ions)
#plt.plot([0,1,2,3], [0,1,2,3])


c2 = 2.925
c3 = 2e6
c4 = 0.35626482
c5 = 0.75371245
c6 = 2e-15
c7 = 27.27
c8 = 0.796

plt.show()
