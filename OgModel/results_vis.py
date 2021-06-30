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

def cutAndTransformToArray(l, ind):
    l = l[ind:]
    l = np.array(l)
    return l

def convert_to_fraction(array):
    array = array / np.sum(array)

    return array 


def parametrize(parameters):

    overflow = False
    c = parameters

    IDL.run("cd, '/home/hmorenom/SSW_Files/OgModel'  \n")
    IDL.run("restore,'fast_wind.save'\n", stdout=True)

    height = IDL.HEIGHT
    bfield = IDL.BFIELD
    alfven = IDL.ALFVEN

    original_mass = IDL.mass
    original_temp = IDL.temp
    original_vel = IDL.VELOCITY 

    print(height.shape)
    print(original_vel.shape)

    mass = []

    temp = []

    velocity = []

    for r in height:
        
        m = (c[0]/(r**c[1])) + (c[5]/(r**c[6]))
        
        T = c[2]*(((c[3]*r)**c[7])/(r-c[4]))*math.exp(-(c[3]/(r-c[4]))**3)
        
        mass.append(m)
        temp.append(T)

    ind_20k = 0

    for index, elem in enumerate(temp):
        if elem > 2e4:
            ind_20k = index
            break
    
    r_0 = height[ind_20k]

    print('r_0 = {}'.format(r_0))

    x_min =c[8]/(1.00925 ** c[9]) + c[10]/(1.00925 ** c[11])

    print('x_min = {}'.format(x_min))

    for r in height:
        
        x = c[8]/(r ** c[9]) + c[10]/(r ** c[11])
        v = x_min - x
        velocity.append(v)
    

    ind_vel = 0

    for index, elem in enumerate(velocity):
        if elem > 1:
            ind_vel = index
            break

    index = max(ind_20k, ind_vel) 

    print('temp_0 before cutting = {}'.format(temp[0]))
    print('vel_0 before cutting = {}'.format(velocity[0]))
    

    height = cutAndTransformToArray(height, index)
    mass = cutAndTransformToArray(mass, index)
    temp = cutAndTransformToArray(temp, index)
    velocity = cutAndTransformToArray(velocity, index)
    bfield = cutAndTransformToArray(bfield, index)
    alfven = cutAndTransformToArray(alfven, index)

    print('temp_0 after cutting = {}'.format(temp[0]))
    print('vel_0 after cutting = {}'.format(velocity[0]))
    

    original_temp = cutAndTransformToArray(original_temp, index)
    original_mass = cutAndTransformToArray(original_mass, index)
    original_vel = cutAndTransformToArray(original_vel, index)
    
    IDL.mass = mass
    IDL.temp = temp
    IDL.height = height
    IDL.velocity = velocity
    IDL.bfield = bfield
    IDL.alfven = alfven
    
    #IDL.run("ratio = proton_dens(alog10(temp), /hydrogen, abund_file=\'/usr/local/ssw/packages/chianti/dbase/abundance/sun_photospheric_2011_caffau.abund\', ioneq_file = \'/usr/local/ssw/packages/chianti/dbase/ioneq/chianti.ioneq\')\n")

    IDL.run("save, units, height, mass, velocity, temp, bfield, alfven, filename= 'paramvis_fast_wind.save', /verb ")
    #ratio = IDL.RATIO


    return np.array([height, original_mass, original_temp, original_vel, mass, temp, velocity]), index
    #IDL.run("save, units, height, mass, velocity, temp, bfield, alfven, filename= 'param_fast_wind.save', /verb ")

    #logging.info('Saved param_fast_wind.save file.')


def get_predictions():

    IDL.run("cd, '/home/hmorenom/SSW_Files/OgModel/temp_vis'  \n")
    IDL.run("restore,'pred_c.save'\n", stdout=True)

    height = IDL.HHH_USED
    temp = IDL.TTT_USED
    density = IDL.DDD_USED
    velocity = IDL.VVV_USED
    
    return np.array([height, temp, density, velocity])

def get_carbon_ion():


    IDL.run("cd, '/home/hmorenom/SSW_Files/FastWindData'  \n")
    IDL.run("restore,'fast_wind_measurements_ace.save', /v\n", stdout=True)
    IDL.run("measured_ion = carbon(1,*)")

    measured_ion = IDL.MEASURED_ION

    measured_ion = convert_to_fraction(measured_ion)

    IDL.run("cd, '/home/hmorenom/SSW_Files/OgModel/temp_vis'  \n")
    IDL.run("restore,'pred_c.save'\n", stdout=True)
    IDL.run("predicted_ion = ioneq_evol(*,-1)")

    predicted_ion = IDL.PREDICTED_ION

    return np.column_stack((np.array(measured_ion), np.array(predicted_ion))).T

def get_oxygen_ion():


    IDL.run("cd, '/home/hmorenom/SSW_Files/FastWindData'  \n")
    IDL.run("restore,'fast_wind_measurements_ace.save', /v\n", stdout=True)
    IDL.run("measured_ion = oxygen(1,*)")

    measured_ion = IDL.MEASURED_ION

    measured_ion = convert_to_fraction(measured_ion)

    IDL.run("cd, '/home/hmorenom/SSW_Files/OgModel/temp_vis'  \n")
    IDL.run("restore,'pred_o.save'\n", stdout=True)
    IDL.run("predicted_ion = ioneq_evol(*,-1)")

    predicted_ion = IDL.PREDICTED_ION

    return np.column_stack((np.array(measured_ion), np.array(predicted_ion))).T

def get_iron_ion():


    IDL.run("cd, '/home/hmorenom/SSW_Files/FastWindData'  \n")
    IDL.run("restore,'fast_wind_measurements_ace.save', /v\n", stdout=True)
    IDL.run("measured_ion = iron(1,*)")


    measured_ion = IDL.MEASURED_ION

    measured_ion = convert_to_fraction(measured_ion)

    IDL.run("cd, '/home/hmorenom/SSW_Files/OgModel/temp_vis'  \n")
    IDL.run("restore,'pred_fe.save'\n", stdout=True)
    IDL.run("predicted_ion = ioneq_evol(*,-1)")

    predicted_ion = IDL.PREDICTED_ION

    return np.column_stack((np.array(measured_ion), np.array(predicted_ion))).T
    
def plot_temp(height, temp):

    fig_dens, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
    ax1.plot(height, temp, color='blue',
             label='Final Parametrized Temperature')
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

def plot_vel(height, vel):

    fig_dens, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
    ax1.plot(height, vel, color='blue',
             label='Final Parametrized Velocity')
    ax1.hlines(y=c[-2], xmin = height[0], xmax = height[-1], color='r', linestyle='--')
    ax1.set_xlabel('Distance from the Sun [solar radii] (log scale)')
    ax1.set_ylabel('Velocity (log scale)')
    ax1.set(title='Velocity vs. Distance from the Sun')
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


def plot_density(height, density):

    fig, (ax1) = plt.subplots(1, 1, figsize=(10, 10))
    ax1.plot(height, density, color='blue',
             label='Final Parametrized Electron Density')
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

    ions_fe_cut = [ions_fe[0,5:17], ions_fe[1,5:17]]
    Y =  np.arange(5,17)
    ax3.bar(Y - 0.20, ions_fe_cut[0], color = 'b', width = 0.4, label='Measured Fractions')
    ax3.bar(Y + 0.20, ions_fe_cut[1], color = 'g', width = 0.4, label='Predicted Fractions')
    
    ax3.set_xlabel('Ionization Stages [+e]')
    ax3.set_ylabel('Fraction of Total Number of Fe Atoms')	
    ax3.set(title='Fraction of Measured and Predicted Iron Ions')
    ax3.legend()
    ax3.set_xticks(Y)


c = [2.36196e-17, 3.63, 2000000.0, 0.4, 0.75, 2.42e-15, 21.87, 0.7128]

c =[4.2000000000000005e-17, 2.3213428125, 2000000.0, 0.361, 0.748125, 2e-15, 30, 0.8, 2000, 600, 687.66667, 0.4]

c = [4.4e-17, 2.43, 2000000.0, 0.324, 0.75, 1.8e-15, 39.93, 0.72, 2420.0, 393.65999999999997, 582.33333, 0.4]

starting_c = [4e-17, 3, 2000000.0, 0.4, 0.75, 2e-15, 30, 0.8, 2000, 600, 687.66667, 0.4]

starting_c = [4e-17, 3, 2e6, 0.4, 0.75, 2e-15, 30, 0.8, 2000, 600, 750, 0.4]



values, index = parametrize(c)

ssw.run('run_wind_vis')

#start_vals, start_index = parametrize(starting_c)

pred = get_predictions()

#print(pred)

plot_temp(pred[0], pred[1])

plot_density(pred[0], pred[2])

plot_vel(pred[0], pred[3])


c_ions = get_carbon_ion()
o_ions = get_oxygen_ion()
fe_ions = get_iron_ion()

print(c_ions)


#print(fe_ions)

plot_ions(c_ions, o_ions, fe_ions)

plt.show()
