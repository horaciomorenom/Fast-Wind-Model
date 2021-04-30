#!/usr/bin/env python3
import os
from datetime import datetime
from idlpy import *
import numpy as np
import matplotlib.pyplot as plt
import hissw
import math
import logging

ssw = hissw.Environment(
    ssw_packages=['packages/chianti/'], ssw_paths=['chianti'])


def createDirectory(path):

    now = datetime.now()
    date_time = now.strftime("%m/%d/%Y $H:$M:$S")
    path = path + "/model_{}".format(now)
    os.mkdir(path)
    return path


def get_observations():

    IDL.run("cd, 'FastWindData' \n")
    IDL.run("restore,'fast_wind_measurements.save'\n", stdout=True)
    c_obs = IDL.carbon
    n_obs = IDL.nitrogen
    o_obs = IDL.oxygen
    ne_obs = IDL.neon
    mg_obs = IDL.magnesium
    si_obs = IDL.silicon
    s_obs = IDL.sulphur
    fe_obs = IDL.iron

    observations = [c_obs, n_obs, o_obs, ne_obs, mg_obs, si_obs, s_obs, fe_obs]

    return observations


# takes in a list and an index and returns a numpy array of all elements in the list after the index
def cutAndTransformToArray(l, index):

    l = l[index:]
    l = np.array(l)
    return(l)


def parametrize(parameters):

    overflow = False
    c = parameters

    IDL.run("cd, '/home/hmorenom/SSW_Files/OgModel'  \n")
    IDL.run("restore,'fast_wind.save'\n", stdout=True)

    height = IDL.HEIGHT
    velocity = IDL.VELOCITY
    bfield = IDL.BFIELD
    alfven = IDL.ALFVEN

    mass = []

    temp = []

    for r in height:
        try:
            m = (c[0]/(r**c[1])) + (c[5]/(r**c[6]))
        except OverflowError:
            logging.info(
                'Overflow Error in mass calculation, skipping set of parameters.')
            overflow = True
        try:
            T = c[2]*(((c[3]*r)**c[7])/(r-c[4]))*math.exp(-(c[3]/(r-c[4]))**3)
        except OverflowError:
            logging.info(
                'Overflow Error in temperature calculation, skipping set of parameters.')
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

    IDL.mass = mass
    IDL.temp = temp
    IDL.height = height
    IDL.velocity = velocity
    IDL.bfield = bfield
    IDL.alfven = alfven

    IDL.run("save, units, height, mass, velocity, temp, bfield, alfven, filename= 'param_fast_wind.save', /verb ")

    logging.info('Saved param_fast_wind.save file.')

    return overflow


# takes array with filenames and returns array that includes all predictions
def getPredictions(files):

    path = '/home/hmorenom/SSW_Files/OgModel/temp/'

    IDL.run("cd, '" + path + "'  \n")

    predictions = []

    for x in filenames:
        IDL.run("restore, '"+x + "'\n")
        IDL.run("pred = ioneq_evol(*,-1) \n")
        pred = IDL.pred
        predictions.append(pred)

    return np.array(predictions)


# takes array that includes obbservations and array of filenames
def chiSquared(observations, predictions):

    chi_2_errors = []

    for x in range(len(predictions)):
        error = 0
        for i in range(len(predictions[x])):
            if observations[x][i] != 0:
                error = error + \
                    ((predictions[x][i]-observations[x][i])
                     ** 2/observations[x][i])
            else:
                if predictions[x][i] > 1e-3:
                    logging.info(
                        "Element {}, ionization stage {} has predted value > 10^3".format(x, i))
                    error = error + ((predictions[x][i]-1e-3)**2/1e-3)
        chi_2_errors.append(error)

    return chi_2_errors


# takes parameter values and observation ioneq values, returns an array containing (params, chi2)
def run_iteration(params, obs, filenames):

    overflow = parametrize(params)
    if overflow != True:
        try:
            logging.info('Running run_wind program')
            ssw.run('run_wind')
        except Exception as e:
            print('Error when running run_wind, error: {}'.format(e))
            logging.info("Failed to run run_wind program. Error: {}".format(e))
        predictions = getPredictions(filenames)

        chi_2_val = chiSquared(obs, predictions)

        return [params, chi_2_val, overflow]
    else:
        return [0, 0, overflow]


def randomize_parameters(parameters, factor=0.1):  # randomizes a random parameter by 10%

    params = parameters[:]
    index = 2

    while (index == 2):
        index = np.random.randint(len(params))

    sign = np.random.choice([-1, 1])

    params[index] += sign*factor*params[index]

    return params


def apply_Metropolis(current_iteration, new_iteration):

    u = np.random.rand()  # generate random float between 0 and 1
    A = min(
        1, np.sum(np.array(new_iteration[1])) / np.sum(np.array(current_iteration)[1]))

    if (u > A):
        return new_iteration
    else:
        return current_iteration


def run_MCMC(obs, initial_params, iterations, filenames, log_directory, factor=0.1):

    logging.info('Runing MCMC for {} iterations. Initial parameters: {}'.format(
        iterations, initial_params))

    update_progress(0)

    best_iteration = run_iteration(initial_params, obs, filenames)
    old_params = initial_params[:]

    chi2_values = []

    for x in range(iterations):

        logging.info(
            '-------------------Begin iteration {} / {}.----------------------'.format(x + 1, iterations))

        new_params = randomize_parameters(old_params, factor=factor)
        logging.info('New randomized parameters: {}'.format(new_params))

        new_iteration = run_iteration(new_params, obs, filenames)

        # Check if there was not an overflow error in the calculation
        if new_iteration[2] != True:

            logging.info('Running Metropolis algorithm...')

            best_iteration = apply_Metropolis(best_iteration, new_iteration)

            logging.info('Ran Metropolis algorithm. Best set of parameters is {} with chi2 error of {}.'.format(
                best_iteration[0], np.sum(np.array(best_iteration[1]))))

            old_params = best_iteration[0]

        else:

            logging.info('Overflow error detected, skipping this value')

        chi2_values.append(np.sum(np.array(best_iteration[1])))

        # if we have found a better set of parameters
        if best_iteration[0] == new_iteration[0]:

            save_chi2(iterations, chi2_values, log_directory)

        update_progress((x+1)/iterations)

    return best_iteration, chi2_values


def plot_chi2(iterations, chi2_vals, filename, display=False):

    plt.plot(np.arange(1, iterations + 1), chi2_vals, color='r')
    plt.grid()
    plt.title('Error vs. Number of Iterations')
    plt.xlabel('Number of Iterations', fontsize=14)
    plt.ylabel(r'Sum of $\chi^2$ Values For All Elements')

    plt.savefig(filename)

    logging.info('Saved error plot.')

    if display:
        plt.show()


def save_chi2(iterations, chi2_vals, path):

    vals = np.array([np.arange(1, iterations + 1), np.array(chi2_vals)])
    np.save(path + 'chi2_vals.npy', vals)
    logging.info('Saved chi 2 values at ()'.format(path))


def update_progress(progress):
    barLength = 10  # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format(
        "#"*block + "-"*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()
# 3

############# LOGGER SETUP ###################


LOG_DIRECTORY = datetime.now().strftime(
    '/home/hmorenom/SSW_Files/OgModel/log/log_%d_%m_%Y_%H_%M_%S/')

os.mkdir(LOG_DIRECTORY)

LOG_FILENAME = LOG_DIRECTORY + 'logfile.log'
IMG_FILENAME = LOG_DIRECTORY + 'chi2_plot.png'

logging.basicConfig(filename=LOG_FILENAME, format='%(asctime)s %(levelname)-8s %(message)s',
                    level=logging.INFO, datefmt='%m-%d-%Y %H:%M:%S')

logging.info('Start program')

# parameters from run 19:17:56 on 3/16/21
c1 = 2.9e-17
c2 = 5.314683
c3 = 2e6
c4 = 0.36
c5 = 0.7425
c6 = 2.2e-15
c7 = 27
c8 = 0.8

iterations = 1000

initial_params = [c1, c2, c3, c4, c5, c6, c7, c8]

filenames = ['pred_c.save', 'pred_fe.save']
#filenames = ['pred_c.save', 'pred_n.save', 'pred_o.save', 'pred_ne.save','pred_mg.save','pred_si.save', 'pred_s.save', 'pred_fe.save']

obs = get_observations()

obs = [obs[0], obs[-1]]  # isolate carbon for testing

best_iteration, chi2_vals = run_MCMC(
    obs, initial_params, iterations, filenames, LOG_DIRECTORY, factor=0.10)

logging.info('Finished iterating. Final parameters: {}. Final Chi^2 value: {}.'.format(
    best_iteration[0], np.sum(np.array(best_iteration[1]))))

save_chi2(iterations, chi2_vals, LOG_DIRECTORY)

plot_chi2(iterations, chi2_vals, IMG_FILENAME)

logging.info('End of program.')
