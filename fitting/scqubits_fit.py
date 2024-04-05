import scqubits as scq
from scipy.optimize import curve_fit
import numpy as np

'''
The idea here is to write something that will take in arrays of array for the freq/ flux 
values from a qubit scan and fit the transmon curve to it:
- we expect the g-e transition stuff to be the first set of arrays
- we might get a few points from the g-f/2 transition: here we probably need to make sure
we multiply the frequency values by 2 in the original code so we can fit that as well 
- probably long term we should write something that takes a list of the transitions we're
providing (like e,f,etc.) and grabs the right ones from scqubits 

- current problem: the fitting function really doesn't like ragged arrays so we have to do some massaging around
'''

def qubit_energies(param_dict, qubit_type, flux_vals):
    if qubit_type=='Transmon':
        qubit_params = np.array([param_dict['EJ_sum'], param_dict['EJ_diff'], param_dict['EC']])
        qubit_function = transmon
    return qubit_function(flux_vals, *qubit_params)    

def fit_function(param_dict, qubit_type, flux_vals, freq_vals):
    if qubit_type=='Transmon':
         params = ['EJ_sum', 'EJ_diff', 'EC']
         init_vals = np.array([param_dict['EJ_sum'], param_dict['EJ_diff'], param_dict['EC']])
         qubit_function = transmon
    
    bounds = (0.5*init_vals, 1.5*init_vals)
    popt, pcov = curve_fit(qubit_function, flux_vals, freq_vals, method='trf', 
                            p0 = init_vals, bounds = bounds, maxfev = 10000)
    fit_params = {}
    for ind, p in enumerate(params):
         fit_params[p] = popt[ind]
    return fit_params

def transmon(phi, EJ_sum, EJ_diff, EC):
    d = EJ_diff/EJ_sum
    transmon = scq.TunableTransmon(EJmax=EJ_sum, EC=EC, d=d, flux=0, ng=0.5, ncut=30)
    spectrum = transmon.get_spectrum_vs_paramvals(param_name='flux', param_vals=phi, evals_count=2, subtract_ground=True)
    energies = spectrum.energy_table.transpose()[1]
    return energies   
 
def transmon_levels(phi, EJ_sum, EJ_diff, EC, num_levels):
    d = EJ_diff/EJ_sum
    energies = []
    transition_energy = []
    transmon = scq.TunableTransmon(EJmax=EJ_sum, EC=EC, d=d, flux=0, ng=0.5, ncut=30)
    spectrum = transmon.get_spectrum_vs_paramvals(param_name='flux', param_vals=phi, evals_count=num_levels+2, subtract_ground=True)
    energies = spectrum.energy_table.transpose()[1:num_levels]
    return energies

def transmon_fitting_wrapper(flux_vals, freq_vals, EJ_sum, EJ_diff, EC, num_points):
    flux_points = np.linspace(-1,1, 501)
    interp_energies = np.zeros(np.sum(num_points))
    offset = 0
    energies = transmon_levels(flux_points, EJ_sum, EJ_diff, EC, len(num_points)+1)

    for flux_list, energy in zip(phi, energies):
        evaluated_vals = np.interp(flux_list, flux_vals, energy)
        interp_energies[0+offset:len(evaluated_vals)+offset] = evaluated_vals
        offset+= len(evaluated_vals)
    return interp_energies
def bare_transitions_fit(flux_vals, freq_vals, guess):
        
        flattened_flux = np.concatenate(flux_vals)
        flattened_freqs = np.concatenate(freq_vals)
        num_points = np.array([len(a) for a in flux_vals])

        bounds = (0.5*guess, 1.5*guess)
        popt, pcov = curve_fit(transmon_levels, flux_vals, freq_vals, method='trf', 
                                p0 = guess, bounds = bounds, maxfev = 10000)
        return {'Ec':popt[2],'EJ_sum':popt[0],'EJ_diff':popt[1]}
        