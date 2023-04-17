# -*- coding: utf-8 -*-

"""
PHYS20161 Assignment 2: Thickness of BN

This programme reads in a data file of electron energy values, transmission coefficients,
and their errors (energy errors are optional) from an experiment on quantum
tunnelling in BN. It then proceeds to calculate the thickness of the BN, d, by
using a hill climbing algorithm with a varying step size to minimise the reduced
chi squared with respect to d. It also estimates the error on d using the
statistcal property of a difference in 1 in the chi squared corresponding to
a 68% confidence level (1 standard deviation). The approximate number of 
layers is also calculated. 

Kiran Jonathan 19/11/19
"""

import numpy as np
import matplotlib.pyplot as plt
import tkinter, tkinter.filedialog

"""Physical Constants (all unites relative to electron-volts and Angstrom)"""

V_0 = 3. # (eV)
relative_permittivity = 4 # (dimensionless)
permittivity_of_free_space = 0.00553 # (eV * Angstrom ^-1)
mass_planck_constant = 0.512317 # (eV ^-1/2 * Angstrom ^-1) Equal to sqrt(2m) / h_bar where m is the electron mass
electron_charge = 1 # (eV)
lambda_times_d =  ( np.square(electron_charge) * np.log(2) ) \
                    / ( 8 * np.pi * relative_permittivity * permittivity_of_free_space ) # (eV * Angstrom)

"""Data Specifics"""

approximate_single_layer_thickness = 3 
file_name = 'Tunnelling_data_BN.csv'

physical_data_lower_bounds = [0,0,0,0] #Physical lower bounds fo T, E, T error, E error
physical_data_upper_bounds = [1,V_0,1,V_0] #Physical upper bounds fo T, E, T error, E error
"""Transmission coefficient boundaries defined 0 to 1 as it is a fraction.
Electron energy boundaries defined to ensure physicality of the fit (energy values
greater than V_0 will produce imaginary numbers for expected_fit() due the root term."""
 

def file_exists(file_directory):
    """
    Checks if a file exists at the specified directory by trying to open it.
    
    file_directory (string)

    Kiran Jonathan 19/11/19
    """
    
    try:
        temp_file = open(file_directory,'r')
        
        temp_file.close()
        
        return True
    
    except:
        
        return False

def get_file_location():
    """
    Asks the user to choose a file they would like to read in as the data from
    a file explorer.

    Kiran Jonathan 19/11/19
    """
    
    temporary_root = tkinter.Tk()

    temporary_root.attributes("-topmost", True)
        
    temporary_root.withdraw()
    
    return tkinter.filedialog.askopenfilename(parent = temporary_root, title = 
                                              'Please Select File', filetypes = 
                                              (("csv files","*.csv"),("all files","*.*")))
 

def is_valid(variable_array, lower_bound_array, upper_bound_array):
    """
    Checks if each element in an array is a float and between two boundaries
    (specified by corresponding values in the boundary arrays).
    
    variable_array (array)
    lower_bound_array (array)
    upper_bound_array (array)

    Kiran Jonathan 19/11/19
    """
    
    for i in range(len(variable_array)):
        
        try:       
            numeric_variable = float(variable_array[i])
            
            if lower_bound_array[i] <= numeric_variable < upper_bound_array[i]: continue
            
            else: return False
            
        except: return False
    
    return True
    

def read_data(file_directory, data_lower_bounds, data_upper_bounds):
    """
    Attempts to open and read a file, appending all validated data points to a 
    multidimensional numpy array which it returns. If the specified file does not
    exist or is not valid it will use get_file_location() to obtain a directory 
    of a valid file.
    
    file_directory (string)
    data_lower_bounds (array)
    data_upper_bounds (array)

    Kiran Jonathan 19/11/19
    """
    
    maximum_attempts = 3 #Maximum attempts to open a valid file
    
    attempts = 0 #Current number of attempts
    
    while attempts < maximum_attempts:
        
        if file_exists(file_directory):
    
            file = open(file_directory,'r')
                
            invalid_counter = 0
            
            first_line = True
            
            for line in file:
            
                split_line = line.split(',')
                
                if first_line: 
                    
                    data = np.empty((0,len(split_line)))
                    
                    first_line = False
                
                if split_line[0][0] == '%': continue
                
                elif is_valid(split_line, data_lower_bounds, data_upper_bounds):
                    
                    data = np.vstack((data,list(map(float,split_line))))
                
                else: invalid_counter += 1
                
            file.close()
            
            if len(split_line) >= 3:
            
                print('\n{0} data points were determined to be invalid '.format(invalid_counter) +
                  'and were excluded, this does not include comments (denoted by a % sign at the start of the line).\n')
            
                return data 
            
            else:
                
                file_directory = ''
        
        else:
            
            attempts += 1
            
            try:
        
                print("\nFile not found/invalid, please select the correct data file" +
                      "and make sure it has at least three columns of data.")

                file_directory = get_file_location()
            
            except:
            
                print("\nNo file could be opened. Please move" + file_name 
                  + "to the same directory as this script and run it again.")
                
                return []

    return []


def expected_fit(electron_energy_array, d_fit):
    """
    Calculates an array of transmission coefficients corresponding to the data's
    electron energy values for a given value of d.    
    
    electron_energy_array (array)
    d_fit (float)

    Kiran Jonathan 19/11/19
    """
    
    d_1 = (1.2 * lambda_times_d) / V_0
    
    d_2 = d_fit - d_1
    
    average_potential = V_0 - ( ( (1.15 * lambda_times_d) / (d_2-d_1) ) * np.log( ( (d_fit - d_1) * d_2 ) / ( (d_fit - d_2) * d_1) ) )
    
    return np.exp( ((-2) * (d_2-d_1) * mass_planck_constant) * np.sqrt(average_potential - electron_energy_array) )


def chi_squared(observed_data, observed_data_error, expected_data):
    """
    Calculates the chi squared value for a trial fit from the observed data.
    
    observed_data (array)
    observed_data_error (array)
    expected_data (array)

    Kiran Jonathan 19/11/19
    """
    
    return np.sum( np.square( (observed_data - expected_data) / observed_data_error ) )


def estimate_d(x_data, y_data, y_error, single_layer_thickness, maximum_layers):
    """
    Finds an estimate for d by calculating the chi squared for multiples of the
    approximate single layer thickness and returning the value corresponding to
    the lowest chi squared.
    
    x_data (array)
    y_data (array)
    y_error (array)
    single_layer_thickness (float)
    maximum_layers (int)

    Kiran Jonathan 19/11/19
    """
    
    d_test_values = single_layer_thickness * np.arange(maximum_layers+1)
    
    d_chi_values = np.array([])
    
    for d_test in d_test_values:
        
        d_chi_values = np.append(d_chi_values, chi_squared(y_data, y_error, expected_fit(x_data, d_test)))
    
    return d_test_values[np.argmin(d_chi_values[np.isfinite(d_chi_values)])] 
    #Returns the value of d corresponding to the lowest non-nan value of chi squared


def hill_climb_chi_squared(x_data, y_data, y_errors, d_start, chi_squared_target, tolerance = 0.0001):
    """
    Uses a hill climbing algorithm with a varying step size to minimise the
    the difference between the test chi squared and a target chi squared with
    respect to the free parameter d. Tolerance implicitly defined as 0.0001 as 
    a difference of two fitted parameters must be known to 3 decimal places.
    
    x_data (array)
    y_data (array)
    y_errors (array)
    d_start (float)
    chi_squared_target (float)
    tolerance (float)
    
    Kiran Jonathan 19/11/19
    """
    
    d_step = 0.2
    
    d_current = d_start
    
    while True:
        
        chi_d_difference = np.abs(chi_squared(y_data, y_errors, expected_fit(x_data, d_current)) - chi_squared_target)
        
        chi_d_difference_plus = np.abs(chi_squared(y_data, y_errors, expected_fit(x_data, d_current + d_step)) - chi_squared_target)
        
        chi_d_difference_minus = np.abs(chi_squared(y_data, y_errors, expected_fit(x_data, d_current  -d_step)) - chi_squared_target)
        
        if chi_d_difference_minus < chi_d_difference: d_current -= d_step
            
        elif chi_d_difference_plus < chi_d_difference: d_current += d_step

        elif d_step >= tolerance: d_step /= 2
            
        else: 
            
            return np.round(d_current, 4) , chi_squared(y_data, y_errors, expected_fit(x_data, round(d_current,4)))
            #d_current is rounded to the precision to which the tolerance ensures it is known.           
 
    
def parameter_error(x_data, y_data, y_errors, parameter, parameter_chi_squared):
    """
    Calculates the error on a parameter calculated from a chi squared minimisation 
    fit by exploiting the statistical property that one above the chi squared
    corresponds to a 68% confidence level, i.e. to one standard deviation.
    
    x_data (array)
    y_data (array)
    y_errors (array)
    fit_parameter (float)
    fit_parameter_chi_squared (float)
    
    Kiran Jonathan 19/11/19
    """
    
    parameter_minus_sigma = hill_climb_chi_squared(x_data,  y_data, y_errors, 0.9 * parameter,
                                                 parameter_chi_squared + 1)[0]
    #parameter is multiplied by 0.9 to apply a bias in the left direction
    

    parameter_plus_sigma = hill_climb_chi_squared(x_data, y_data, y_errors, parameter + 
                                      np.abs(parameter - parameter_minus_sigma), parameter_chi_squared + 1)[0]
    #starting point is defined based on the expected position of d_plus_sigma if the chi squared distribution is symmetric


    return np.abs(parameter_plus_sigma - parameter_minus_sigma)/2
    #returns the average value of sigma in case the chi squared distribution is antisymmetric

def plot_with_errors(x_array, y_array, y_errors, x_errors = 0, fit_x_linspace = [], 
                     fit_y = [], fit_label = "Calculated Fit", plot_title = "Plot", 
                     x_label = "x", y_label = "y", save_name = ''):
    """
    Plots a customised y error or x-y error graph with the given data and labels.
    
    x_array (array)
    y_array (array)
    y_errors (array)
    x_errors (array)
    fit (array)
    fit_label  (string)
    plot_title (string)
    x_label (string)
    y_label (string)
    save_name (string)

    Kiran Jonathan 19/11/19
    """
    
    plt.figure(plot_title, figsize=(8,5))
    
    plt.title(plot_title, fontsize = 18)
        
    plt.xlabel(x_label, fontsize = 14)
    
    plt.ylabel(y_label, fontsize = 14)
    
    plt.grid(alpha = 0.4, linestyle = '--')
    
    plt.errorbar(x_array,y_array, xerr = x_errors, yerr = y_errors, 
                 capsize = 3, fmt='D', c='black', ecolor = '#d24d4d', marker = 'x')
    
    if len(fit_y) != 0 and len(fit_x_linspace) != 0: 
        
        plt.plot(fit_x_linspace, fit_y, linestyle = '-', label = fit_label, 
                 alpha = 0.7, c = 'cornflowerblue')
    
        plt.legend(fontsize = 12)
    
    if save_name != '': plt.savefig(save_name, dpi = 300)
    
    plt.show()
    
    pass



"""Main Body:"""
         
raw_data = read_data(file_name,physical_data_lower_bounds,physical_data_upper_bounds)

if len(raw_data) != 0: #Checks if a valid file was selected within the maximum number of attempts

    transmission_data = raw_data[:,0]
    electron_energy_data = raw_data[:,1]
    electron_energy_data_errors = 0
    transmission_errors = raw_data[:,2]
    
    if len(raw_data[0,:]) == 4: electron_energy_data_errors = raw_data[:,3]
    """If a column of electron energy errors is present it will be read in and plotted,
    though it has no effect on the calculated values."""
            
    d_estimate = estimate_d(electron_energy_data, transmission_data, transmission_errors, approximate_single_layer_thickness, 100)
    
    d, chi_squared_fit = hill_climb_chi_squared(electron_energy_data, transmission_data, transmission_errors, d_estimate, 0)
    #target_chi_squared is zero as this finds the d value for which chi squared is minmised.
    
    reduced_chi_squared = chi_squared_fit / (len(transmission_data)-1)
    
    d_error = parameter_error(electron_energy_data, transmission_data, transmission_errors, d, chi_squared_fit)
    
    print(u'Barrier thickness, d = ({0:.3f} \u00B1 {1:.3f})\u212B.'.format(d,d_error))
    print('Reduced \u03C7 Squared = {0:.2f}.'.format(reduced_chi_squared))
    print('There are approximately {0:.0f} layer(s)'.format(d / approximate_single_layer_thickness))
    print('\nA graph is now being generated...')  
        
    electron_energy_linspace = np.linspace(min(electron_energy_data), max(electron_energy_data), 500)
    
    plot_with_errors(electron_energy_data, transmission_data, transmission_errors, x_errors = 
                     electron_energy_data_errors, fit_x_linspace = electron_energy_linspace,
                     fit_y = expected_fit(electron_energy_linspace,d), 
                     fit_label = u'Optimised fit: \nd = ({0:.3f} $\pm$ {1:.3f}) \u212B '.format(d,d_error) +
                     '\n$\chi^2_{reduced}$' + ' = {0:.2f}'.format(reduced_chi_squared), 
                     plot_title = 'Transmission Coefficient Vs. Electron Energy', 
                     x_label = 'Electron Energy (eV)', y_label = 'Transmission Coefficient',
                     save_name = 'T_Against_E.png'.format(d,reduced_chi_squared))

else:
    
    print('\nMaximum attempts to open data file exceeded. Please check your data file and run the programme again.')

print('\nThank you for using this programme.')



