# -*- coding: utf-8 -*-
"""
PHYS20161 Assignment 3: 79Rb Decay

This programme reads in several csv files with data from an Sr/Rb nuclear decay
experiment in the format [time, activity, activity errors]. It then combines the 
data files and corrects the units. An initial fit is performed using a minimisation
chi squared fit of the expected activity function and outliers from the fit are 
removed. The remaining validated data is refitted and the calculated values for
the decay constant of Sr and Rb are used to create chi squared contour plots, 
from which the uncorrelated errors on the parameters are estimated. Finally, 
the decay constants, half lives, and their respective errors are quoted and 
plotted with all of the data.

Kiran Jonathan 3/12/19
"""

import numpy as np
import matplotlib.pyplot as plt
import tkinter, tkinter.filedialog
from scipy import constants as pc
from scipy.optimize import fmin
from mpl_toolkits.mplot3d import Axes3D
#Warning is given that Axes3D is not used but it must be imported for the 3d projection for the 3d contour plot

"""--------------------{ Data Specifics }--------------------"""

initial_moles_Sr = 1e-6 #(moles) given in script
initial_number_of_nuclei_Sr = initial_moles_Sr * pc.Avogadro #(unitless)
expected_decay_constant_Sr = 0.005 #(per second) given in script
expected_decay_constant_Rb = 0.0005 #(per second) given in script

file_names = ['Nuclear_data_1.csv', 'Nuclear_data_2.csv']
number_of_fitted_parameters = 2
data_lower_bounds = [0,0,0] #all values must be positive
data_upper_bounds = [2147483647,2147483647,2147483647] #upper bounds defined based on maximum integer for 32-bit system
data_units = [pc.hour, pc.tera, pc.tera]


"""--------------------{ Functions }--------------------"""


def file_exists(file_directory):
    """
    Checks if a file exists at a specified directory by trying to open it.
    
    file_directory (string)

    Kiran Jonathan 12/12/19
    """
    
    try:
        
        temp_file = open(file_directory,'r')
            
        temp_file.close()
            
        return True
        
    except:
            
        return False
            
    
def get_file_location(file_directory, desired_file):
    """
    Opens a file explorer and let's the user choose a file; returning the file's
    directory.
    
    file_directory (string)
    desired_file (string)

    Kiran Jonathan 12/12/19
    """
    
    temporary_root = tkinter.Tk()
                
    temporary_root.attributes("-topmost", True)
                        
    temporary_root.withdraw()
                    
    file_directory = tkinter.filedialog.askopenfilename(parent = temporary_root, 
                                                        title = 'Please Select The File: ' 
                                                        + desired_file, filetypes =
                                                        (("csv files","*.csv"),("all files","*.*")))
                    
    temporary_root.destroy()
            
    return file_directory


def is_valid(variable_array, lower_bound_array, upper_bound_array):
    """
    Checks if each element in an array is a float, finite, and between two boundaries
    (specified by corresponding values in the boundary arrays).
    
    variable_array (array)
    lower_bound_array (array)
    upper_bound_array (array)

    Kiran Jonathan 12/12/19
    """
    
    for i in range(len(variable_array)):
        
        try: 
            
            numeric_variable = float(variable_array[i])
            
            if lower_bound_array[i] < numeric_variable < upper_bound_array[i] \
            and np.isfinite(numeric_variable): continue
            
            else: return False
            
        except: return False
    
    return True


def read_data(file_directory, lower_bounds_array, upper_bounds_array):
    """
    Attempts to read in a data file and select only the physical points. If the file
    doesn't exist or doesn't have 3 columns, they will be given 3 attempts to select
    a valid file. Returns a 2d array of the data.
    
    file_directory (string)
    lower_bound_array (array)
    upper_bound_array (array)

    Kiran Jonathan 12/12/19
    """
    
    attempts = 0
    
    max_attempts = 3
    
    desired_file = file_directory
    
    while attempts < max_attempts:
    
        if file_exists(file_directory):
            
            data = np.genfromtxt(file_directory, delimiter = ',')
            
            if len(data[0]) == 3:
                
                data = [data_points for data_points in data if  is_valid(data_points, lower_bounds_array, upper_bounds_array)]
                
                return data, file_directory
            
            else: file_directory = ''
            
        else:
            
            attempts += 1
            
            print('\nPlease select a valid csv file with three columns of data for: ' + desired_file)
            
            file_directory = get_file_location(file_directory, desired_file)
            
    return [] , ''


def get_compiled_data(file_directories, units_array, lower_bounds_array, upper_bounds_array):
    """
    Reads in data files at specified directories and combines them into a single,
    sorted (by column 0) 2d array with corrected units.
    
    file_directory (string)
    units_array (array)
    lower_bound_array (array)
    upper_bound_array (array)

    Kiran Jonathan 12/12/19
    """
    
    data = np.empty((0,3))

    files_added = []
    
    for file_name in file_names:
        
        temp_data, file_name = read_data(file_name, lower_bounds_array, upper_bounds_array)
        
        if len(temp_data) != 0 and not file_name in files_added:
            
            data = np.vstack((data, temp_data))
            
            files_added.append(file_name)
            
        else:
            
            return []
        
    return data * units_array
    


def validate_data(data, fit_array):
    """
    Removes outliers in a data set which are more than three standard deviations
    away from the given fit.
    
    data (array)
    fit_array (array)

    Kiran Jonathan 12/12/19
    """
    
    indexes_to_delete = []
    
    for i in range(len(data)):
        
        if np.absolute(fit_array[i] - data[i,1]) / data[i,2] > 3:
            
            indexes_to_delete.append(i)
    
    deleted = 0
    
    indexes_to_delete = list(dict.fromkeys(indexes_to_delete))
    
    for index in indexes_to_delete:
        
        data = np.delete(data, index - deleted, 0)
        
        deleted +=1   
    
    return data


def expected_fit(time_values, decay_constant_Sr, decay_constant_Rb):
    """
    Returns the expected activity at the given times given the initial number 
    of nuclei and the two decay constants.
    
    time_values (array)
    decay_constant_Sr (float)
    decay_constant_Rb (float)

    Kiran Jonathan 12/12/19
    """
    
    return initial_number_of_nuclei_Sr * (decay_constant_Rb * decay_constant_Sr) / \
            (decay_constant_Rb - decay_constant_Sr) * (np.exp(-decay_constant_Sr \
            * time_values) - np.exp(-decay_constant_Rb * time_values))


def chi_squared_of_fit(decay_constants, time_array, activity_array, activity_error_array):
    """
    Returns the chi squared value for the fit associated with the given decay
    constants for the data given.
    
    decay_constants (array)
    time_array (array)
    activity_array (array)
    activity_errors_array (array)
    
    Kiran Jonathan 12/12/19
    """
    
    chi_squared = 0
    
    for i in range(len(time_array)):
        
        difference = (expected_fit(time_array[i], decay_constants[0], decay_constants[1]) - activity_array[i])
        
        chi_squared += np.square(difference/activity_error_array[i])
    
    return chi_squared


def calculate_decay_constants(decay_constant_Sr_start, decay_constant_Rb_start, args_array):
    """
    Uses the fmin function from scipy with the given starting points and required
    arguments for the chi squared function to find the decay constants which give the
    fit with the lowest associated chi squared.
    
    decay_constant_Sr_start (float)
    decay_constant_Rb_start (float)
    arg_array (array)

    Kiran Jonathan 12/12/19
    """
    
    return fmin(chi_squared_of_fit, (decay_constant_Sr_start, decay_constant_Rb_start), 
                args = args_array , full_output = True, disp = False)
            

def calculate_half_life(decay_constant):
    """
    Returns the half life associated with the given decay constant.
    
    decay_constant (float)

    Kiran Jonathan 12/12/19
    """
    
    return np.log(2) / decay_constant


def find_decay_constant_errors(decay_constants_array):
    """
    Returns the uncorrelated error for two decay constants given an array of pairs
    of decay constants which, when fitted, give a chi squared value one larger
    than the minimum chi squared. The uncorrelated error calculated by this function
    is equal to half of the difference between the largest decay constant value
    and the smallest decay constant value which lie on the minimum chi squared plus
    one ellipse/contour in the chi squared contour map. i.e. for the decay constant
    which lies along the x axis in the contour plot it is equal to half of the 
    difference between the largest x coordinate and smallest x coordinate of points
    lying on the minimum chi squared plus one contour.
    
    decay_constants_array (array)

    Kiran Jonathan 12/12/19
    """
    
    decay_constant_Sr_uncertainty = (np.max(decay_constants_array[:,0]) - np.min(decay_constants_array[:,0])) / 2
    
    decay_constant_Rb_uncertainty = (np.max(decay_constants_array[:,1]) - np.min(decay_constants_array[:,1])) / 2
    
    return decay_constant_Sr_uncertainty, decay_constant_Rb_uncertainty


def calculate_half_life_errors(half_life, decay_constant, decay_constant_uncertainty):
    """
    Propagates the error on the decay constant into an error for the half life
    by using the fact that they should have the same fractional uncertainty
    (quotient relationship with constant multiplier).
    
    half_life (float)
    decay_constant (float)
    decay_constant_uncertainty (float)

    Kiran Jonathan 12/12/19
    """
    
    return half_life * decay_constant_uncertainty / decay_constant


def formatted_value_and_uncertainty(value, uncertainty, significant_figures, units):
    """
    Returns a string of the value rounded to the specified number of significant 
    figures 'plus or minus' the uncertainty rounded to the same number of decimal
    places as the value, along with the units given.
    
    value (float)
    uncertainty (float)
    significant_figures (int)
    units (string)

    Kiran Jonathan 12/12/19
    """
    
    formatted_value = ('{0:.' + str(significant_figures) + 'g}').format(value)
    
    while (formatted_value[-3] == '0' or formatted_value[-3] == '.') \
            and (formatted_value[-4] == '0' or formatted_value[-4] == '.'): 
        
        formatted_value += '0'

    """While loop ensures the correct number of significant figures is shown for
    the value (and therefore uncertainty) as the .g formatter will not display
    the zeros at the end of any formatted float"""

    consistent_decimal_places = len(formatted_value) - 1 - formatted_value.index('.')
    
    formatted_uncertainty = ('{0:.' + str(consistent_decimal_places) + 'f}').format(uncertainty)
    
    return '(' + formatted_value + u' \u00B1 ' + formatted_uncertainty + ')' + units


def formatted_standard_deviations_away(value, value_error, expected_value):
    """
    Returns a formatted string of how many standard deviations a given value
    is from the expected value.
    
    value (float)
    value_error (float)
    expected_value (float)

    Kiran Jonathan 12/12/19
    """
    
    return 'That is {0:.2f} standard deviations from the accepted value of {1:.3g}.'\
          .format( (np.abs(value - expected_value) / value_error) , expected_value)


def plot_with_errors(x_array, y_array, y_errors, x_errors = 0, fit_x_linspace = [], 
                     fit_y = [], fit_label = "Calculated Fit", plot_title = "Plot", 
                     x_label = "x", y_label = "y", save_name = '', residuals_array = [],
                     annotation = ''):
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
    residuals_array (array)
    annotation (string)

    Kiran Jonathan 19/11/19
    """
    
    figure = plt.figure(plot_title, figsize=(10,9))
    
    figure_grid = figure.add_gridspec(ncols = 1, nrows = 3, height_ratios = [8,1,2])
    
    error_plot = figure.add_subplot(figure_grid[0,0])
    
    error_plot.set_title(plot_title, fontsize = 16)
        
    error_plot.set_xlabel(x_label, fontsize = 14)
    
    error_plot.set_ylabel(y_label, fontsize = 14)
    
    error_plot.grid(alpha = 0.4, linestyle = '--')
    
    error_plot.errorbar(x_array,y_array, xerr = x_errors, yerr = y_errors, 
                 capsize = 3, fmt='D', c='black', ecolor = '#d24d4d', marker = 'x')
    
    if len(fit_y) != 0 and len(fit_x_linspace) != 0: 
        
        error_plot.plot(fit_x_linspace, fit_y, linestyle = '-', label = fit_label, 
                 alpha = 0.7, c = 'cornflowerblue')
    
        error_plot.legend(fontsize = 12)
        
    if annotation != '': error_plot.annotate(annotation, (0,-0.09), xytext = (0,0), 
                                             va = 'top', xycoords = 'axes fraction', textcoords = 'offset points')
        
    if len(residuals_array) != 0:
        
        residuals = figure.add_subplot(figure_grid[2,0])    
        
        residuals.set_title('Residuals', fontsize = 16)
            
        residuals.set_xlabel(x_label, fontsize = 12)
        
        residuals.set_ylabel('Residual\n' + y_label, fontsize = 12)
        
        residuals.errorbar(x_array, residuals_array, yerr = y_errors,
                           capsize = 3, fmt='D', c='black', ecolor = '#d24d4d', marker = 'x')
                           
        residuals.plot(x_array, np.zeros((len(x_array),)), alpha = 0.7)
    
    if save_name != '': plt.savefig(save_name, dpi = 300)
    
    plt.show()
    
    pass


def chi_squared_parameters_contour_plot_meshes(time_array, activity_array, activity_errors_array, 
                                               chi_squared_fit_results, fractional_range = 0.1):
    """
    Generates the meshes for the two parameters being used for the x and y coordinates
    of the chi squared contour plot (and 3d chi squared surface plot), as well
    as the corresponding chi squared mesh from the expected fit corresponding to
    the two parameters.
    
    time_array (array)
    activity_array (array)
    activity_errors_array (array)
    chi_squared_fit_results (array)
    fractional_range (float)

    Kiran Jonathan 12/12/19
    """
    
    decay_constant_Sr_local, decay_constant_Rb_local = chi_squared_fit_results[0]
    
    decay_constant_Sr_values = np.linspace(decay_constant_Sr_local * (1-fractional_range), 
                                           decay_constant_Sr_local * (1+fractional_range), 500)
    
    decay_constant_Rb_values = np.linspace(decay_constant_Rb_local * (1-fractional_range), 
                                           decay_constant_Rb_local * (1+fractional_range), 500)
        
    decay_constant_Sr_values_mesh , decay_constant_Rb_values_mesh = np.meshgrid(decay_constant_Sr_values, decay_constant_Rb_values)
        
    chi_squared_mesh = chi_squared_of_fit((decay_constant_Sr_values_mesh, decay_constant_Rb_values_mesh), 
                                           time_array, activity_array, activity_errors_array)
    
    return decay_constant_Sr_values_mesh, decay_constant_Rb_values_mesh, chi_squared_mesh    
    

def chi_squared_parameters_contour_plot(x_mesh, y_mesh, height_mesh, minimum_chi_squared, 
                             minimum_chi_sqaured_xy, title = "Chi Squared Contour Plot",
                             plot_title = "Chi Squared Contour Plot",
                             x_label = "x", y_label = "y", save_name = ''):
    """
    Produces a contour plot with the chi squared value associated with the fit
    from the x and y parameters as the contour height. It plots contour lines
    correcponding to the minimum chi squared values plus the increments
    corresponding to n standard deviations. It returns the plot of the minimum
    chi squared plus one contour.
    
    x_mesh (array)
    y_mesh (array)
    height_mesh (array)
    minimum_chi_squared (float)
    minimum_chi_sqaured_xy (array)
    title (string)
    plot_title (string)
    x_label (string)
    y_label (string)
    save_name (string)

    Kiran Jonathan 12/12/19
    """
    
    contour_figure = plt.figure(plot_title, figsize=(10,8))
        
    parameters_contour_plot = contour_figure.add_subplot(111)
        
    x_min_max = [x_mesh[0,0],x_mesh[0,-1]]
    y_min_max = [y_mesh[0,0],y_mesh[-1,0]]
    
    parameters_contour_plot.set_xlim(x_min_max[0], x_min_max[1])
    parameters_contour_plot.set_ylim(y_min_max[0], y_min_max[1])
    
    parameters_contour_plot.set_xticks(np.linspace(x_min_max[0], x_min_max[1], 5))
    parameters_contour_plot.set_yticks(np.linspace(y_min_max[0], y_min_max[1], 5))
    
    parameters_contour_plot.set_title(title, fontsize = 16)
    parameters_contour_plot.set_xlabel(x_label, fontsize = 15)
    parameters_contour_plot.set_ylabel(y_label, fontsize = 15)
    
    parameters_contour_plot.contourf(x_mesh, y_mesh, height_mesh, cmap = 'YlGnBu', levels = 10)
    
    chi_plus_one_contour = parameters_contour_plot.contour(x_mesh,y_mesh, height_mesh, 
                                                           levels = [minimum_chi_squared + 1.00], 
                                                           linestyles = 'dashed', colors = 'darkblue')
    
    chi_differences = np.array([2.3,5.99,9.21])
    
    chi_levels = chi_differences + minimum_chi_squared
    
    contour_plot = parameters_contour_plot.contour(x_mesh, y_mesh, height_mesh, levels = chi_levels, 
                                                   linestyles = '-', cmap = 'tab20b')
    
    parameters_contour_plot.clabel(contour_plot, fmt = "%.2f")

    minimum_point = parameters_contour_plot.scatter(minimum_chi_sqaured_xy[0], minimum_chi_sqaured_xy[1], 
                                                    marker = 'x', color = 'darkgoldenrod')

    labels = [r'$\chi^2_{{\mathrm{{min.}}}}+2.30$',r'$\chi^2_{{\mathrm{{min.}}}}+5.99$',
                                                   r'$\chi^2_{{\mathrm{{min.}}}}+9.21$']

    box = parameters_contour_plot.get_position()
    
    parameters_contour_plot.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    minimum_point.set_label(r'$\chi^2_{{\mathrm{{min.}}}}$ = '+'{0:.2f}'.format(minimum_chi_squared))
    
    chi_plus_one_contour.collections[0].set_label(r'$\chi^2_{{\mathrm{{min.}}}}+1.00$')

    for i in range(len(labels)):
        
        contour_plot.collections[i].set_label(labels[i])
        
    parameters_contour_plot.legend(loc='center left', bbox_to_anchor=(1, 0.5), 
                                   fontsize = 14)

    if len(save_name) != 0: plt.savefig(save_name, dpi = 300)
    
    plt.show()
    
    return chi_plus_one_contour


def chi_squared_parameters_3d_plot(x_mesh, y_mesh, height_mesh, plot_title = "Chi Squared Surface Plot",
                                   title = "Chi Squared 3d Plot",
                                   x_label = "x", y_label = "y", z_label = "z", save_name = ''):
    """
    Produces a contour surface plot with the chi squared value associated with the fit
    from the x and y parameters as the contour height for demonstative purposes to
    visualize the shape of the contour plot.
    
    x_mesh (array)
    y_mesh (array)
    height_mesh (array)
    plot_title (string)
    title (string)
    x_label (string)
    y_label (string)
    z_label (string)
    save_name (string)

    Kiran Jonathan 12/12/19
    """
    
    surface_figure = plt.figure(plot_title, figsize=(10,8))
        
    parameters_3d_plot = surface_figure.add_subplot(111, projection = '3d')
    
    parameters_3d_plot.set_title(title, fontsize = 16)
    parameters_3d_plot.set_xlabel(x_label, fontsize = 15)
    parameters_3d_plot.set_ylabel(y_label, fontsize = 15)
    parameters_3d_plot.set_zlabel(z_label, fontsize = 15)

    parameters_3d_plot.contour3D(meshes[0], meshes[1], meshes[2], cmap = 'plasma', levels = 50)

    if len(save_name) != 0: plt.savefig(save_name, dpi = 300)
    
    plt.show()
    
    pass


"""--------------------{ Main Body }--------------------"""
         

raw_data = get_compiled_data(file_names, data_units, data_lower_bounds, data_upper_bounds)
    
if len(raw_data) != 0: #Checks if a valid file was selected within the maximum number of attempts
    
    """--------------------{ Data Validation }--------------------"""
    
    fit_results = calculate_decay_constants(expected_decay_constant_Sr, expected_decay_constant_Rb, 
                                            (raw_data[:,0], raw_data[:,1], raw_data[:,2]))
    
    validated_data = validate_data(raw_data, expected_fit(raw_data[:,0], fit_results[0][0], fit_results[0][1]))
    
    time_data = validated_data[:,0]
    activity_data = validated_data[:,1]
    activity_errors = validated_data[:,2]
    
    fit_results = calculate_decay_constants(fit_results[0][0], fit_results[0][1], (time_data, activity_data, activity_errors))
    
    decay_constant_Sr, decay_constant_Rb = fit_results[0]
        
    """--------------------{ Graph Plotting }--------------------"""
    
    print("\nGenerating contour plots...")

    meshes = chi_squared_parameters_contour_plot_meshes(time_data, activity_data, activity_errors, fit_results)

    chi_squared_plot = chi_squared_parameters_contour_plot(meshes[0], meshes[1], meshes[2], fit_results[1], fit_results[0],
                                                           title = "$\chi^2$ Contour Plot for Values of $^{Sr}\lambda$ and $^{Rb}\lambda$",
                                                           x_label = "$^{Sr}\lambda$ (s$^{-1}$)", y_label = "$^{Rb}\lambda$ (s$^{-1}$)",
                                                           save_name = 'chi_squared_contour_plot.png')
    
    chi_squared_parameters_3d_plot(meshes[0], meshes[1], meshes[2], title = "$\chi^2$ 3D Contour Plot", x_label = "$^{Sr}\lambda$",
                                   y_label = "$^{Rb}\lambda$", z_label = "$\chi^2$")

    print("\nContour plots generated.")

    """--------------------{ Value and Error Calculation }--------------------"""

    reduced_chi_squared = fit_results[1] / (len(time_data) - number_of_fitted_parameters)
    
    half_life_Sr = calculate_half_life(decay_constant_Sr) / pc.minute
    half_life_Rb = calculate_half_life(decay_constant_Rb) / pc.minute

    decay_constant_Sr_error, decay_constant_Rb_error = find_decay_constant_errors(chi_squared_plot.allsegs[0][0])
    
    half_life_Sr_error = calculate_half_life_errors(half_life_Sr, decay_constant_Sr, decay_constant_Sr_error)
    half_life_Rb_error = calculate_half_life_errors(half_life_Rb, decay_constant_Rb, decay_constant_Rb_error)
    
    """--------------------{ Printing Results }--------------------"""
    
    formatted_decay_constant_Sr = formatted_value_and_uncertainty(decay_constant_Sr, decay_constant_Sr_error, 3, ' /s')
    formatted_decay_constant_Rb = formatted_value_and_uncertainty(decay_constant_Rb, decay_constant_Rb_error, 3, ' /s')
    formatted_half_life_Sr = formatted_value_and_uncertainty(half_life_Sr, half_life_Sr_error, 3, ' Minutes')
    formatted_half_life_Rb = formatted_value_and_uncertainty(half_life_Rb, half_life_Rb_error, 3, ' Minutes')
    
    print('\nSr Decay Constant = ' + formatted_decay_constant_Sr)
    print(formatted_standard_deviations_away(decay_constant_Sr, decay_constant_Sr_error, expected_decay_constant_Sr))
    print('Sr Half Life = ' + formatted_half_life_Sr)
    print('\nRb Decay Constant = ' + formatted_decay_constant_Rb)
    print(formatted_standard_deviations_away(decay_constant_Rb, decay_constant_Rb_error, expected_decay_constant_Rb))
    print('Rb Half Life = ' + formatted_half_life_Rb)
    print('\nReduced \u03C7 Squared = {0:.2f}'.format(reduced_chi_squared))
    
    """--------------------{ Final Plot }--------------------"""
    
    print("\nGenerating final plot...")
    
    fit_linspace = np.linspace(min(time_data), max(time_data), 500)
    
    plot_with_errors(time_data, activity_data, activity_errors, fit_x_linspace= fit_linspace,
                     fit_y= expected_fit(fit_linspace, decay_constant_Sr, decay_constant_Rb),
                     plot_title = 'Activity Vs. Time', x_label = 'Time (s)', y_label = 'Activty (Bq/s)', 
                     fit_label= "Minimised $\chi^2$ fit"+ "\nReduced $\chi^2$ = {0:.2f}".format(reduced_chi_squared),
                     save_name = 'activity_against_time.png', 
                     residuals_array = activity_data - expected_fit(time_data, decay_constant_Sr, decay_constant_Rb),
                     annotation = '\n$^{Sr}\lambda = $' + formatted_decay_constant_Sr
                     + '\t\t\t$^{Sr}t_{1/2} = $' + formatted_half_life_Sr
                     + '\n$^{Rb}\lambda = $' + formatted_decay_constant_Rb
                     + '\t\t$^{Rb}t_{1/2} = $' + formatted_half_life_Rb)
    
    print("\nFinal plot generated.")
    
else:
    
    print('\nMaximum attempts to open data file exceeded. Please check you are adding different data files with three columns and try again.')

print('\nThank you for using this programme.')

