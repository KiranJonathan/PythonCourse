# -*- coding: utf-8 -*-
"""
PHYS20161 Assignment 1: Bouncy Ball

This programme asks the user for an initial height at which a ball is dropped, 
a minimum height of interest over which the ball must bounce, and a coefficient 
of efficiency (eta) which relates to the ball’s conservation of gravitational 
potential energy. The user is asked to decide whether bounces with peak heights equal 
to the minimum height of interest should be included. It then proceeds to calculate 
the number of bounces above the minimum height and the time taken for these bounces 
to occur. This time is the time from the dropping of the ball at the initial height 
to when the ball touches the ground after the last bounce above or equal to the 
minimum height. The programme will generate and show a plot of the ball’s height 
as a function of time at the user’s discretion.

This programme makes two assumptions:
        
1)	The acceleration due to gravity acting on the ball is constant and equal 
    to the standard acceleration due to gravity, g = -9.81 m/s^2. That is the 
    nominal value for the gravitational acceleration of an object in a vacuum 
    near the surface of the Earth. Therefore, it is assumed that the ball is 
    dropped near the Earth’s surface and restrictions have been placed on the 
    initial height to keep g accurate to a variation of less than 1%.
    
2)	When the ball touches the ground, its energy is multiplied by the 
    coefficient of efficiency, eta. This, along with the near surface condition,
    gives rise to the following relation for the gravitational potientially 
    energy at two neighbouring bounce peaks:
        
        mg(h_n) = eta * mg(h_n-1) 

Kiran Jonathan 15/10/19
"""

import numpy as np

import matplotlib.pyplot as plt

def validated_numeric_input(variable_name, units = '', lower_bound = 0., upper_bound = 2147483647.):

    """
    Requests numeric input from the user and validates it within set boundaries.

    variable_name (string)
    units (string)
    lower_bound (float)
    upper_bound (float)

    upper_bound implicitly defined based on maximum integer value for a 32bit system.

    Kiran Jonathan 15/10/19
    """
        
    numeric_input = lower_bound - 1 
    """Defined to initiate the while loop. A boolean could be used as an
    extra condition of the while loop and its value could be changed once input 
    conditions are met but that would reduce efficiency."""
    
    while numeric_input <= lower_bound or numeric_input >= upper_bound:
            
        input_prompt = ('Please enter a numeric value for'\
                        ' {0} which is greater than {2}{1} but is less'\
                        ' than {3}{1}: ').format(variable_name, units, lower_bound, upper_bound) 
            
        try: numeric_input = float(input(input_prompt))
            
        except: print('\nPlease make sure your input is a number.')
            
    return numeric_input
 
    
def validated_yes_no_input(input_prompt):
    
    """
    Asks the user a yes or no question, input_prompt, and returns a boolean 
    corresponding to the users answer.
    
    input_prompt (string)
    
    Kiran Jonathan 17/10/19
    """
    
    positive_inputs = ['Yes', 'yes', 'Y', 'y']
    
    negative_inputs = ['No', 'no', 'N', 'n']
    
    while True:
        
        user_input = input(input_prompt)
        
        if user_input in positive_inputs:
            
            return True
        
        elif user_input in negative_inputs:
        
            return False
        
        else:
        
            print('\nPlease input yes or no as an answer.')
        

def number_of_bounces(initial_height_local, minimum_height_local, efficiency_coefficient, minimum_height_counts_local):
        
    """
    Calculates the number of bounces above a minimum height given the initial 
    height of release and a coefficient of efficieny. This function will count
    a bounce with a peak height of exactly the minimum height as a bounce over
    the minimum height.

    initial_height_local (float)
    minimum_height_local (float)
    efficiency_coefficient (float)
    minimum_height_counts_local (boolean)

    Kiran Jonathan 15/10/19
    """
    
    current_height = initial_height_local * efficiency_coefficient

    peak_heights_local = []
    
    while current_height >= minimum_height_local:
        
        if minimum_height_counts_local == False and current_height == minimum_height_local:
            
            return len(peak_heights_local), peak_heights_local
        
        peak_heights_local.append(current_height)
        
        current_height = current_height * efficiency_coefficient
        
    number_of_bounces = len(peak_heights_local)
           
    return number_of_bounces, peak_heights_local
 

def total_time_and_distance(peak_heights_local, initial_height_local, gravitational_acceleration = 9.81):
    
    """
    Calculates the time taken for an object to fall from a given starting height 
    and then bounce up and down for a set of given peak heights.

    peak_heights_local (array)
    initial_height_local (float)
    gravitational_acceleration (float)
    
    gravitational_acceleration implicitly defined based on the the standard 
    acceleration due to gravity, g. It is the nominal value for the gravitational 
    acceleration of an object in a vacuum near the surface of the Earth.

    Kiran Jonathan 15/10/19
    """  
    
    total_distance_local = initial_height_local
    
    total_time_local = np.sqrt(2 * initial_height_local / gravitational_acceleration)
    
    peak_times_local = np.array([total_time_local])
    
    for element in peak_heights_local:
        
        bounce_time = 2 * np.sqrt(2 * element / gravitational_acceleration)
    
        total_time_local += bounce_time
        
        peak_times_local = np.append(peak_times_local, total_time_local - 0.5 * bounce_time)
    
        total_distance_local += 2 * element
        
    return total_time_local, total_distance_local, peak_times_local
 
    
def print_results(minimum_height_local, initial_height_local, number_of_bounces_local, \
                  total_time_local, total_distance_local, peak_heights_local, peak_times_local):
    
    """
    This function prints all calculated results. It also calls for the generation
    and showing of a plot of the ball's height against time at the user's discretion.
    
    minimum_height_local (float)
    initial_height_local (float)
    number_of_bounces_local (int)
    total_time_local (float)
    total_distance_local (float)
    peak_heights_local (array)
    peak_times_local (array)
    
    Kiran Jonathan 17/10/19
    """
    
    if number_of_bounces_local != 0:

        print(('\n\nThe ball bounces {0} time(s) above {1}m.'\
               ).format(number_of_bounces_local,minimum_height_local))
    
        print(('\nThe time taken for these bounces is {0:.2f}s.').format(total_time_local))
    
        print('\nThe total distance travelled in this time is {0:.2f}m.'.format(total_distance_local))
    
        if minimum_height_local in peak_heights_local:
        
            print('\nNote: This includes a bounce with a peak exactly at the minimum height.')
    
        if validated_yes_no_input("\nWould you like to see a graph of the ball's "\
                                  "height as a function of time? "):

            if total_time_local/5000 <= 0.1: 
                """Checks if the time divisions of the plot would produce a precise 
                enough to graph. Plots with more points than 5000 were tested and 
                deemed to take too long to produce."""
            
                print("\n\nHere is the plot of the ball's height as a function of time. "\
                      'A reference line of the minimum height is included.')

                plot(total_time_local, peak_heights_local, \
                     peak_times_local, initial_height_local, minimum_height_local)
            
            else:
            
                print('\nUnfortunately, due to the timescale of your values, a'\
                      ' precise plot would take too long to produce. Many apologies.')

    else:
    
        print('\nThe ball will not bounce above the minimum height.')
            
    pass


def plot(total_time_local, peak_heights_local, \
         peak_times_local, initial_height_local, minimum_height_local, gravitational_acceleration = 9.81):

    """
    Produces a plot of the balls height as a function of time.
    
    total_time_local (float)
    peak_heights_local (array)
    peak_times_local (array)
    initial_height_local (float)
    minimum_height_local (float)
    gravitational_acceleration (float)
    
    gravitational_acceleration implicitly defined based on the the standard 
    acceleration due to gravity, g. It is the nominal value for the gravitational 
    acceleration of an object in a vacuum near the surface of the Earth.
    
    Kiran Jonathan 17/10/19
    """
    
    try:
        
        t = np.arange(0, total_time_local, total_time_local/5000)
        
        y = np.array([])
        
        i = 0
        
        for element in t:
            
            while  i <= len(peak_heights_local):
            
                if i == 0:
        
                    y_element = initial_height_local \
                    - 0.5 * gravitational_acceleration * np.power(element,2)
                
                else:
                
                    y_element = peak_heights_local[i-1] \
                    - 0.5 * gravitational_acceleration * np.power(element - peak_times_local[i],2)
                
                if y_element > 0:
                
                    y = np.append(y,y_element)
                    
                    break
                    
                else:
                    
                    y = np.append(y,0)
                    
                    i += 1   
                    
                    break
         
        plt.grid(True)
        plt.plot(t,y)
        plt.plot(t, minimum_height_local * np.power(t,0))
        plt.ylabel('Height / m')
        plt.xlabel('Time / s')
        plt.title("A graph of height as a function of time for the ball:")
        plt.axis([0,total_time_local,0,initial_height_local])
        plt.show()
        
        pass
    
    except: 
        
        print('Unfortunately a plot could not be made.')
        
        pass


print('\nWelcome to PHYS20161 Assignment 1: Bouncy Ball!')
print('\nThis programme calculates the number of times a ball bounces over a '\
      'minimum height given the initial height it is dropped from and the '\
      'coeffecient of efficency for the conservation of gravitational potential '\
      'energy, η. It also calculates the time taken to complete these bounces '\
      'from the dropping of the ball to the touching of the ground after the '\
      'final bounce, as well as the total distance travelled over this time.\n')

minimum_height_counts = validated_yes_no_input('Would you like to include bounces'\
                                               ' with a peak height equal to the minimum height? ')

eta = validated_numeric_input('η', upper_bound = 1.)

initial_height = validated_numeric_input('the initial height', units = 'm', upper_bound = 30000.)
"""upper_bound defined to ensure a variation of less than 1% in g to allow
the assumption of the acceleration due to gravity being constant to remain valid"""

minimum_height = validated_numeric_input('the minimum height', units = 'm', upper_bound = initial_height)

number_of_bounces, peak_heights = number_of_bounces(initial_height, minimum_height, eta, minimum_height_counts)

time_taken, total_distance, peak_times = total_time_and_distance(peak_heights, initial_height)

print_results(minimum_height, initial_height, number_of_bounces, time_taken, \
              total_distance, peak_heights, peak_times)
    
print('\nThank you for using this programme.')
        
        
    