import numpy as np 
from Calculation import TDMAsolver


### Calculation of the matrix for the two linear systems created by implicit scheme ###

def system_matrix(Nx, Ny, x_fourrier_coeff, y_fourrier_coeff):

    x_a_lign = np.concatenate([np.ones((Nx-1))*(-x_fourrier_coeff/2), np.array([0])])
    x_b_lign = np.concatenate([np.array([1]), np.ones((Nx-2))*(1+x_fourrier_coeff), np.array([1])])
    x_c_lign = np.concatenate([np.array([0]), np.ones((Nx-1))*(-x_fourrier_coeff/2)])

    y_a_lign = np.concatenate([np.ones((Ny-1))*(-y_fourrier_coeff/2), np.array([0])])
    y_b_lign = np.concatenate([np.array([1]), np.ones((Ny-2))*(1+y_fourrier_coeff), np.array([1])])
    y_c_lign = np.concatenate([np.array([0]), np.ones((Ny-1))*(-y_fourrier_coeff/2)])

    return x_a_lign, x_b_lign, x_c_lign, y_a_lign, y_b_lign, y_c_lign


### Calculation of the numerical solution ###

def solution(simulation_time, Nx, Ny, dt, x_fourrier_coeff, y_fourrier_coeff, x_list, y_list, temperature_list, temperature_T1, temperature_T2, 
             x_a_lign, x_b_lign, x_c_lign, y_a_lign, y_b_lign, y_c_lign):

    #Initialization
    time = 0
    Nt = 0
    numerical_solution_total = []

    while time < simulation_time:

        # Print current time
        print(f"Time = {time:.2f} s")

        # Boundary conditions
        temperature_list[:, 0] = temperature_T1
        temperature_list[:, -1] = temperature_T2
        temperature_list[0, :] = temperature_T2
        temperature_list[-1, :] = temperature_T2

        # Iteration on the grid
        for y in range(1, Ny-1): #t = n + 1/2

            x_dy_lign = [(1-y_fourrier_coeff)*temperature_list[x, y] + 
                        y_fourrier_coeff*(temperature_list[x, y-1] + temperature_list[x, y+1])/2 for x in range(len(x_list))]
            
            temperature_list[:, y] = TDMAsolver(x_a_lign.copy(), x_b_lign.copy(), x_c_lign.copy(), x_dy_lign.copy())

        for x in range(1, Nx-1): #t = n + 1

            y_dx_lign = [(1-x_fourrier_coeff)*temperature_list[x, y] + 
                        x_fourrier_coeff*(temperature_list[x-1, y] + temperature_list[x+1, y])/2 for y in range(len(y_list))]
            
            temperature_list[x, :] = TDMAsolver(y_a_lign.copy(), y_b_lign.copy(), y_c_lign.copy(), y_dx_lign.copy())

        # Update time
        time += dt
        Nt += 1
        numerical_solution_total.append(temperature_list.copy())

    return numerical_solution_total

### Estimation of the time to reach the stationnary state ###

def caracteristic_time (all_temperature, dt):

    previous_temperature = all_temperature[0]

    for index, temperature in enumerate(all_temperature[1:]):

        if np.mean(np.abs(previous_temperature - temperature)) < 10**(-5):

            time = index*dt
            
            return time

        previous_temperature = temperature
    
    return None

