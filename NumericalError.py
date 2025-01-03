import numpy as np 
import matplotlib.pyplot as plt
import datetime as dtm
import AnalyticalSolution
import NumericalSolution
import Calculation


### Initialization ###

beginning_date_and_time = dtm.datetime.now()

## Discretization ##

Lx, Ly = 10*10**(-2), 10*10**(-2)

Nx_min, Nx_max = 5, 100
Nx_step = 20
Nx_list = np.linspace(Nx_min, Nx_max, Nx_step, dtype = int)

dx_list, error_list = [], []

simulation_time = 45
dt = 0.01


## Physical parameters ##

temperature_init = 50
temperature_T1 = 150
temperature_T2 = 100

# For a steel plate
thermal_conductivity = 30 # (in W.m^-1.K^-1)
density = 800  # (in kg.m^-3)
heat_capacity = 520  # (in J.kg^-1.K^-1)

thermal_diffusivity = thermal_conductivity / (density * heat_capacity)  # (in m^2.s^-1)


### Calculation for each meshgrid ###

for Nx in Nx_list: 

    print(Nx)

    Ny = Nx
    dx, dy = Lx/Nx, Ly/Ny

    temperature_init_list = np.zeros((Nx, Ny))*temperature_init

    ## Meshgrid ##

    x_list, y_list = np.linspace(0, Lx, Nx), np.linspace(0, Ly, Ny)


    ## Boundary conditions ##

    temperature_init_list[:, 0] = temperature_T1
    temperature_init_list[:, -1] = temperature_T2
    temperature_init_list[0, :] = temperature_T2
    temperature_init_list[-1, :] = temperature_T2


    ## Simulation parameters ##

    x_fourrier_coeff = thermal_diffusivity*dt/(dx**2)
    y_fourrier_coeff = thermal_diffusivity*dt/(dy**2)


    ### Solutions ###

    ## Analytical solutions ##
    temperature_analytical_list = AnalyticalSolution.solution(x_list, y_list, Lx, Ly, temperature_T1, temperature_T2, temperature_init_list.copy())

    ## Numerical solutions ##
    x_a_lign, x_b_lign, x_c_lign, y_a_lign, y_b_lign, y_c_lign = NumericalSolution.system_matrix(Nx, Ny, x_fourrier_coeff, y_fourrier_coeff)
    numerical_solution_total = NumericalSolution.solution(simulation_time, Nx, Ny, dt, x_fourrier_coeff, y_fourrier_coeff, x_list, y_list, temperature_init_list.copy(), 
                                                        temperature_T1, temperature_T2, x_a_lign, x_b_lign, x_c_lign, y_a_lign, y_b_lign, y_c_lign)

    ## Error ##
    temperature_error = np.linalg.norm(temperature_analytical_list - numerical_solution_total[-1])

    dx_list.append(dx)
    error_list.append(temperature_error)

plt.figure(1)
plt.scatter(dx_list, error_list, marker = "x", color = "black")
plt.title("Erreur en fonction du maillage")
plt.xlabel("h")
plt.ylabel("error")

ln_of_dx_list = []
ln_of_error_list = []

for element in dx_list:
    ln_of_dx_list.append(np.log(element))

for element in error_list:
    ln_of_error_list.append(np.log(element))

plt.figure(2)
plt.scatter(ln_of_dx_list, ln_of_error_list, marker = "x", color = "black")
plt.title("Erreur en fonction du maillage (courbe linéarisée)")
plt.xlabel("ln(h)")
plt.ylabel("ln(error)")



### Runtime ###
Calculation.runtime_program(beginning_date_and_time)

plt.show()