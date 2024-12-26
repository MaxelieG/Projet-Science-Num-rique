import numpy as np 
import matplotlib.pyplot as plt
import datetime as dtm
import AnalyticalSolution
import NumericalSolution
import Calculation
import Graph


### Initialization ###

beginning_date_and_time = dtm.datetime.now()

## Discretization ##

Lx, Ly = 10*10**(-2), 10*10**(-2)
Nx, Ny = 100, 100
dx, dy, = Lx/Nx, Ly/Ny

simulation_time = 0.5
dt = 0.01


## Physical parameters ##

temperature_init = 50
temperature_T1 = 150
temperature_T2 = 100
temperature_init_list = np.zeros((Nx, Ny))*temperature_init

# For a steel plate
thermal_conductivity = 30 # (in W.m^-1.K^-1)
density = 800  # (in kg.m^-3)
heat_capacity = 520  # (in J.kg^-1.K^-1)

thermal_diffusivity = thermal_conductivity / (density * heat_capacity)  # (in m^2.s^-1)


## Meshgrid ##

x_list, y_list = np.linspace(0, Lx, Nx), np.linspace(0, Ly, Ny)


## Boundary conditions ##

temperature_init_list[:, 0] = temperature_T1
temperature_init_list[:, -1] = temperature_T2
temperature_init_list[0, :] = temperature_T2
temperature_init_list[-1, :] = temperature_T2


##Simulation parameters ##

x_fourrier_coeff = thermal_diffusivity*dt/(dx**2)
y_fourrier_coeff = thermal_diffusivity*dt/(dy**2)

### Solutions ###

## Analytical solutions ##

temperature_analytical_list = AnalyticalSolution.solution(x_list, y_list, Lx, Ly, temperature_T1, temperature_T2, temperature_init_list.copy())


## Numerical solutions ##

x_a_lign, x_b_lign, x_c_lign, y_a_lign, y_b_lign, y_c_lign = NumericalSolution.system_matrix(Nx, Ny, x_fourrier_coeff, y_fourrier_coeff)
temperature_numerical_list = NumericalSolution.solution(simulation_time, Nx, Ny, dt, x_fourrier_coeff, y_fourrier_coeff, x_list, y_list, temperature_init_list.copy(), 
                                                        temperature_T1, temperature_T2, x_a_lign, x_b_lign, x_c_lign, y_a_lign, y_b_lign, y_c_lign)


## Difference between analytical and numerical solution ##

temperature_error = np.abs(temperature_analytical_list - temperature_numerical_list)


### Plotting ###
figure = plt.figure("result_heat_2D")

Graph.heatmap_plot(figure, temperature_analytical_list, "Solution analytique", 1)
Graph.heatmap_plot(figure, temperature_numerical_list, "Solution numérique", 2)
Graph.heatmap_plot(figure, temperature_error, "|Analytique - Numérique|", 3)

plt.suptitle("Conduction instationnaire en 2D dans une plaque d'acier (simulation de " + str(simulation_time)+"s)", fontsize = 16)

plt.subplots_adjust(left=0.05, right=0.95, top=0.55, bottom=0.15)

### Runtime ###
Calculation.runtime_program(beginning_date_and_time)



plt.show()