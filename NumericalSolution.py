import numpy as np 
import matplotlib.pyplot as plt
from TDMA import TDMAsolver


### Initialization ###

## Discretization parameters ##

Lx, Ly = 100, 100
Nx, Ny = 101, 101
dx, dy, = Lx/Nx, Ly/Ny
simulation_time = 600
dt = 0.1

## Physical parameters ##

temperature_init = 50
temperature_T1 = 150
temperature_T2 = 100

# For a steel plate
thermal_conductivity = 3000000  # (in W.m^-1.K^-1)
density = 500  # (in kg.m^-3)
heat_capacity = 520  # (in J.kg^-1.K^-1)

thermal_diffusivity = thermal_conductivity / (density * heat_capacity)  # (in m^2.s^-1)

temperature_numerical_list = np.ones((Nx, Ny))*temperature_init

## Meshgrid ##

x_list, y_list = np.linspace(0, Lx, Nx, dtype=int), np.linspace(0, Ly, Ny, dtype=int)

## Boundary conditions ##

temperature_numerical_list[:, 0] = temperature_T1
temperature_numerical_list[:, -1] = temperature_T2
temperature_numerical_list[0, :] = temperature_T2
temperature_numerical_list[-1, :] = temperature_T2

## Definition of the matrix ##

x_fourrier_coeff = thermal_diffusivity*dt/(dx**2)
y_fourrier_coeff = thermal_diffusivity*dt/(dy**2)

x_a_lign = np.concatenate([np.ones((Nx-1))*(-x_fourrier_coeff/2), np.array([0])])
x_b_lign = np.concatenate([np.array([1]), np.ones((Nx-2))*(1+x_fourrier_coeff), np.array([1])])
x_c_lign = np.concatenate([np.array([0]), np.ones((Nx-1))*(-x_fourrier_coeff/2)])

y_a_lign = np.concatenate([np.ones((Ny-1))*(-y_fourrier_coeff/2), np.array([0])])
y_b_lign = np.concatenate([np.array([1]), np.ones((Ny-2))*(1+y_fourrier_coeff), np.array([1])])
y_c_lign = np.concatenate([np.array([0]), np.ones((Ny-1))*(-y_fourrier_coeff/2)])

### Iteration loop ###

time = 0
Nt = 0


while time < simulation_time:

    # Print current time
    print(f"Time = {time:.2f} s")

    for y in range(1, len(y_list)-1):

        x_dy_lign = [(1-y_fourrier_coeff)*temperature_numerical_list[x, y] + 
                     y_fourrier_coeff*(temperature_numerical_list[x, y-1] + temperature_numerical_list[x, y+1])/2 for x in range(len(x_list))]
        
        temperature_numerical_list[:, y] = TDMAsolver(x_a_lign.copy(), x_b_lign.copy(), x_c_lign.copy(), x_dy_lign.copy())

    for x in range(1, len(x_list)-1):

        y_dx_lign = [(1-x_fourrier_coeff)*temperature_numerical_list[x, y] + 
                     x_fourrier_coeff*(temperature_numerical_list[x-1, y] + temperature_numerical_list[x+1, y])/2 for y in range(len(y_list))]
        
        temperature_numerical_list[x, :] = TDMAsolver(y_a_lign.copy(), y_b_lign.copy(), y_c_lign.copy(), y_dx_lign.copy())

    # Update time
    time += dt
    Nt += 1

    
plt.matshow(temperature_numerical_list, cmap = "hot")
plt.title("Solution numérique")
plt.colorbar()

plt.show()

