import numpy as np 
import matplotlib.pyplot as plt


### Initialization ###

## Discretization parameters ##

Lx, Ly = 100, 100
Nx, Ny = 101, 101
dx, dy, = Lx/Nx, Ly/Ny
dt = 0

## Physical parameters ##

temperature_init = 50
temperature_T1 = 100
temperature_T2 = 0
thermal_diffusivity = 0
temperature_list = np.zeros((Nx, Ny))

## Meshgrid ##

x_list, y_list = np.linspace(0, Lx, Nx, dtype=int), np.linspace(0, Ly, Ny, dtype=int)

X, Y = np.meshgrid(x_list, y_list)

for x in x_list:

    for y in y_list:



        temperature = temperature_T2 + 4*(temperature_T1 - temperature_T2)/np.pi

        sum_value = 0
        rest = 10000
        index = 0
        stop_condition = 1

        while stop_condition > 10**(10):

            stop_condition = abs((1/(2*index + 1))*np.sin((2*index + 1)*np.pi*x/Lx)*np.sinh((Ly - y)*(2*index + 1)*np.pi/Lx)/np.sinh((2*index + 1)*np.pi*Ly/Lx))
            sum_value += (1/(2*index + 1))*np.sin((2*index + 1)*np.pi*x/Lx)*np.sinh((Ly - y)*(2*index + 1)*np.pi/Lx)/np.sinh((2*index + 1)*np.pi*Ly/Lx)
            index += 1

        temperature += sum_value

        temperature_list[x, y] = temperature



plt.matshow(temperature_list)
plt.show()