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
temperature_T2 = 150
thermal_diffusivity = 0
temperature_list = np.zeros((Nx, Ny))

## Meshgrid ##

x_list, y_list = np.linspace(0, Lx, Nx, dtype=int), np.linspace(0, Ly, Ny, dtype=int)

X, Y = np.meshgrid(x_list, y_list)

## Definition of the infinite sum ##

def term_the_sum (n, x, y):

    n_term = (1/(2*n + 1))*np.sin((2*n + 1)*np.pi*x/Lx)*np.sinh((Ly - y)*(2*n + 1)*np.pi/Lx)/np.sinh((2*n + 1)*np.pi*Ly/Lx)

    return n_term

def infinite_sum (x, y):

    index = 0
    sum = 0

    while abs(term_the_sum(index, x, y) - term_the_sum(index + 1, x, y)) > 10**(-10):

        sum += term_the_sum(index, x, y)
        index += 1
    
    return sum

### Calculation of the analytical solution ###

for x in x_list:

    for y in y_list:

        temperature = temperature_T2 + 4*(temperature_T1 - temperature_T2)*infinite_sum(x, y)/np.pi 

        temperature_list[x, y] = temperature

plt.matshow(temperature_list)
plt.show()