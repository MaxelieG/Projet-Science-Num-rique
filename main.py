import numpy as np 
import matplotlib.pyplot as plt


### Initialization ###

## Discretization parameters ##

Lx, Ly = 1, 1
Nx, Ny = 100, 100
dx, dy, = Lx/Nx, Ly/Ny
dt = 0

## Physical parameters ##

temperature_init = 50
temperature_T1 = 150
temperature_T2 = 100
thermal_diffusivity = 0
temperature_list = np.zeros((Nx, Ny))

## Meshgrid ##

x_list, y_list = np.linspace(0, Lx, Nx), np.linspace(0, Ly, Ny)


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

for i, x in enumerate(x_list):

    for j, y in enumerate(y_list):

        temperature = temperature_T2 + 4*(temperature_T1 - temperature_T2)*infinite_sum(x, y)/np.pi 

        temperature_list[i, j] = temperature


plt.matshow(temperature_list, cmap = "hot")
plt.colorbar()
plt.show()