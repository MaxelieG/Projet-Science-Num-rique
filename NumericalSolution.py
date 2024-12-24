import numpy as np 
import matplotlib.pyplot as plt
from TDMA import TDMAsolver
import datetime as dtm


### Initialization ###

beginning_date_and_time = dtm.datetime.now()

## Runtime managing ##

def runtime_program (beginning_date_and_time):

    """Calculates and prints program runtime based on start time.

    Parameters
    ----------
    beginning_date_and_time : datetime.datetime
    
    Returns
    -------
    None
    """

    ending_date_and_time=dtm.datetime.now()
    runtime = ending_date_and_time - beginning_date_and_time

    #type(runtime) =  datetime.timedelta, it should be convert in a string to be easily manipulated and printed
    string_runtime = str(runtime)

    hour = string_runtime.split(":")[0]
    minute = string_runtime.split(":")[1]
    second = str(round(float(string_runtime.split(":")[2]),2))

    # These different cases are made to avoid printing 00h00min00,4s for the shortest program (just for the visual aspect)
    if hour == "0":

        if minute == "00":

            print("Runtime : "+second+"s")

        else:

            print("Runtime : "+minute+"min"+second+"s")
    
    else:

        print("Runtime : "+hour+"h"+minute+"min"+second+"s")
    
    return None


## Discretization parameters ##

Lx, Ly = 10*10**(-2), 10*10**(-2)
Nx, Ny = 100, 100
dx, dy, = Lx/Nx, Ly/Ny
simulation_time = 60
dt = 0.01

## Physical parameters ##

temperature_init = 50
temperature_T1 = 150
temperature_T2 = 100

# For a steel plate
thermal_conductivity = 30 # (in W.m^-1.K^-1)
density = 8000  # (in kg.m^-3)
heat_capacity = 520  # (in J.kg^-1.K^-1)

thermal_diffusivity = thermal_conductivity / (density * heat_capacity)  # (in m^2.s^-1)

temperature_numerical_list = np.ones((Nx, Ny))*temperature_init

## Meshgrid ##

x_list, y_list = np.linspace(0, Lx, Nx), np.linspace(0, Ly, Ny)

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

    ## Boundary conditions ##

    temperature_numerical_list[:, 0] = temperature_T1
    temperature_numerical_list[:, -1] = temperature_T2
    temperature_numerical_list[0, :] = temperature_T2
    temperature_numerical_list[-1, :] = temperature_T2

    for y in range(1, Ny-1):

        x_dy_lign = [(1-y_fourrier_coeff)*temperature_numerical_list[x, y] + 
                     y_fourrier_coeff*(temperature_numerical_list[x, y-1] + temperature_numerical_list[x, y+1])/2 for x in range(len(x_list))]
        
        temperature_numerical_list[:, y] = TDMAsolver(x_a_lign.copy(), x_b_lign.copy(), x_c_lign.copy(), x_dy_lign.copy())

    for x in range(1, Nx-1):

        y_dx_lign = [(1-x_fourrier_coeff)*temperature_numerical_list[x, y] + 
                     x_fourrier_coeff*(temperature_numerical_list[x-1, y] + temperature_numerical_list[x+1, y])/2 for y in range(len(y_list))]
        
        temperature_numerical_list[x, :] = TDMAsolver(y_a_lign.copy(), y_b_lign.copy(), y_c_lign.copy(), y_dx_lign.copy())

    # Update time
    time += dt
    Nt += 1

print(temperature_numerical_list)
plt.matshow(temperature_numerical_list, cmap = "hot")
plt.title("Solution numérique")
plt.colorbar()


#Runtime
runtime_program(beginning_date_and_time)

plt.show()

