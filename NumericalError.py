import numpy as np 
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score
import datetime as dtm
import AnalyticalSolution
import NumericalSolution
import Calculation


### Initialization ###

beginning_date_and_time = dtm.datetime.now()

## Discretization ##

Lx, Ly = 10*10**(-2), 10*10**(-2)

Nx_analytical = 240
Nx_list = [2, 3, 4, 5, 6, 8, 10, 12, 15, 16, 20, 24, 30, 40, 48, 60, 80, 120]


dx_list, error_list = [], []
ln_of_dx_list, ln_of_error_list = [], []

simulation_time = 45
dt = 0.01


## Physical parameters ##

temperature_init = 50
temperature_T1 = 150
temperature_T2 = 100

# For a steel plate
thermal_conductivity = 30 # (in W.m^-1.K^-1)
density = 8000  # (in kg.m^-3)
heat_capacity = 520  # (in J.kg^-1.K^-1)

thermal_diffusivity = 10 * thermal_conductivity / (density * heat_capacity)  # (in m^2.s^-1)


### Calculation of the "continuous" analytical solution (ie. 100x100)

x_list_analytical, y_list_analytical = np.linspace(0, Lx, Nx_analytical), np.linspace(0, Ly, Nx_analytical)
temperature_init_list_analytical = np.zeros((Nx_analytical, Nx_analytical))*temperature_init
temperature_analytical_list = AnalyticalSolution.solution(x_list_analytical, y_list_analytical, Lx, Ly, temperature_T1, temperature_T2, temperature_init_list_analytical.copy())


### Calculation of the numerical solution for each meshgrid ###

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


    ## Solutions ##

    x_a_lign, x_b_lign, x_c_lign, y_a_lign, y_b_lign, y_c_lign = NumericalSolution.system_matrix(Nx, Ny, x_fourrier_coeff, y_fourrier_coeff)
    numerical_solution_total = NumericalSolution.solution(simulation_time, Nx, Ny, dt, x_fourrier_coeff, y_fourrier_coeff, x_list, y_list, temperature_init_list.copy(), 
                                                        temperature_T1, temperature_T2, x_a_lign, x_b_lign, x_c_lign, y_a_lign, y_b_lign, y_c_lign)
    
    ## Matrix size adaptation ##

    final_temperature = numerical_solution_total[-1]
    reshaped_temperature = np.zeros_like(temperature_analytical_list)

    ratio = int(Nx_analytical/Nx)
    
    for i in range(len(reshaped_temperature)):
        for j in range(len(reshaped_temperature[0])):
            
            i_final = i//ratio
            j_final = j//ratio

            reshaped_temperature[i][j] = final_temperature[i_final][j_final]

    ## Error ##
    temperature_error = Calculation.norm_L2(temperature_analytical_list - reshaped_temperature)

    dx_list.append(dx)
    error_list.append(temperature_error)
    ln_of_dx_list.append(np.log(dx))
    ln_of_error_list.append(np.log(temperature_error))

### Linear Regression ###

X = np.array(ln_of_dx_list).reshape(-1, 1)
Y = np.array(ln_of_error_list)
scaler = StandardScaler()
#X = scaler.fit_transform(X)
#Y = scaler.fit(Y)
model = LinearRegression()
model.fit(X, Y)

a = model.coef_[0]
b = model.intercept_

y_pred = a*X +b

Y = Y.flatten()
y_pred = y_pred.flatten()

R2 = r2_score(Y, y_pred)

### Plot ###

plt.figure(1)
plt.scatter(dx_list, error_list, marker = "x", color = "black")
plt.title("Erreur en fonction du maillage")
plt.xlabel("h")
plt.ylabel("error")

plt.figure(2)
plt.scatter(ln_of_dx_list, ln_of_error_list, marker = "x", color = "black")
plt.plot(ln_of_dx_list, a*np.array(ln_of_dx_list) + b, color = "red", label = "ln(error) = " + str(int(a*100)/100) + " * ln(h) + " + \
         str(int(b*100)/100) + " (R^2 = " + str(int(R2*100)/100) + ")")
plt.title("Erreur en fonction du maillage (courbe linéarisée)")
plt.xlabel("ln(h)")
plt.ylabel("ln(error)")
plt.legend(title = "Linear Regression")


### Runtime ###
Calculation.runtime_program(beginning_date_and_time)

plt.show()