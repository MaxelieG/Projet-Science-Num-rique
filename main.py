import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import datetime as dtm
import AnalyticalSolution
import NumericalSolution
import Calculation
import Graph


### Initialization ###

beginning_date_and_time = dtm.datetime.now()
SAVE_ANIMATION = False

## Discretization ##

Lx, Ly = 10*10**(-2), 10*10**(-2)
Nx, Ny = 100, 100
dx, dy, = Lx/Nx, Ly/Ny

simulation_time = 1
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


## Difference between analytical and numerical solution ##

temperature_error = np.abs(temperature_analytical_list - numerical_solution_total[-1])

## Estimation of the caracteristic time of evolution of the system ##




### Plotting ###

## Animation initialization ##
print("Création de l'animation")
# Création de la figure
fig, ax = plt.subplots()
heatmap = ax.imshow(numerical_solution_total[0], cmap='hot', interpolation='nearest')
ax.set_title("Évolution de la température")
cbar = plt.colorbar(heatmap, ax=ax)  # Ajouter une barre de couleur
cbar.set_label('Température (°C)')

# Fonction de mise à jour pour l'animation
def update(frame):
    data = numerical_solution_total[frame]
    heatmap.set_array(data)
    
    # Mise à jour de l'échelle de la heatmap (et donc de la colorbar)
    heatmap.set_clim(vmin=data.min(), vmax=data.max())
    ax.set_title(f"Évolution de la température (time =  {(frame*dt):.2f}s)")
    return [heatmap]

ani = FuncAnimation(fig, update, frames=int(simulation_time/dt), blit=True, interval=10)

if SAVE_ANIMATION:
    print("Enregistrement de l'animation")
    ani.save('evolution_temperature.mp4', fps=10, extra_args=['-vcodec', 'libx264'])
    print("Animation sauvegardée")

# ----------------------------- #

figure = plt.figure("result_heat_2D", constrained_layout=True)

Graph.heatmap_plot(figure, temperature_analytical_list, "Solution analytique", 1)
Graph.heatmap_plot(figure, numerical_solution_total[-1], "Solution numérique", 2)
Graph.heatmap_plot(figure, temperature_error, "|Analytique - Numérique|", 3)

plt.suptitle("Conduction instationnaire en 2D dans une plaque d'acier (simulation de " + str(simulation_time)+"s)", fontsize = 16)

plt.subplots_adjust(left=0.05, right=0.95, top=0.55, bottom=0.15)


### Runtime ###
Calculation.runtime_program(beginning_date_and_time)

plt.show()