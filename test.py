import numpy as np
import matplotlib.pyplot as plt
from TDMA import TDMAsolver  

# Physical parameters
thermal_conductivity = 30  # (in W.m^-1.K^-1)
density = 8000  # (in kg.m^-3)
heat_capacity = 520  # (in J.kg^-1.K^-1)

thermal_diffusivity = thermal_conductivity / (density * heat_capacity)  # (in m^2.s^-1)

# Discretization
Lx = 10 * 10**(-2)  # Length of the domain
Nx = 100  # Number of grid points
dx = Lx / Nx  # Grid spacing

# Temperatures
temperature_init = 80 + 273  # Initial temperature (in Kelvin)
temperature_water = 20 + 273  # Boundary temperature (in Kelvin)


### Initialization ###

temperature = np.ones((Nx))*temperature_init
temperature[0], temperature[-1] = temperature_water, temperature_water

simulation_time = 5 #Time of simulation
dt = 0.0005 #Time steps to respect stability criteria
print(f"Time step: {dt}")


fourrier_coeff = thermal_diffusivity*dt/(dx**2)

c = np.concatenate([np.array([0]), np.ones((Nx-1))*(-fourrier_coeff)])
b = np.concatenate([np.array([1]), np.ones((Nx-2))*(1+2*fourrier_coeff), np.array([1])])
a = np.concatenate([np.ones((Nx-1))*(-fourrier_coeff), np.array([0])])

time = 0
Nt = 0


# Iteration loop

plt.figure("Profil de température")
plt.plot(np.linspace(0, Nx, Nx), temperature, label = str(time))
while time < simulation_time:
    # Print current time
    print(f"Time = {time:.2f} s")


    #print(Nt)

    temperature = TDMA.TDMAsolver(a.copy(), b.copy(), c.copy(), temperature.copy())
    #temperature[0], temperature[-1] = temperature_water, temperature_water
    


    # Update time
    time += dt
    Nt += 1


    if Nt%150 == 0:
        plt.plot(np.linspace(0, Nx, Nx), temperature, label = str(time))


plt.show()
