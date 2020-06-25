import numpy as np

U_0 =  1   # [m/s]  Velocty of the uniform flow
L   =  1   # [m]    Length of the domain in each direction
Nx  = 16   # []     Number of cells in each direction
#FOS = 10  #        Factor of safety for particles to pass though cell
N_flow_through_times = 1.5 # Number of times particles will traverse domain
restarts_per_flow_though = 30

# Fluid Properties
gamma = 1.4
specificGasConstant = 287.058  # [J/(kg*K)]
T_max = 300.0                  # [K]

# Speek of sound
c = np.sqrt(gamma*specificGasConstant*T_max)
dt = (L/Nx)/c

## Time step based on particles
#dt = tau/FOS # [s]

dx = L/Nx
t_flow_through = L/U_0
t_final = t_flow_through*N_flow_through_times 

# Time scale for flow and particles to pass though cell
tau_p = dx/U_0 # [s] 

copy_every_time_steps = int(tau_p/dt)

Nt = int(t_final/dt)
restartEveryTimeSteps = int((t_flow_through/dt)/restarts_per_flow_though) 

print('maxIter = {}'.format(Nt))
print('fixedDeltaTime = {}'.format(dt))
print('restartEveryTimeSteps = {}'.format(restartEveryTimeSteps))
print('copyEveryTimeSteps = {}'.format(copy_every_time_steps))
