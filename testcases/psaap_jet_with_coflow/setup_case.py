import argparse
import numpy as np
import json
import sys

parser = argparse.ArgumentParser()
parser.add_argument('json_template', type=argparse.FileType('r'))
args = parser.parse_args()

json_template = json.load(args.json_template)

###############################################################################
# User input
###############################################################################

# Domain
N_nodes_hit = 2
N_nodes_inflow_outflow = 2
N_partitions_x_per_node = 4

# Domain
L = 0.16   # [m]
W_interior = 0.0462 # [m]
streamwise_buffer_percentage = 0.50 # Add this much space on each side to mimic the coflow

# Fluid properties
rho_0 = 1.2                 # [kg/m^3]
mu = 1.5*np.power(10,-5.0)  # [Pa*s] a.k.a [(kg*m)/s]
c_p_fluid = 1012            # [J/ (kg K)]
c_v_fluid =  723            # [J/ (kg K)]
T_0 = 300.0                 # [K]
Pr = 0.707                  # []
#g = 9.81                   # [m/s^2]
g = 0.0                     # [m/s^2] # Gravity off to test particle outrunning stencil

# Flow
rhoU_0 = 9.6               # [kg*m/s]
# From monana and gianluca paper
u_rms_to_U_0_ratio = 0.12
#u_rms_to_U_0_ratio = 0.02 # Less intense fluctutations for debugging
# From monana and gianluca paper
l_to_W_ratio = 0.45

# Particle properties
d_p = 12e-6               # [m]
#d_p = 12.197e-6           # [m]
rho_p = 8900.0              # [kd/m^3]
# from "core" of soleil-mpi simulation
particle_number_density = 27044485758 # [particles/m^3]
restitution_coefficient = 1.0
Nu = 2.0
c_v_particle = 440.0      # [J/(kg*K)]

maxSkew = 1.05
escapeRatioPerDir = 0.05

# For determining the time step
# Maximum expected temperature
T_max = 1000.0             # [K]
# CFL for speed of sound
CFL = 0.50

# Time Integration
N_flow_through_times = 1.0
restarts_per_flow_though_time = 20

## Max step size in mesh in Soleil-MPI
#dx = 4.5*np.power(10.0,-4.0) # [m]
#dy = 4.5*np.power(10.0,-4.0) # [m]
#dz = 4.5*np.power(10.0,-4.0) # [m]

# Average step size in mesh in Soleil-MPI 
dx = 5.0*np.power( 10.0,-4.0) # [m]
dy = 1.25*np.power(10.0,-4.0) # [m]
dz = 1.25*np.power(10.0,-4.0) # [m]

# Average step size in mesh in Soleil-MPI double x resolution
#dx = 2.5*np.power( 10.0,-4.0) # [m]
#dy = 1.25*np.power(10.0,-4.0) # [m]
#dz = 1.25*np.power(10.0,-4.0) # [m]

## Average step size in mesh in Soleil-MPI quadruple x resolution
#dx = 1.25*np.power(10.0,-4.0) # [m]
#dy = 1.25*np.power(10.0,-4.0) # [m]
#dz = 1.25*np.power(10.0,-4.0) # [m]

# Radiation
Q_a = 0.396
Q_s = 0.745
N_angles = 50
# Total amount of radiation entering the domain
Q = 10.0*np.power(10.0,3.0) # [W]
#Q = 1.0 # [W]
DOM_coarsening_factor = 1

###############################################################################
# Dependent Parameters
###############################################################################

N_partitions_x_hit = int(N_nodes_hit*N_partitions_x_per_node)
N_partitions_x_inflow_outflow = int(N_nodes_inflow_outflow*N_partitions_x_per_node)

Ny_interior = int(W_interior/dy)
Nz_interior = int(W_interior/dz)

#while True:
#  if (Ny_interior %  DOM_coarsening_factor) == 0:
#    break
#  else:
#    Ny_interior+=1


Ny_buffer_per_side = int((W_interior*streamwise_buffer_percentage)/dy)
W_buffer_per_side = Ny_buffer_per_side*dy
W_inflow_outflow = W_interior + 2*W_buffer_per_side 

# Domain
HIT_fluid_volume = W_interior*W_interior*W_interior
inflow_outflow_fluid_volume = L*W_inflow_outflow*W_inflow_outflow

# Grid
Nx_hit = int(W_interior/dx)
Ny_hit = int(W_interior/dy)
Nz_hit = int(W_interior/dz)

Nx_inflow_outflow = int(L/dx)
Ny_inflow_outflow = int(W_inflow_outflow/dy)
Nz_inflow_outflow = int(W_inflow_outflow/dz)

## Hard coded
#Nx_hit = 92
#Nx_inflow_outflow = 92*4
#Ny = 92
#Nz = 92
#dx = L/Nx_inflow_outflow
#dy = W/Ny
#dz = W/Nz

# Add cells to the x direction until evenly divisible by the number of partitions
while True:
  if (Nx_inflow_outflow % N_partitions_x_inflow_outflow) != 0:
    Nx_inflow_outflow += 1
  else:
    break


# Fluid
nu = mu/rho_0
gamma = c_p_fluid/c_v_fluid
specificGasConstant = c_v_fluid*(gamma-1.0)
k_fluid = c_p_fluid*mu/Pr
# Initial Pressure so that initial temperature is what you asked for
P_0 = T_0*specificGasConstant*rho_0  # [Pa]

# Flow
# Inlet velocity 
U_0 = rhoU_0/rho_0 
# Bulk Reynolds Number
Re_W = rho_0*U_0*W_interior/mu

# Mean rms velcoity fluctuation (in HIT section)
u_rms = u_rms_to_U_0_ratio*U_0

# Large-eddy length scale (HIT)
l_0 = l_to_W_ratio*W_interior   # [m]

# Large eddy reynolds number
Re_L = rho_0*u_rms*l_0/mu

# Taylor Renyolds Number (from Pope chapter 6 equation 6.64 in the )
Re_lambda = np.sqrt((20.0/3.0)*Re_L)

# Taylor Length Scale
lambda_ = ( Re_lambda*mu )/( rho_0*u_rms )

# Mean KTE dissipation rate (HIT)
epsilon = ( 15.0*mu*( u_rms**2.0 ) )/( rho_0*( lambda_**2.0 ) )                         # [m^2/s^3]

# Kolmogorov length scale (HIT)
eta = ( ( mu**3.0 )/( ( rho_0**3.0 )*epsilon ) )**( 1.0/4.0 )                           # [m]

# Kolmogorov time scale
tau_eta = np.sqrt( nu/epsilon )                                               # [s]

# Kolmogorov velocity scale
u_eta = np.power(nu*epsilon, 1.0/4.0)                                               # [m/s]

dx_dns = (np.pi*np.sqrt(2)/3)*eta
dy_dns = (np.pi*np.sqrt(2)/3)*eta
dz_dns = (np.pi*np.sqrt(2)/3)*eta
Nx_hit_dns = int(W_interior/dx_dns)
Ny_hit_dns = int(W_interior/dy_dns)
Nz_hit_dns = int(W_interior/dz_dns)

Nx_inflow_outflow_dns = int(L/dx_dns)
Ny_inflow_outflow_dns = int(W_inflow_outflow/dy_dns)
Nz_inflow_outflow_dns = int(W_inflow_outflow/dz_dns)

## Large-eddy velocity (HIT)
#u_0 = ( ( mu/rho_0 )*( Re_lambda**2.0 ) )/( 15.0*l_0 )                                  # [m/s]
u_0 = u_rms # APPROXIMATION GO BACK AND FIX LATER                                 # [m/s]

## Large-eddy time scale (HIT)
t_0 = l_0/u_0                                                                           # [s]

## Turbulent kinetic energy (HIT)
k_0 = ( 3.0/2.0 )* np.power(u_0,2.0)                                                           # [m^2/s^2]

# Speed of sound in air
#c = 331.5 + 0.6*( T_max - 273.15 )
#print('c_old = {}'.format(c))
c = np.sqrt(gamma*specificGasConstant*T_max)

# Flow acoustic time step 
delta_t_c = CFL*(min(min(dx,dy),dz)/c)

# Particle relaxation time: time required for a particle to adapt to the flow eddies
tau_p = (rho_p*np.power(d_p,2.0))/(18.0*mu)                                                                  # [s]

# Particle thermal relaxation time: time required for the thermal energy of a particle to adapt to the flow
tau_th = ( 3.0*tau_p*c_v_particle*Pr )/( Nu*c_p_fluid )                                              # [s]

# Ratio between minimal tau scales and acoustic time step
FoS = 1.0/5.0 # factor of safty in ratio of scales
ratio_delta_t = FoS*( min( tau_eta, tau_p, tau_th )/delta_t_c )                # [-]

## Frequency copy between domains in number of flow time step
#ratio_copy = (dx/(U_0))/delta_t_c                                                # [-]
ratio_copy = (dx/(U_0+u_rms))/delta_t_c                                                # [-]

Re_p = rho_0*u_rms*d_p/mu

## Churchill-Bernstein Equation
#Nu = 0.3 + (0.62*np.power(Re_p,0.5)*np.power(Pr,1.0/3.0))/np.power(1 + np.power(0.4/Pr, 2.0/3.0),1.0/4.0) * np.power(1 + np.power(Re_p/282000,5.0/8.0),4.0/5.0)
#print('Pr*Re = {}'.format(Pr*Re_p))
#print('Nu = {}'.format(Nu))

#########################
# Particles
#########################
# particle volume
V_p = (np.pi*np.power(d_p,3.0))/6.0

# particle mass
m_p = rho_p*V_p

# Number of particles
Np_HIT = int(particle_number_density * HIT_fluid_volume)
Np_HIT_max = int(1.3*Np_HIT)

Np_inflow_outflow = int(particle_number_density * inflow_outflow_fluid_volume)
Np_inflow_outflow_max = int(1.3*Np_inflow_outflow)

# Add particles until there can be an even partitioning of them
while True:
  if (Np_HIT  % N_partitions_x_per_node) == 0:
    break
  else:
    Np_HIT += 1

while True:
  if (Np_HIT_max  % N_partitions_x_per_node) == 0:
    break
  else:
    Np_HIT_max += 1

# Add particles until there can be an even partitioning of them
while True:
  if (Np_inflow_outflow % N_partitions_x_per_node) == 0:
    break
  else:
    Np_inflow_outflow += 1

while True:
  if (Np_inflow_outflow_max % N_partitions_x_per_node) == 0:
    break
  else:
    Np_inflow_outflow_max += 1

# particle convection coefficient 
h = k_fluid*Nu/d_p # Based on the definition of the Nusselt number


#########################
# Time
#########################
dt = delta_t_c 
t_flow_through = L/U_0
Nt_per_flow_through_time = int(np.ceil(t_flow_through/dt))
t_final = N_flow_through_times*t_flow_through
Nt_between_restarts = int(np.ceil(Nt_per_flow_through_time/restarts_per_flow_though_time))
Nt = int(np.ceil(t_final/dt))

#########################
# Particle Spread for Co-flow Buffer
#########################
W_spread = u_rms * t_flow_through 
W_spread_percent = W_spread / W_interior

#########################
# Radiation
#########################
# Radiative Heat Flux at Wall
q_0 = Q/(L*W_interior)
# 'Intensity' at wall this is not really intensity 
Intensity = q_0

#########################
# Memory
#########################
doubles_per_fluid_cell = 70
doubles_per_dom_cell   = 12
doubles_per_particle   = 35

#########################
# Copy
#########################
### Ratio between minimal tau scales and acoustic time step
#ratio_delta_t = FoS*( min( tau_eta, tau_p, tau_th, tau_rad )/delta_t_c )                # [-]
#output( 'Ratio scales to c time step: \t  ratio     =', ratio_delta_t, '\t \t     [-]' )

print('ratio_delta_t = {}'.format(ratio_delta_t))
print('ratio_copy = {}'.format(ratio_copy))
# Adjust stagger factor to be compatible with copy frequency
if ratio_delta_t > ratio_copy:
    ratio_copy = int(ratio_copy)
    ratio_delta_t = ratio_copy
else:
    factor = math.ceil(ratio_copy / ratio_delta_t)
    ratio_delta_t = math.ceil(ratio_copy / factor)
    ratio_copy = ratio_delta_t * factor
print('ratio_delta_t = {}'.format(ratio_delta_t))
print('ratio_copy = {}'.format(ratio_copy))

#print('staggerFactor = {}'.format(ratio_delta_t))
#print('copyEveryTimeSteps= {}'.format(ratio_copy))

#########################
# Hand tweak stuff
#########################


#if (Nx_hit % 2) != 0:
# Nx_hit = Nx_hit+1
#
#if (Nx_inflow_outflow % 2) != 0:
# Nx_inflow_outflow = Nx_inflow_outflow+1
#
### round the number of cells to the nearest paritionsize
##if Nx_inflow_outflow % N_partitions_x != 0:
##  Nx_inflow_outflow  = int(np.ceil(Nx_inflow_outflow / N_partitions_x )) * N_partitions_x 
##  print('Rounded to the nearest number cells divisible by the number of partitions')
#
#if (Ny % 2) != 0:
# Ny = Ny+1
#if (Nz % 2) != 0:
# Nz = Nz+1
#
#Nx_DOM = int(Nx_inflow_outflow/2)
#Ny_DOM = int(Ny/2)
#Nz_DOM = int(Nz/2)

if (Nx_inflow_outflow %  DOM_coarsening_factor) != 0:
  print('ERROR')
  print('Nx_inflow_outflow ({}) is not evenly divisible by DOM_coarsening_factor ({})'.format(Nx_inflow_outflow, DOM_coarsening_factor))
  sys.exit()

if (Ny_inflow_outflow %  DOM_coarsening_factor) != 0:
  print('ERROR')
  print('Ny_inflow_outflow ({}) is not evenly divisible by DOM_coarsening_factor ({})'.format(Ny_inflow_outflow, DOM_coarsening_factor))
  sys.exit()

if (Nz_inflow_outflow %  DOM_coarsening_factor) != 0:
  print('ERROR')
  print('Nz_inflow_outflow ({}) is not evenly divisible by DOM_coarsening_factor ({})'.format(Nz_inflow_outflow, DOM_coarsening_factor))
  sys.exit()


Nx_DOM = int(Nx_inflow_outflow/DOM_coarsening_factor)
Ny_DOM = int(Ny_inflow_outflow/DOM_coarsening_factor)
Nz_DOM = int(Nz_inflow_outflow/DOM_coarsening_factor)

# Round number of particles to the nearest number divisible by the number of partitions
Np_inflow_outflow     = int(np.ceil(Np_inflow_outflow     / N_partitions_x_inflow_outflow)) * N_partitions_x_inflow_outflow
Np_inflow_outflow_max = int(np.ceil(Np_inflow_outflow_max / N_partitions_x_inflow_outflow)) * N_partitions_x_inflow_outflow

# Check that the number of grid cells in the x direction is divisible by the number of partitions
if (Nx_inflow_outflow %  N_partitions_x_inflow_outflow) != 0:
  print('ERROR')
  print('Nx_inflow_outflow ({}) is not evenly divisible by N_partitions_x_inflow_outflow ({})'.format(Nx_inflow_outflow, N_partitions_x_inflow_outflow))
  sys.exit()

# Check that the number of grid cells in the x direction is divisible by the number of partitions
if (Nx_DOM % N_partitions_x_inflow_outflow) != 0:
  print('ERROR')
  print('Nx_DOM ({}) is not evenly divisible by N_partitions_x_inflow_outflow ({})'.format(Nx_inflow_outflow, N_partitions_x_inflow_outflow))
  sys.exit()


print('###############################################################################')
print('#                                     Summary                                 #')
print('###############################################################################')
print('###############################################################################')
print('# Grid Properties')
print('###############################################################################')
print('dx = {} [m]'.format(dx))
print('dy = {} [m]'.format(dy))
print('dz = {} [m]'.format(dz))
print('dx_dns = {} [m]'.format(dx_dns))
print('dy_dns = {} [m]'.format(dy_dns))
print('dz_dns = {} [m]'.format(dz_dns))

print('')
print('###############################################################################')
print('# Fluid Properties')
print('###############################################################################')
print('rho = {} [kg/m^3]'.format(rho_0))
print('mu  = {} [Pa*s]'.format(mu))
print('nu = {} [m^2/s]'.format(nu))
print('T   = {} [K]'.format(T_0))
print('c_p = {} [J/(kg*K)]'.format(c_p_fluid))
print('c_v = {} [J/(kg*K)]'.format(c_v_fluid))
print('gamma = {} []'.format(gamma))
print('specificGasConstant = {} [J/(kg*K)]'.format(specificGasConstant))
print('k_fluid = {} [W/(m^2*K)]'.format(k_fluid))
print('P_0 = {} [Pa]'.format(P_0))
print('Pr = {} []'.format(Pr))
print('')

print('')
print('###############################################################################')
print('# Flow Scales')
print('###############################################################################')
print('U_0 = {} [m/s]'.format(U_0))
print('u_rms = {} [m/s]'.format(u_rms))
print('Re_W = {}'.format(Re_W))
print('Re_L = {}'.format(Re_L))
print('Re_lambda = {}'.format(Re_lambda))
print('Taylor Length Scale: lambda = {}'.format(lambda_))
print('Mean TKE dissipation rate: epsilon = {} [m^2/s^3]'.format(epsilon))
print('Kolmogorov length scale: eta = {} [m]'.format(eta))
print('Kolmogorov time scale: : tau_eta = {} [s]'.format(tau_eta))
print('Kolmogorov velocity scale: : u_eta = {} [m/s]'.format(u_eta))
print('dns grid factor = {} * eta'.format(np.pi*np.sqrt(2)/3))
print('c = {} [m/s] speed of sound'.format(c))
print('delta_t_c  = {} [s]'.format(delta_t_c))

print('')
print('###############################################################################')
print('# Particle Properties')
print('###############################################################################')
print('Re_p = {}'.format(Re_p))
print('Particle relaxation time: tau_p  = {}'.format(tau_p))
print('Particle thermal relaxation time: tau_th  = {}'.format(tau_th))
print('W_spread = {}'.format(W_spread))
print('W_spread_percent = {}'.format(W_spread_percent))

print('')
print('###############################################################################')
print('# Time Scale')
print('###############################################################################')
print( 'Frequency copy between domains: \t delta_t = {}'.format(ratio_copy))
print('Ratio scales to c time step: \t ratio = {}'.format(ratio_delta_t))

print('')
print('###############################################################################')
print('The HIT Scales')
print('###############################################################################')
print('l_0  = {} [m]'.format(l_0))
print('u_0  = {} [m/s]'.format(u_0))
print('t_0  = {} [s]'.format(t_0))
print('')


######################
# print output for soleil x input file
######################
print('##########################################')
print('           HIT Generator Section')
print('##########################################')
print('###########')
print('   Grid')
print('###########')
print('xNum = {}'.format(Nx_hit))
print('yNum = {}'.format(Ny_hit))
print('zNum = {}'.format(Nz_hit))
print('N_cells = {}'.format(Nx_hit*Ny_hit*Nz_hit))
print('xNum_dns = {}'.format(Nx_hit_dns))
print('yNum_dns = {}'.format(Ny_hit_dns))
print('zNum_dns = {}'.format(Nz_hit_dns))
print('N_cells_dns = {}'.format(Nx_hit_dns*Ny_hit_dns*Nz_hit_dns))
print('xWidth = {}'.format(W_interior))
print('yWidth = {}'.format(W_interior))
print('zWidth = {}'.format(W_interior))
print('')

print('#############')
print('    Flow')
print('#############')
print('turbForcing')
print('  meanVelocity = [{}, 0.0, 0.0]'.format(U_0))
print('  G   = {}'.format(67))
print('  t_0 = {}'.format(t_0))
print('  k_0 = {}'.format(k_0))
print('')

print('#############')
print('  Particles')
print('#############')
print('initNum = {}'.format(Np_HIT))
print('maxNum = {}'.format(Np_HIT_max))

print('##########################################')
print('           Inflow Outflow Section')
print('##########################################')
print('###########')
print('   Grid')
print('###########')
print('xNum = {}'.format(Nx_inflow_outflow))
print('yNum = {}'.format(Ny_inflow_outflow))
print('zNum = {}'.format(Nz_inflow_outflow))
print('N_cells = {}'.format(Nx_inflow_outflow*Ny_inflow_outflow*Nz_inflow_outflow))
print('xNum_dns = {}'.format(Nx_inflow_outflow_dns))
print('yNum_dns = {}'.format(Ny_inflow_outflow_dns))
print('zNum_dns = {}'.format(Nz_inflow_outflow_dns))
print('N_cells_dns = {}'.format(Nx_inflow_outflow_dns*Ny_inflow_outflow_dns*Nz_inflow_outflow_dns))
print('xWidth = {}'.format(L))
print('yWidth = {}'.format(W_inflow_outflow))
print('zWidth = {}'.format(W_inflow_outflow))
print('')

print('#############')
print(' Integrator')
print('#############')
print('maxIter = {}'.format(Nt))
print('fixedDeltaTime = {}'.format(dt))
print('staggerFactor = {}'.format(ratio_delta_t))
print('')

print('#############')
print('    Flow')
print('#############')
print('gasConstant = {}'.format(specificGasConstant))
print('gamma = {}'.format(gamma))
print('prandtl = {}'.format(Pr))
print('mu = {}'.format(mu))
print('initParams = [{}, {}, {}, 0.0, 0.0, 0.0] = [rho, P, u, v, w, pertub_magnitude]'.format(rho_0, P_0, U_0))
print('bodyForce = [{}, 0.0, 0.0]'.format(g))
print('')

print('#############')
print('  Particles')
print('#############')
print('initNum = {}'.format(Np_inflow_outflow))
print('maxNum = {}'.format(int(1.10*Np_inflow_outflow)))
print('restitutionCoeff = {}'.format(restitution_coefficient))
print('convectiveCoeff = {}'.format(h))
print('heatCapacity = {}'.format(c_v_particle))
print('initTemperature = {}'.format(T_0))
print('density = {}'.format(rho_p))
print('diameterMean = {}'.format(d_p))
print('bodyForce = [{}, 0.0, 0.0]'.format(g))
print('maxSkew = {}'.format(maxSkew))
print('escapeRatioPerDir = {}'.format(escapeRatioPerDir))
print('staggerFactor= {}'.format(ratio_delta_t))
print('')

print('#############')
print('     IO')
print('#############')
print('restartEveryTimeSteps = {}'.format(Nt_between_restarts))
print('')


print('#############')
print('  Configs')
print('#############')
print('toCell = [0, {}, {}]'.format(Ny_hit,Nz_hit))
print('copyEveryTimeSteps= {}'.format(ratio_copy))
print('')
print('')
print('')


print('#############')
print('  Data Size')
print('#############')
N_fluid_cells = int(Nx_inflow_outflow*Ny_inflow_outflow*Nz_inflow_outflow)
N_DOM_cells   = int(Nx_DOM*Ny_DOM*Nz_DOM)


print('grid size = [{}, {}, {}]'.format(Nx_inflow_outflow,Ny_inflow_outflow,Nz_inflow_outflow))
print('number cells = {}'.format(N_fluid_cells))
print('number cells per tile = {}'.format(N_fluid_cells/N_partitions_x_inflow_outflow))
print('number cells per node = {}'.format(N_fluid_cells/N_nodes_inflow_outflow))

doubles_per_fluid_cell = 70
doubles_per_dom_cell   = 12
doubles_per_particle   = 35
fluid_gigabytes = (N_fluid_cells*doubles_per_fluid_cell*8)/(1.0*np.power(10,9))
dom_gigabytes = (N_DOM_cells*doubles_per_dom_cell*8)/(1.0*np.power(10,9))
particle_gigabytes = (Np_inflow_outflow_max*doubles_per_particle*8)/(1.0*np.power(10,9))
total_gigabytes = fluid_gigabytes + dom_gigabytes + particle_gigabytes
print('fluid gigabytes required = {}'.format(fluid_gigabytes))
print('dom gigabytes required = {}'.format(dom_gigabytes))
print('particle gigabytes required = {}'.format(particle_gigabytes))
print('total gigabytes required = {}'.format(total_gigabytes))

gigs_per_gpu = 16
N_gpu_needed = np.ceil(total_gigabytes/gigs_per_gpu)
N_gpu_per_node = 4
N_node_needed = np.ceil(N_gpu_needed/N_gpu_per_node)
print('gpus needed = {}'.format(N_gpu_needed))
print('nodes needed = {}'.format(N_node_needed))



json_template['configs'][1]['Mapping']['tiles'][0] = N_partitions_x_hit 
json_template['configs'][1]['Mapping']['tilesPerRank'][0] = N_partitions_x_per_node 
json_template['configs'][0]['Grid']['xNum'] = Nx_hit
json_template['configs'][0]['Grid']['yNum'] = Ny_hit
json_template['configs'][0]['Grid']['zNum'] = Nz_hit
json_template['configs'][0]['Grid']['xWidth'] = W_interior
json_template['configs'][0]['Grid']['yWidth'] = W_interior
json_template['configs'][0]['Grid']['zWidth'] = W_interior
json_template['configs'][0]['Integrator']['maxIter'] = Nt
json_template['configs'][0]['Integrator']['fixedDeltaTime'] = dt
json_template['configs'][0]['Flow']['gasConstant'] = specificGasConstant 
json_template['configs'][0]['Flow']['gamma'] = gamma
json_template['configs'][0]['Flow']['prandtl'] = Pr
json_template['configs'][0]['Flow']['viscosityModel']['viscosity'] = mu 
json_template['configs'][0]['Flow']['initParams'][0] = rho_0
json_template['configs'][0]['Flow']['initParams'][1] = P_0
json_template['configs'][0]['Flow']['initParams'][2] = U_0
json_template['configs'][0]['Flow']['initParams'][5] = u_rms
json_template['configs'][0]['Flow']['bodyForce'][0] = g
json_template['configs'][0]['Flow']['turbForcing']['meanVelocity'][0] = U_0
json_template['configs'][0]['Flow']['turbForcing']['t_o'] = t_0
json_template['configs'][0]['Flow']['turbForcing']['K_o'] = k_0
json_template['configs'][0]['Particles']['initNum'] = Np_HIT 
json_template['configs'][0]['Particles']['maxNum'] = Np_HIT_max
json_template['configs'][0]['Particles']['restitutionCoeff'] = restitution_coefficient
json_template['configs'][0]['Particles']['convectiveCoeff'] = h 
json_template['configs'][0]['Particles']['heatCapacity'] = c_v_particle
json_template['configs'][0]['Particles']['initTemperature'] = T_0
json_template['configs'][0]['Particles']['density'] = rho_p
json_template['configs'][0]['Particles']['diameterMean'] = d_p
json_template['configs'][0]['Particles']['bodyForce'][0] = g
json_template['configs'][0]['Particles']['maxSkew'] = maxSkew
json_template['configs'][0]['Particles']['escapeRatioPerDir'] = escapeRatioPerDir
json_template['configs'][0]['Particles']['staggerFactor'] = ratio_delta_t
json_template['configs'][0]['IO']['restartEveryTimeSteps'] = Nt_between_restarts

json_template['configs'][1]['Mapping']['tiles'][0] = N_partitions_x_inflow_outflow
json_template['configs'][1]['Mapping']['tilesPerRank'][0] = N_partitions_x_per_node 
json_template['configs'][1]['Grid']['xNum'] = Nx_inflow_outflow 
json_template['configs'][1]['Grid']['yNum'] = Ny_inflow_outflow 
json_template['configs'][1]['Grid']['zNum'] = Nz_inflow_outflow 
json_template['configs'][1]['Grid']['xWidth'] = L
json_template['configs'][1]['Grid']['yWidth'] = W_inflow_outflow
json_template['configs'][1]['Grid']['zWidth'] = W_inflow_outflow
json_template['configs'][1]['BC']['xBCLeft']['VelocityProfile']['initialVelocity'][0] = U_0
json_template['configs'][1]['BC']['xBCLeft']['VelocityProfile']['addedVelocity'] = 0.0
json_template['configs'][1]['BC']['xBCLeft']['TemperatureProfile']['initialTemperature'] = T_0 
json_template['configs'][1]['BC']['xBCRight']['P_inf'] = P_0
json_template['configs'][1]['BC']['yBCLeft']['Velocity'][0]  = U_0
json_template['configs'][1]['BC']['yBCRight']['Velocity'][0] = U_0
json_template['configs'][1]['BC']['zBCLeft']['Velocity'][0]  = U_0
json_template['configs'][1]['BC']['zBCRight']['Velocity'][0] = U_0
json_template['configs'][1]['Integrator']['maxIter'] = Nt
json_template['configs'][1]['Integrator']['fixedDeltaTime'] = dt
json_template['configs'][1]['Flow']['gasConstant'] = specificGasConstant 
json_template['configs'][1]['Flow']['gamma'] = gamma
json_template['configs'][1]['Flow']['prandtl'] = Pr
json_template['configs'][1]['Flow']['viscosityModel']['viscosity'] = mu 
json_template['configs'][1]['Flow']['initParams'][0] = rho_0
json_template['configs'][1]['Flow']['initParams'][1] = P_0
json_template['configs'][1]['Flow']['initParams'][2] = U_0
#json_template['configs'][1]['Flow']['initParams'][5] = u_rms
json_template['configs'][1]['Flow']['initParams'][5] = 0.0
json_template['configs'][1]['Flow']['bodyForce'][0] = g
#json_template['configs'][1]['Particles']['initNum'] = Np_inflow_outflow
json_template['configs'][1]['Particles']['initNum'] = int(0)
json_template['configs'][1]['Particles']['maxNum'] = Np_inflow_outflow_max
json_template['configs'][1]['Particles']['restitutionCoeff'] = restitution_coefficient
json_template['configs'][1]['Particles']['convectiveCoeff'] = h 
json_template['configs'][1]['Particles']['heatCapacity'] = c_v_particle
json_template['configs'][1]['Particles']['initTemperature'] = T_0
json_template['configs'][1]['Particles']['density'] = rho_p
json_template['configs'][1]['Particles']['diameterMean'] = d_p
json_template['configs'][1]['Particles']['bodyForce'][0] = g
json_template['configs'][1]['Particles']['maxSkew'] = maxSkew
json_template['configs'][1]['Particles']['escapeRatioPerDir'] = escapeRatioPerDir
json_template['configs'][1]['Particles']['staggerFactor'] = ratio_delta_t
json_template['configs'][1]['Radiation']['qa'] = Q_a 
json_template['configs'][1]['Radiation']['qs'] = Q_s
json_template['configs'][1]['Radiation']['xNum'] = Nx_DOM
json_template['configs'][1]['Radiation']['yNum'] = Ny_DOM
json_template['configs'][1]['Radiation']['zNum'] = Nz_DOM
json_template['configs'][1]['Radiation']['angles'] = N_angles 
json_template['configs'][1]['Radiation']['yLoIntensity'] = Intensity 
json_template['configs'][1]['IO']['restartEveryTimeSteps'] = Nt_between_restarts
json_template['configs'][1]['IO']['probes'][0]['fromCell'][0] = Nx_inflow_outflow
json_template['configs'][1]['IO']['probes'][0]['uptoCell'][0] = Nx_inflow_outflow
json_template['configs'][1]['IO']['probes'][0]['uptoCell'][1] = Ny_inflow_outflow 
json_template['configs'][1]['IO']['probes'][0]['uptoCell'][2] = Nz_inflow_outflow 

N_ghost_per_side = 1

json_template['copySrc']['uptoCell'][1] = Ny_hit-1
json_template['copySrc']['uptoCell'][2] = Nz_hit-1
json_template['copyTgt']['fromCell'][1] = N_ghost_per_side + Ny_buffer_per_side 
json_template['copyTgt']['fromCell'][2] = N_ghost_per_side + Ny_buffer_per_side 
json_template['copyTgt']['uptoCell'][1] = N_ghost_per_side + Ny_buffer_per_side + (Ny_interior - 1)
json_template['copyTgt']['uptoCell'][2] = N_ghost_per_side + Ny_buffer_per_side + (Nz_interior - 1)
json_template['copyEveryTimeSteps'] = ratio_copy

#json.dump(json_template, sys.stdout, indent=4)
f = open("test_output.json", "w")
json.dump(json_template, f, indent=4)
f.close()

