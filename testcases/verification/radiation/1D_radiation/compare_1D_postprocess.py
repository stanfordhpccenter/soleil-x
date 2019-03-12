import sys
import numpy as np
import matplotlib.pyplot as plt

# Reads in numpy compressed array files (.npz) with the simulation and
# analytical solutions, plots them, and computes the L2 norm of the error.

# --------------------------------------------------------------------------- #
#                            Read Command Line Input                          #
# --------------------------------------------------------------------------- #

if len(sys.argv) == 2:
    # Simulations data to postprocess
    simulation_solution_filename = sys.argv[1]

    # Analytical solution to compare with
    analytical_solution_filename = '1D_analytical_solution.npz'

elif len(sys.argv) == 3:
    # Simulations data to postprocess
    simulation_solution_filename = sys.argv[1]

    # Analytical solution to compare with
    analytical_solution_filename = sys.argv[2]

else :
    print('*************************************')
    print('Error: Not enough command line arguments')
    print('Useage: Pick one of the following')
    print('$ python {} simulation_solution_filename'.format(sys.argv[0]))
    print('$ python {} simulation_solution_filenamean alytical_solution_filename'\
                                                       .format(sys.argv[0]))
    print('*************************************')
    sys.exit(2)

# --------------------------------------------------------------------------- #
#                         Load the Simulation Results                         #
# --------------------------------------------------------------------------- #

data_simulation = np.load(simulation_solution_filename)
# These are assumed to be 3D numpy arrays
x_simulation  = data_simulation['x']
y_simulation  = data_simulation['y']
z_simulation  = data_simulation['z']
dx_simulation = data_simulation['dx']
dy_simulation = data_simulation['dy']
dz_simulation = data_simulation['dz']
G_simulation  = data_simulation['G']

# -------------------------------------------------------------------------- #
#                            Load the Analytical Solution                    #
# -------------------------------------------------------------------------- #

data_analytical = np.load(analytical_solution_filename)
# These are assumed to be 1D numpy arrays
x_analytical = data_analytical['x_analytical']
G_analytical = data_analytical['G_analytical']
q_analytical = data_analytical['q_analytical']

# -------------------------------------------------------------------------- #
#                               Print Data Summary                           #
# -------------------------------------------------------------------------- #

print('Simulation Data:')
print('shape  : {}'.format(G_simulation.shape))
print('min dx : {}'.format(np.min(dx_simulation)))
print('max dx : {}'.format(np.max(dx_simulation)))
#print('min dy : {}'.format(np.min(dy_simulation)))
#print('max dy : {}'.format(np.max(dy_simulation)))
#print('min dz : {}'.format(np.min(dz_simulation)))
#print('max dz : {}'.format(np.max(dz_simulation)))
print('')
print('Simulation Data:')
print('shape  : {}'.format(x_analytical.shape))
print('dx : {}'.format(x_analytical[1]-x_analytical[0]))
print('')
if np.min(dx_simulation) < (x_analytical[1]-x_analytical[0]):
  print('WARNING: Simulation solution is at a higher resolution than the analytical one')
print('')

# -------------------------------------------------------------------------- #
#              Interpolate Analtyical Solution to Simulation Points          #
# -------------------------------------------------------------------------- #
# Find size of simulation data
(Nx, Ny, Nz) = np.shape(G_simulation)

# Find the analtyical solution at the simulation data points
G_analytical_at_simulation_points = np.empty([Nx, Ny, Nz])
for i, G_analytical_interp_at_simulation_x in enumerate(np.interp(x_simulation[:,0,0], x_analytical, G_analytical)):
  G_analytical_at_simulation_points[i,:,:] = G_analytical_interp_at_simulation_x 

# -------------------------------------------------------------------------- #
#                                Compute Error                               #
# -------------------------------------------------------------------------- #

# The error at each simulation point
G_error = G_simulation - G_analytical_at_simulation_points

# -------------------------------------------------------------------------- #
#                   Find Simulation Points at Analytical x                   #
# -------------------------------------------------------------------------- #
G_simulation_at_analytical_x = np.empty([len(x_analytical[1:-1])])
for i_ana, i_sim in enumerate(np.searchsorted(x_simulation[:,0,0]+0.5*dx_simulation[:,0,0],x_analytical[1:-1])):
    G_simulation_at_analytical_x[i_ana] = G_simulation[i_sim,int(Ny/2),int(Nz/2)]

# --------------------------------------------------------------------------- #
#                                     Plot G                                  #
# --------------------------------------------------------------------------- #
plt.figure()
plt.plot(x_analytical, G_analytical , '-k', label=r'Analytical')
plt.plot(x_simulation[:, int(Ny/2), int(Nz/2)], 
         G_simulation[:, int(Ny/2), int(Nz/2)], 
                      'ob', fillstyle='full', label=r'Simulation Cell Center')
plt.plot(x_analytical[1:-1], 
         G_simulation_at_analytical_x,
         '-b', fillstyle='full', label=r'Simulation Cell')
plt.xlabel(r'$x \ [m]$', fontsize = 20)
plt.ylabel(r'$G \ \left[ \frac{W}{m^2} \right]$', fontsize = 20)
plt.legend(loc = 'best')
plt.savefig('G_Nx_{}.pdf'.format(Nx), bbox_inches='tight')

# --------------------------------------------------------------------------- #
#                                 Plot Error                                  # 
# --------------------------------------------------------------------------- #

plt.figure()
plt.plot(x_simulation[:, int(Ny/2), int(Nz/2)], 
         G_error[:, int(Ny/2), int(Nz/2)], 
         '-ok', fillstyle='full')
plt.xlabel(r'$x \ [m]$', fontsize = 20)
plt.ylabel(r'$G - G_{analytical} \ \left[ \frac{W}{m^2} \right]$', fontsize = 20)
plt.savefig('G_Error_Nx_{}.pdf'.format(Nx), bbox_inches='tight')

# --------------------------------------------------------------------------- #
#                         Print Output to Terminal                            #
# --------------------------------------------------------------------------- #

G_L2_error  = np.sqrt(np.mean(np.square(G_error)))
G_L2_error_center_line  = np.sqrt(np.mean(np.square(G_error[:,int(Ny/2),int(Nz/2)])))
print('Nx : {}'.format(Nx))
print('G L2 error global      : {}'.format(repr(G_L2_error)))
print('G L2 error center line : {}'.format(repr(G_L2_error_center_line)))

plt.show()
