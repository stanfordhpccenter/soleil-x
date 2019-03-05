import sys
import numpy as np
import matplotlib.pyplot as plt

# Reads in numpy compressed array files (.npz) files with the mcrt solution and
# and the analytical solution and plots them data.

# --------------------------------------------------------------------------- #
#                            Read Command Line Input                          #
# --------------------------------------------------------------------------- #

if len(sys.argv) == 2:
    # MCRT data to postprocess
    simulation_solution_filename = sys.argv[1]

    # Analytical solution to compare with
    analytical_solution_filename = '1D_analytical_solution.npz'

elif len(sys.argv) == 3:
    # MCRT data to postprocess
    simulation_solution_filename = sys.argv[1]  # MCRT data to postprocess

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
#                              Load the MCRT Results                          #
# --------------------------------------------------------------------------- #

data_simulation = np.load(simulation_solution_filename)
x = data_simulation['x']
y = data_simulation['y']
z = data_simulation['z']
G = data_simulation['G']
#q_x = data_simulation['q_x']
#q_y = data_simulation['q_y']
#q_z = data_simulation['q_z']

# -------------------------------------------------------------------------- #
#                            Load the Analytical Solution                    #
# -------------------------------------------------------------------------- #

data_analytical = np.load(analytical_solution_filename)
x_analytical = data_analytical['x_analytical']
G_analytical = data_analytical['G_analytical']
q_analytical = data_analytical['q_analytical']

# -------------------------------------------------------------------------- #
#                              Compute L2 Error                              #
# -------------------------------------------------------------------------- #

# Find size of mcrt simulation data
(Nx, Ny, Nz) = np.shape(G)

# Find step size in each direction... Assume uniform mesh... 
if Nx == 1:
  dx = 2*x[0,0,0]
else:
  dx = x[1,0,0] - x[0,0,0]

if Ny == 1:
  dy = 2.0*y[0,0,0] 
else:
  dy = y[0,1,0] - y[0,0,0]

if Nz == 1:
  dz = 2.0*z[0,0,0] 
else:
  dz = z[0,0,1] - z[0,0,0]

# Points where the L2 error is going to be computed
x_compare = x_analytical
#y_compare = np.linspace(dy,Ny*dy-dy,10) 
#z_compare = np.linspace(dy,Nz*dz-dz,10)
y_compare = np.array([Ny*dy/2.0])
z_compare = np.array([Nz*dz/2.0])

# MCRT data at comparison points
a_at_compare  = np.empty([len(x_compare), len(y_compare), len(z_compare)])
G_at_compare  = np.empty([len(x_compare), len(y_compare), len(z_compare)])
qx_at_compare = np.empty([len(x_compare), len(y_compare), len(z_compare)])

# The error at the comparison points
a_error  = np.empty([len(x_compare), len(y_compare), len(z_compare)])
G_error  = np.empty([len(x_compare), len(y_compare), len(z_compare)])
qx_error = np.empty([len(x_compare), len(y_compare), len(z_compare)])

# Find the same 
for x_idx in xrange(len(x_compare)):
  for y_idx in xrange(len(y_compare)):
    for z_idx in xrange(len(z_compare)):

        # Find what cell the compare point is in
        if x_compare[x_idx] >= Nx*dx:
            i = Nx-1
        else:
            i = int(np.floor(x_compare[x_idx]/dx))
        if y_compare[y_idx] >= Ny*dy:
            j = Ny-1
        else:
            j = int(np.floor(y_compare[y_idx]/dy))
        if z_compare[z_idx] >= Nz*dz:
            k = Nz-1
        else:
            k = int(np.floor(z_compare[z_idx]/dz))

        # Set the compare values from that cell 
        G_at_compare[x_idx,y_idx,z_idx]  = G[i,j,k]
        #qx_at_compare[x_idx,y_idx,z_idx] = q_x[i,k,k]

        # Take the difference in I and G at that x-value
        G_error[ x_idx,y_idx,z_idx] = G_at_compare[x_idx,y_idx,z_idx] - \
                                      G_analytical[x_idx]
        #qx_error[x_idx,y_idx,z_idx] = qx_at_compare[x_idx,y_idx,z_idx]- \
        #                              q_analytical[x_idx]


# --------------------------------------------------------------------------- #
#                                     Plot G                                  #
# --------------------------------------------------------------------------- #
plt.figure(1)
plt.plot(x_analytical, G_analytical , '-k', label=r'Analytical')
plt.plot(x_compare, G_at_compare[:,
                                 int(len(y_compare)/2.0),  \
                                 int(len(z_compare)/2.0)], \
                                 '-b', fillstyle='full', label=r'Simulation')
plt.xlabel(r'$x \ [m]$', fontsize = 20)
plt.ylabel(r'$G \ \left[ \frac{W}{m^2} \right]$', fontsize = 20)
plt.legend(loc = 'best')
plt.savefig('G_Nx_{}.pdf'.format(Nx), bbox_inches='tight')

# --------------------------------------------------------------------------- #
#                                     Plot q_x                                #
# --------------------------------------------------------------------------- #
#plt.figure(2)
#plt.plot(x_analytical, q_analytical , '-k', label=r'Analytical')
#plt.plot(x_compare, qx_at_compare[:, int(len(y_compare)/2), int(len(z_compare)/2)], \
#         '-b', fillstyle='full', label=r'Simulation')

#plt.xlabel(r'$x \ [m]$', fontsize = 20)
#plt.ylabel(r'$q \ \left[ \frac{W}{m^2} \right]$', fontsize = 20)
#plt.legend(loc = 'best')
#plt.savefig('q_Nx_{}.pdf'.format(Nx), bbox_inches='tight')

# --------------------------------------------------------------------------- #
#                                 Plot Error                                  # 
# --------------------------------------------------------------------------- #

plt.figure(3)
plt.plot(x_compare, G_error[:,
                            int(len(y_compare)/2.0),  \
                            int(len(z_compare)/2.0)], \
                            '-b', fillstyle='full')
plt.xlabel(r'$x \ [m]$', fontsize = 20)
plt.ylabel(r'$G - G_{analytical} \ \left[ \frac{W}{m^2} \right]$', fontsize = 20)
plt.savefig('G_Error_Nx_{}.pdf'.format(Nx), bbox_inches='tight')


# --------------------------------------------------------------------------- #
#                         Print Output to Terminal                            #
# --------------------------------------------------------------------------- #

G_L2_error  = np.sqrt(np.mean(np.square(G_error )))
#qx_L2_error = np.sqrt(np.mean(np.square(qx_error)))
print('Nx           : {}'.format(Nx))
print('G   L2 error : {}'.format(repr(G_L2_error)))
#print('q_x L2 error : {}'.format(repr(qx_L2_error)))

plt.show()
