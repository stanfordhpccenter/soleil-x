import numpy as np
import os
import json
from pprint import pprint
import matplotlib.pyplot as plt
import sys
import h5py

##############################################################################
#                                 User Input                                 #
##############################################################################
if len(sys.argv) == 2:
    # Data to postprocess
    filename = sys.argv[1]
else :
    print('*************************************')
    print('Error: Not enough command line arguments')
    print('Useage: Pick one of the following')
    print('$ python {} your_data.hdf'.format(sys.argv[0]))
    print('*************************************')
    sys.exit(2)


dir_name = os.path.join(os.environ['SOLEIL_DIR'], 'verification/couette')

soleil_input_file = os.path.join(dir_name, 'couette.json')
#hdf_filename = os.path.join(dir_name, 'sample0/fluid_iter0000008000/0,0,0-31,33,31.hdf')
hdf_filename = filename


##############################################################################
#                           Read Soleil Input File                           #
##############################################################################

with open(soleil_input_file) as f:
  data = json.load(f)

#pprint(data)

yNum   = data["Grid"]["yNum"]
yWidth = data["Grid"]["yWidth"]
U      = data["BC"]["yBCRightVel"][0]

# this solution assumes that y_min = 0.0
if not (data["Grid"]["origin"] == [0.0,0.0,0.0]):
  sys.exit() 

# this solution assumes that u_bottom wall = 0.0
if not (data["BC"]["yBCLeftVel"] == [0.0,0.0,0.0]):
  sys.exit() 

remove_y_ghost = True
if data["BC"]["yBCLeftVel"] == 'Periodic':
  remove_y_ghost = False


##############################################################################
#                          Read Soleil Output Data                           #
##############################################################################


f = h5py.File(hdf_filename, 'r')

# List all groups
print('Data Sets:')
for k in f.keys() :
  print(k)
print('')

# Get the data
centerCoordinates = f['centerCoordinates']
pressure    = f['pressure']
rho         = f['rho']
velocity    = f['velocity']
temperature = f['temperature']

if remove_y_ghost:
  Ny = rho.shape[1]
  centerCoordinates = centerCoordinates[:,1:Ny-1,:]
  pressure = pressure[:,1:Ny-1,:]
  rho = rho[:,1:Ny-1,:]
  velocity = velocity[:,1:Ny-1,:]
  temperature = temperature[:,1:Ny-1,:]

# Get dimension of data
Nx = rho.shape[2]
Ny = rho.shape[1]
Nz = rho.shape[0]

# Get simulation data along a line (ignore ghost cells)
x_slice_idx = 0
z_slice_idx = 0
#u_slice = velocity[z_slice_idx,:,x_slice_idx][:,0][1:Ny-1]
#rho_slice = rho[z_slice_idx,:,x_slice_idx][1:Ny-1]
#pressure_slice = pressure[z_slice_idx,:,x_slice_idx][1:Ny-1]
#temperature_slice = temperature[z_slice_idx,:,x_slice_idx][1:Ny-1]

y_slice   =   centerCoordinates[z_slice_idx,:,x_slice_idx][:,1]
u_slice   =            velocity[z_slice_idx,:,x_slice_idx][:,0]
rho_slice =                 rho[z_slice_idx,:,x_slice_idx]
pressure_slice =       pressure[z_slice_idx,:,x_slice_idx]
temperature_slice = temperature[z_slice_idx,:,x_slice_idx]

##############################################################################
#                          Compute Analytical Solution                       #
##############################################################################
# Couette profile
def u(y):
  return (U/yWidth)*y

u_slice_analytical = u(y_slice)

##############################################################################
#                          Read Soleil Output Data                           #
##############################################################################

L2_error = np.linalg.norm(u_slice-u_slice_analytical)
print('L2 Error = {}'.format(L2_error))

plt.figure(1)
plt.plot(u_slice, y_slice, 'ok', label='Soleil-X')
plt.plot(u_slice_analytical, y_slice, '-b', label='analytical')
plt.xlabel(r'$u \ \left[ \frac{m}{s} \right]$', fontsize = 20)
plt.ylabel(r'$y \ [m]$', fontsize = 20)
plt.legend()
plt.savefig('couette_solutions.pdf', bbox_inches='tight')

plt.show()
