import argparse
import numpy as np
import json
from pprint import pprint
import matplotlib.pyplot as plt
import sys
import os
import h5py

##############################################################################
#                                 User Input                                 #
##############################################################################
parser = argparse.ArgumentParser()
parser.add_argument('case_json_file', 
                    type=str,
                    help='soleil-x input json file')
parser.add_argument('hdf_file', 
                    type=str,
                    help='fluid restart file')
parser.add_argument('-v', '--verbose', 
                    action='store_true',
                    help='run in verbose mode')
args = parser.parse_args()

#dir_name = os.path.join(os.environ['SOLEIL_DIR'], 'testcases/verification/fluid/poiseuille')
#soleil_input_file = os.path.join(dir_name, 'poiseuille.json')

debug = True

##############################################################################
#                           Read Soleil Input File                           #
##############################################################################

#with open(soleil_input_file) as f:
with open(args.case_json_file) as f:
  data = json.load(f)

yNum = data["Grid"]["yNum"]
yWidth  = data["Grid"]["yWidth"]
g = data["Flow"]["bodyForce"][0]
constantVisc = data["Flow"]["viscosityModel"]["viscosity"]

# this solution assumes that y_min = 0.0
if not (data["Grid"]["origin"] == [0.0,0.0,0.0]):
  sys.exit() 


##############################################################################
#                          Read Soleil Output Data                           #
##############################################################################

f = h5py.File(args.hdf_file, 'r')

if args.verbose:
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

remove_y_ghost = True

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

y_slice =     centerCoordinates[z_slice_idx,:,x_slice_idx][:,1]
u_slice    =           velocity[z_slice_idx,:,x_slice_idx][:,0]
rho_slice  =                rho[z_slice_idx,:,x_slice_idx]
pressure_slice =       pressure[z_slice_idx,:,x_slice_idx]
temperature_slice = temperature[z_slice_idx,:,x_slice_idx]

##############################################################################
#                          Compute Analytical Solution                       #
##############################################################################

def u(y):
  return -(g/constantVisc)*np.power(y,2.0)/2.0 + (g/constantVisc)*(yWidth/2.0) * y

u_slice_analytical = u(y_slice)

##############################################################################
#                         Plot and Compare the Solutions                     #
##############################################################################

L2_error = np.linalg.norm(u_slice-u_slice_analytical)
print('L2 Error = {}'.format(L2_error))

plt.figure()
plt.plot(u_slice, y_slice, 'ok', label='Soleil-X')
plt.plot(u_slice_analytical, y_slice, '-b', label='analytical')
plt.xlabel(r'$u \ \left[ \frac{m}{s} \right]$', fontsize = 20)
plt.ylabel(r'$y \ [m]$', fontsize = 20)
plt.legend()
plt.savefig('poiseuille_solutions.pdf', bbox_inches='tight')

if debug:
  print('Temperature Range: [{}, {}]'.format(np.min(temperature), np.max(temperature)))
  print('Pressure Range: [{}, {}]'.format(np.min(pressure), np.max(pressure)))
  print('Density Range: [{}, {}]'.format(np.min(rho), np.max(rho)))
  plt.figure()
  plt.plot(temperature_slice, y_slice, 'ok', label='Soleil-X')
  plt.xlabel(r'$T \ \left[ K \right]$', fontsize = 20)
  plt.ylabel(r'$y \ [m]$', fontsize = 20)
  plt.legend()
  
  plt.figure()
  plt.plot(rho_slice, y_slice, 'ok', label='Soleil-X')
  plt.xlabel(r'$\rho \ \left[ \frac{kg}{m^3} \right]$', fontsize = 20)
  plt.ylabel(r'$y \ [m]$', fontsize = 20)
  plt.legend()
  
  plt.figure()
  plt.plot(pressure_slice, y_slice, 'ok', label='Soleil-X')
  plt.xlabel(r'$P \ \left[ Pa \right]$', fontsize = 20)
  plt.ylabel(r'$y \ [m]$', fontsize = 20)
  plt.legend()

 

plt.show()
