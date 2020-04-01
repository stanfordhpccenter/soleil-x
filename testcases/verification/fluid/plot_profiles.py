import argparse
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
parser = argparse.ArgumentParser()
#parser.add_argument('case_json_file', 
#                    type=str,
#                    help='soleil-x input json file')
parser.add_argument('hdf_file', 
                    type=str,
                    help='fluid restart file')
#parser.add_argument('wall_normal', 
#                    type=str,
#                    help='wall normal')
#parser.add_argument('shear_direction', 
#                    type=str,
#                    help='wall normal')
parser.add_argument('-v', '--verbose', 
                    action='store_true',
                    help='run in verbose mode')
args = parser.parse_args()

#dir_name = os.path.join(os.environ['SOLEIL_DIR'], 'testcases/verification/fluid/couette')
#soleil_input_file = os.path.join(dir_name, 'couette.json')

debug = False

##############################################################################
#                           Read Soleil Input File                           #
##############################################################################

##with open(soleil_input_file) as f:
#with open(args.case_json_file) as f:
#  data = json.load(f)
#
##pprint(data)
#
#yNum   = data["Grid"]["yNum"]
#yWidth = data["Grid"]["yWidth"]
#U      = data["BC"]["yBCRight"]["Velocity"][0]
#
## this solution assumes that y_min = 0.0
#if not (data["Grid"]["origin"] == [0.0,0.0,0.0]):
#  sys.exit() 
#
## this solution assumes that u_bottom wall = 0.0
#if not (data["BC"]["yBCLeft"]["Velocity"] == [0.0,0.0,0.0]):
#  sys.exit() 
#
#remove_y_ghost = True
#if data["BC"]["yBCLeft"] == 'Periodic':
#  remove_y_ghost = False

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

#if remove_y_ghost:
#  Ny = rho.shape[1]
#  centerCoordinates = centerCoordinates[:,1:Ny-1,:]
#  pressure = pressure[:,1:Ny-1,:]
#  rho = rho[:,1:Ny-1,:]
#  velocity = velocity[:,1:Ny-1,:]
#  temperature = temperature[:,1:Ny-1,:]

# Get dimension of data
Nx = rho.shape[2]
Ny = rho.shape[1]
Nz = rho.shape[0]

# Get simulation data along a line (ignore ghost cells)
x_slice_idx = int(Nx/2.0)
y_slice_idx = int(Ny/2.0)
z_slice_idx = int(Nz/2.0)

x_slice   =   centerCoordinates[z_slice_idx,y_slice_idx,:][:,0]
y_slice   =   centerCoordinates[z_slice_idx,:,x_slice_idx][:,1]
z_slice   =   centerCoordinates[:,y_slice_idx,x_slice_idx][:,2]
u_slice_x =            velocity[z_slice_idx,y_slice_idx,:][:,0]
v_slice_x =            velocity[z_slice_idx,y_slice_idx,:][:,1]
w_slice_x =            velocity[z_slice_idx,y_slice_idx,:][:,2]
u_slice_y =            velocity[z_slice_idx,:,x_slice_idx][:,0]
v_slice_y =            velocity[z_slice_idx,:,x_slice_idx][:,1]
w_slice_y =            velocity[z_slice_idx,:,x_slice_idx][:,2]
u_slice_z =            velocity[:,y_slice_idx,x_slice_idx][:,0]
v_slice_z =            velocity[:,y_slice_idx,x_slice_idx][:,1]
w_slice_z =            velocity[:,y_slice_idx,x_slice_idx][:,2]
rho_slice_x =               rho[z_slice_idx,y_slice_idx,:]
rho_slice_y =               rho[z_slice_idx,:,x_slice_idx]
rho_slice_z =               rho[:,y_slice_idx,x_slice_idx]
pressure_slice_x =     pressure[z_slice_idx,y_slice_idx,:]
pressure_slice_y =     pressure[z_slice_idx,:,x_slice_idx]
pressure_slice_z =     pressure[:,y_slice_idx,x_slice_idx]
temperature_slice_x = temperature[z_slice_idx,y_slice_idx,:]
temperature_slice_y = temperature[z_slice_idx,:,x_slice_idx]
temperature_slice_z = temperature[:,y_slice_idx,x_slice_idx]

##############################################################################
#                          Compute Analytical Solution                       #
##############################################################################
## Couette profile
#def u(y):
#  return (U/yWidth)*y
#
#u_slice_analytical = u(y_slice)

##############################################################################
#                         Plot and Compare the Solutions                     #
##############################################################################

##################
## all 3 velocity components vs slice plots
##################
#plt.figure()
#plt.plot(u_slice_x, x_slice, '-ok', label='u Soleil-X')
#plt.plot(v_slice_x, x_slice, '-xr', label='v Soleil-X')
#plt.plot(w_slice_x, x_slice, '-sb', label='w Soleil-X')
#plt.xlabel(r'$u \ \left[ \frac{m}{s} \right]$', fontsize = 20)
#plt.ylabel(r'$x \ [m]$', fontsize = 20)
#plt.legend()
#
#plt.figure()
#plt.plot(u_slice_y, y_slice, '-ok', label='u Soleil-X')
#plt.plot(v_slice_y, y_slice, '-xr', label='v Soleil-X')
#plt.plot(w_slice_y, y_slice, '-sb', label='w Soleil-X')
#plt.xlabel(r'$u \ \left[ \frac{m}{s} \right]$', fontsize = 20)
#plt.ylabel(r'$y \ [m]$', fontsize = 20)
#plt.legend()
#
#plt.figure()
#plt.plot(u_slice_z, z_slice, '-ok', label='u Soleil-X')
#plt.plot(v_slice_z, z_slice, '-xr', label='v Soleil-X')
#plt.plot(w_slice_z, z_slice, '-sb', label='w Soleil-X')
#plt.xlabel(r'$u \ \left[ \frac{m}{s} \right]$', fontsize = 20)
#plt.ylabel(r'$z \ [m]$', fontsize = 20)
#plt.legend()

#################
# x slice plots
#################
# u 
plt.figure()
plt.plot(u_slice_x, x_slice, '-ok', label='u Soleil-X')
plt.xlabel(r'$u \ \left[ \frac{m}{s} \right]$', fontsize = 20)
plt.ylabel(r'$x \ [m]$', fontsize = 20)
# v 
plt.figure()
plt.plot(v_slice_x, x_slice, '-ok', label='u Soleil-X')
plt.xlabel(r'$v \ \left[ \frac{m}{s} \right]$', fontsize = 20)
plt.ylabel(r'$x \ [m]$', fontsize = 20)
# w
plt.figure()
plt.plot(w_slice_x, x_slice, '-ok', label='u Soleil-X')
plt.xlabel(r'$w \ \left[ \frac{m}{s} \right]$', fontsize = 20)
plt.ylabel(r'$x \ [m]$', fontsize = 20)
# rho
plt.figure()
plt.plot(rho_slice_x, x_slice, '-ok', label='u Soleil-X')
plt.xlabel(r'$\rho \ \left[ \frac{kg}{m^3} \right]$', fontsize = 20)
plt.ylabel(r'$x \ [m]$', fontsize = 20)
# pressure
plt.figure()
plt.plot(pressure_slice_x, x_slice, '-ok', label='u Soleil-X')
plt.ylabel(r'$x \ [m]$', fontsize = 20)
plt.xlabel(r'$P \ \left[ Pa \right]$', fontsize = 20)
# temperature
plt.figure()
plt.plot(temperature_slice_x, x_slice, '-ok', label='u Soleil-X')
plt.ylabel(r'$x \ [m]$', fontsize = 20)
plt.xlabel(r'$T \ \left[ K \right]$', fontsize = 20)

#################
# y slice plots
#################
# u 
plt.figure()
plt.plot(u_slice_y, y_slice, '-ok', label='u Soleil-X')
plt.xlabel(r'$u \ \left[ \frac{m}{s} \right]$', fontsize = 20)
plt.ylabel(r'$y \ [m]$', fontsize = 20)
# v 
plt.figure()
plt.plot(v_slice_y, y_slice, '-ok', label='u Soleil-X')
plt.xlabel(r'$v \ \left[ \frac{m}{s} \right]$', fontsize = 20)
plt.ylabel(r'$y \ [m]$', fontsize = 20)
# w
plt.figure()
plt.plot(w_slice_y, y_slice, '-ok', label='u Soleil-X')
plt.xlabel(r'$w \ \left[ \frac{m}{s} \right]$', fontsize = 20)
plt.ylabel(r'$y \ [m]$', fontsize = 20)
# rho
plt.figure()
plt.plot(rho_slice_y, y_slice, '-ok', label='u Soleil-X')
plt.xlabel(r'$\rho \ \left[ \frac{kg}{m^3} \right]$', fontsize = 20)
plt.ylabel(r'$y \ [m]$', fontsize = 20)
# pressure
plt.figure()
plt.plot(pressure_slice_y, y_slice, '-ok', label='u Soleil-X')
plt.ylabel(r'$y \ [m]$', fontsize = 20)
plt.xlabel(r'$P \ \left[ Pa \right]$', fontsize = 20)
# temperature
plt.figure()
plt.plot(temperature_slice_y, y_slice, '-ok', label='u Soleil-X')
plt.ylabel(r'$y \ [m]$', fontsize = 20)
plt.xlabel(r'$T \ \left[ K \right]$', fontsize = 20)


#################
# z slice plots
#################
# u 
plt.figure()
plt.plot(u_slice_z, z_slice, '-ok', label='u Soleil-X')
plt.xlabel(r'$u \ \left[ \frac{m}{s} \right]$', fontsize = 20)
plt.ylabel(r'$z \ [m]$', fontsize = 20)
# v 
plt.figure()
plt.plot(v_slice_z, z_slice, '-ok', label='u Soleil-X')
plt.xlabel(r'$v \ \left[ \frac{m}{s} \right]$', fontsize = 20)
plt.ylabel(r'$z \ [m]$', fontsize = 20)
# w
plt.figure()
plt.plot(w_slice_z, z_slice, '-ok', label='u Soleil-X')
plt.xlabel(r'$w \ \left[ \frac{m}{s} \right]$', fontsize = 20)
plt.ylabel(r'$z \ [m]$', fontsize = 20)
# rho
plt.figure()
plt.plot(rho_slice_z, z_slice, '-ok', label='u Soleil-X')
plt.xlabel(r'$\rho \ \left[ \frac{kg}{m^3} \right]$', fontsize = 20)
plt.ylabel(r'$z \ [m]$', fontsize = 20)
# pressure
plt.figure()
plt.plot(pressure_slice_z, z_slice, '-ok', label='u Soleil-X')
plt.ylabel(r'$z \ [m]$', fontsize = 20)
plt.xlabel(r'$P \ \left[ Pa \right]$', fontsize = 20)
# temperature
plt.figure()
plt.plot(temperature_slice_z, z_slice, '-ok', label='u Soleil-X')
plt.ylabel(r'$z \ [m]$', fontsize = 20)
plt.xlabel(r'$T \ \left[ K \right]$', fontsize = 20)


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
