import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------- #
#                                   User Input                                #
# --------------------------------------------------------------------------- #
#cell_centers_to_plot  = ['x', 'y', 'z']
#cell_width_to_plot    = ['x', 'y', 'z']
#cell_centers_to_print = ['x', 'y', 'z']
#cell_width_to_print   = ['x', 'y', 'z']

cell_centers_to_plot  = ['z']
cell_width_to_plot    = ['z']
cell_centers_to_print = ['z']
cell_width_to_print   = ['z']

# --------------------------------------------------------------------------- #
#                            Read Command Line Input                          #
# --------------------------------------------------------------------------- #

parser = argparse.ArgumentParser()
parser.add_argument('hdf_file', 
                    help='fluid file to read')
args = parser.parse_args()

# --------------------------------------------------------------------------- #
#                               Read the Data                                 #
# --------------------------------------------------------------------------- #
print(args.hdf_file)
f = h5py.File(args.hdf_file, 'r')

# Get the data
# Get cell center values in each direction
x_slice_idx = 0
y_slice_idx = 0
z_slice_idx = 0
x_centers = f['centerCoordinates'][z_slice_idx,y_slice_idx,   :      ][:,0]
y_centers = f['centerCoordinates'][z_slice_idx,    :     ,x_slice_idx][:,1]
z_centers = f['centerCoordinates'][    :     ,y_slice_idx,x_slice_idx][:,2]
dx_values = f['cellWidth'][z_slice_idx,y_slice_idx,   :      ][:,0]
dy_values = f['cellWidth'][z_slice_idx,    :     ,x_slice_idx][:,1]
dz_values = f['cellWidth'][    :     ,y_slice_idx,x_slice_idx][:,2]

Nx = len(x_centers)
Ny = len(y_centers)
Nz = len(z_centers)

print('Nx = {}'.format(Nx))
print('Ny = {}'.format(Ny))
print('Nz = {}'.format(Nz))

# --------------------------------------------------------------------------- #
#                           Plot and Print the Data                           #
# --------------------------------------------------------------------------- #
# Dictionary with key of direction to value of cell center or cell width in that dirction
center_values     = {'x' : x_centers, 'y' : y_centers, 'z' : z_centers}
cell_width_values = {'x' : dx_values, 'y' : dy_values, 'z' : dz_values}

# Plot cell centers in the requested directions
for direction in cell_centers_to_plot:
  plt.figure()
  plt.plot(center_values[direction], np.zeros(len(center_values[direction])),'o')
  plt.xlabel('{}'.format(direction))
  plt.title('{} values'.format(direction))

# Plot cell widths in the requested  each directions
for direction in cell_width_to_plot:
  plt.figure()
  plt.plot(range(len(cell_width_values[direction])),cell_width_values[direction],'o')
  plt.xlabel('{}_idx'.format(direction))
  plt.ylabel('d{}'.format(direction))
  plt.title('d{} values'.format(direction))


# Print cell centers in the requested directions
for direction in cell_centers_to_print:
  print('{}_center = {}'.format(direction, repr(center_values[direction])))

# Print cell centers in the requested directions
for direction in cell_width_to_print:
  print('d{} = {}'.format(direction, repr(cell_width_values[direction])))

plt.show()
