import sys
import numpy as np
import h5py
#from mayavi import mlab
import matplotlib.pyplot as plt

# list of things to plot
#scalar_data_to_plot = ['rho','pressure','temperature']
#vector_data_to_plot = ['cellWidth','velocity','velocityGradientX','velocityGradientY','velocityGradientZ','temperatureGradient']

scalar_data_to_plot = ['rho','pressure','temperature', 'debug_scalar']
vector_data_to_plot = ['velocity','velocityGradientY', 'debug_vector1', 'debug_vector2', 'debug_vector3']

# --------------------------------------------------------------------------- #
#                            Read Command Line Input                          #
# --------------------------------------------------------------------------- #

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

# --------------------------------------------------------------------------- #

f = h5py.File(filename, 'r')

# List all groups
print('Data Sets:')
for k in f.keys() :
  print(k)
print('')

# Get the data
rho_hdf = f['rho']

# Note that the hdf5 files are dumped in [k,j,i] or [z,y,x] format
Nx = f['rho'].shape[2]
Ny = f['rho'].shape[1]
Nz = f['rho'].shape[0]

print('Nx = {}'.format(Nx))
print('Ny = {}'.format(Ny))
print('Nz = {}'.format(Nz))

x_slice_idx = 0
y_slice_idx = 0
z_slice_idx = 0
x_values = f['centerCoordinates'][z_slice_idx,y_slice_idx,   :      ][:,0]
y_values = f['centerCoordinates'][z_slice_idx,    :     ,x_slice_idx][:,1]
z_values = f['centerCoordinates'][    :     ,y_slice_idx,x_slice_idx][:,2]

## Plot cell centers in each direction
#plt.figure()
#plt.plot(x_values,np.zeros(len(x_values)),'o')
#plt.xlabel('x')
#plt.title('x values')
#
#plt.figure()
#plt.plot(y_values,np.zeros(len(y_values)),'o')
#plt.xlabel('y')
#plt.title('y values')
#
#plt.figure()
#plt.plot(z_values,np.zeros(len(z_values)),'o')
#plt.xlabel('y')
#plt.title('z values')

for scalar_feild_name in scalar_data_to_plot: 
  #scalar_feild_name = 'rho'
  plt.figure()
  plt.plot(f['{}'.format(scalar_feild_name)][z_slice_idx,:,x_slice_idx], y_values, 'ok', label='Soleil-X')
  plt.xlabel('{}'.format(scalar_feild_name), fontsize = 20)
  plt.ylabel(r'$y \ [m]$', fontsize = 20)
  #plt.legend()
  #plt.savefig('{}.pdf'.format(scalar_feild_name), bbox_inches='tight')

for vector_feild_name in vector_data_to_plot: 
  for vector_component in ['x','y','z']:
    idx = {'x':0, 'y':1, 'z':2}
    plt.figure()
    plt.plot(f['{}'.format(vector_feild_name)][z_slice_idx,:,x_slice_idx][:,idx[vector_component]], y_values, 'ok', label='Soleil-X')

    x_label = ''
    if vector_feild_name == 'velocityGradientX':
      if vector_component == 'x':
        x_label = r'$\frac{{\partial {}}}{{ \partial {} }}$'.format('u','x')
      elif vector_component == 'y':
        x_label = r'$\frac{{\partial {}}}{{ \partial {} }}$'.format('v','x')
      elif vector_component == 'z':
        x_label = r'$\frac{{\partial {}}}{{ \partial {} }}$'.format('w','x')
    elif vector_feild_name == 'velocityGradientY':
      if vector_component == 'x':
        x_label = r'$\frac{{\partial {}}}{{ \partial {} }}$'.format('u','y')
      elif vector_component == 'y':
        x_label = r'$\frac{{\partial {}}}{{ \partial {} }}$'.format('v','y')
      elif vector_component == 'z':
        x_label = r'$\frac{{\partial {}}}{{ \partial {} }}$'.format('w','y')
    elif vector_feild_name == 'velocityGradientZ':
      if vector_component == 'x':
        x_label = r'$\frac{{\partial {}}}{{ \partial {} }}$'.format('u','Z')
      elif vector_component == 'y':
        x_label = r'$\frac{{\partial {}}}{{ \partial {} }}$'.format('v','Z')
      elif vector_component == 'z':
        x_label = r'$\frac{{\partial {}}}{{ \partial {} }}$'.format('w','Z')
    else:
        x_label = '{}_{}'.format(vector_feild_name, vector_component)

    plt.title(x_label)
    plt.xlabel(x_label, fontsize = 20)
    plt.ylabel(r'$y \ [m]$', fontsize = 20)
    #plt.legend()
    #plt.savefig('{}.pdf'.format(vector_feild_name), bbox_inches='tight')

plt.show()

###############################################################################
#for feild_name in data_to_plot:
#
#  feild = f[feild_name]
#  # check to see if the feild is a vector or a scalar
#
#  figure = mlab.figure('{} in Plane'.format(feild_name))
#  mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(data_dict[feild_name]),
#                                   plane_orientation='x_axes',
#                                   slice_index=0)
#  mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(data_dict[feild_name]),
#                                   plane_orientation='y_axes',
#                                   slice_index=0)
#  mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(data_dict[feild_name]),
#                                   plane_orientation='z_axes',
#                                   slice_index=0)
#
#  mlab.outline()
#  mlab.axes(ranges=[x_min, x_max, y_min, y_max, z_min, z_max])
#  mlab.xlabel('x [m]')
#  mlab.ylabel('y [m]')
#  mlab.zlabel('z [m]')
#  mlab.title('{}'.format(feild_name))
#  mlab.colorbar(orientation='vertical')

#mlab.show()
###############################################################################
