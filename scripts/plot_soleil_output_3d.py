import numpy as np
import sys
import h5py
from mayavi import mlab

# list of things to plot
scalar_data_to_plot = ['rho','pressure','temperature']
vector_data_to_plot = ['velocity']
#vector_data_to_plot = ['cellWidth','velocity']
#vector_data_to_plot = ['cellWidth','velocity','velocityGradientX']

#vector_data_to_plot = ['cellWidth','velocity','velocityGradientX','velocityGradientY','velocityGradientZ','temperatureGradient']

#scalar_data_to_plot = ['rho','pressure','temperature','debug_scalar']
#vector_data_to_plot = ['velocity','velocityGradientY','debug_vector1', 'debug_vector2', 'debug_vector3']

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

# Note that the hdf5 files are dumped in [k,j,i] or [z,y,x] format
Nx = f['rho'].shape[2]
Ny = f['rho'].shape[1]
Nz = f['rho'].shape[0]

print('Nx = {}'.format(Nx))
print('Ny = {}'.format(Ny))
print('Nz = {}'.format(Nz))

## Get cell center values in each direction
#x_slice_idx = 0
#y_slice_idx = 0
#z_slice_idx = 0
#x_values = f['centerCoordinates'][z_slice_idx,y_slice_idx,   :      ][:,0]
#y_values = f['centerCoordinates'][z_slice_idx,    :     ,x_slice_idx][:,1]
#z_values = f['centerCoordinates'][    :     ,y_slice_idx,x_slice_idx][:,2]
#
#x_min = min(x_values) - 0.5*f['cellWidth'][0,0,0][0]
#y_min = min(y_values) - 0.5*f['cellWidth'][0,0,0][1]
#z_min = min(z_values) - 0.5*f['cellWidth'][0,0,0][2]
#x_max = max(x_values) + 0.5*f['cellWidth'][Nz-1,Ny-1,Nx-1][0]
#y_max = max(y_values) + 0.5*f['cellWidth'][Nz-1,Ny-1,Nx-1][1]
#z_max = max(z_values) + 0.5*f['cellWidth'][Nz-1,Ny-1,Nx-1][2]

# --------------------------------------------------------------------------- #
#                              Plot Slices of Results                         #
# --------------------------------------------------------------------------- #


for scalar_feild_name in scalar_data_to_plot:

  figure = mlab.figure('{} in Plane'.format(scalar_feild_name))
  mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(np.array(f['{}'.format(scalar_feild_name)])),
                                   plane_orientation='x_axes',
                                   slice_index=0)
  mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(np.array(f['{}'.format(scalar_feild_name)])),
                                   plane_orientation='y_axes',
                                   slice_index=0)
  mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(np.array(f['{}'.format(scalar_feild_name)])),
                                   plane_orientation='z_axes',
                                   slice_index=0)

  mlab.outline()
#  mlab.axes(ranges=[z_min, z_max, y_min, y_max, x_min, x_max])
  mlab.xlabel('z [m]')
  mlab.ylabel('y [m]')
  mlab.zlabel('x [m]')
  mlab.title('{}'.format(scalar_feild_name))
  mlab.colorbar(orientation='vertical')

for vector_feild_name in vector_data_to_plot:
  vector = np.array(f['{}'.format(vector_feild_name)][:,:,:])
  for vector_component in ['x','y','z']:
    idx = {'x':0, 'y':1, 'z':2}

    figure = mlab.figure('{}_{} in Plane'.format(vector_feild_name,vector_component))
    mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(vector[:,:,:,idx[vector_component]]),
                                     plane_orientation='x_axes',
                                     slice_index=0)
    mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(vector[:,:,:,idx[vector_component]]),
                                     plane_orientation='y_axes',
                                     slice_index=0)
    mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(vector[:,:,:,idx[vector_component]]),
                                     plane_orientation='z_axes',
                                     slice_index=0)

    mlab.outline()
#    mlab.axes(ranges=[z_min, z_max, y_min, y_max, x_min, x_max])
    mlab.xlabel('z [m]')
    mlab.ylabel('y [m]')
    mlab.zlabel('x [m]')
    mlab.title('{}_{}'.format(vector_feild_name,vector_component))
    mlab.colorbar(orientation='vertical')


mlab.show()

################################################################################
#figure = mlab.figure('velocity in Plane')
#velocity = np.array(f['{}'.format('velocity')][:,:,:])
#src = mlab.pipeline.vector_field(velocity[0], velocity[1], velocity[2])
#mlab.pipeline.vector_cut_plane(src, scale_factor=3)
##mlab.pipeline.vector_cut_plane(src, mask_points=1000000, scale_factor=5)
##mlab.axes(ranges=[x_min, x_max, y_min, y_max, z_min, z_max])
##mlab.axes(ranges=[x_min, x_max, y_min, y_max, z_min, z_max], extent=[1, mesh.Nx, 1, mesh.Ny, 1, mesh.Nz])
#mlab.xlabel('x [m]')
#mlab.ylabel('y [m]')
#mlab.zlabel('z [m]')
#mlab.title('velocity')
#mlab.vectorbar(orientation='vertical')
################################################################################
