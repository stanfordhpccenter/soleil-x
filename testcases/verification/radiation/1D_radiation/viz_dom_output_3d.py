import numpy as np
import sys
from mayavi import mlab

# list of things to plot
scalar_data_to_plot = ['G']
vector_data_to_plot = []

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
    print('$ python {} your_data.npz'.format(sys.argv[0]))
    print('*************************************')
    sys.exit(2)

# --------------------------------------------------------------------------- #
#                               Load the Results                              #
# --------------------------------------------------------------------------- #

data = np.load(filename)
x = data['x']
y = data['y']
z = data['z']
G = data['G']

data_dict = {'x' : x,
             'y' : y, 
             'z' : z,
             'G' : G}

# --------------------------------------------------------------------------- #
#                            Pull out the meta data                           #
# --------------------------------------------------------------------------- #

# Find size of simulation data
(Nx, Ny, Nz) = np.shape(G)

print('Nx = {}'.format(Nx))
print('Ny = {}'.format(Ny))
print('Nz = {}'.format(Nz))

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

# Use this for axes limits
x_min = min(x[:,0,0]) - dx/2.0
y_min = min(y[0,:,0]) - dy/2.0
z_min = min(z[0,0,:]) - dz/2.0
x_max = max(x[:,0,0]) + dx/2.0
y_max = max(y[0,:,0]) + dy/2.0
z_max = max(z[0,0,:]) + dz/2.0

# --------------------------------------------------------------------------- #
#                              Plot Slices of Results                         #
# --------------------------------------------------------------------------- #


for scalar_feild_name in scalar_data_to_plot:

  figure = mlab.figure('{} in Plane'.format(scalar_feild_name))
  mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(np.array(data_dict['{}'.format(scalar_feild_name)])),
                                   plane_orientation='x_axes',
                                   slice_index=0)
  mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(np.array(data_dict['{}'.format(scalar_feild_name)])),
                                   plane_orientation='y_axes',
                                   slice_index=0)
  mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(np.array(data_dict['{}'.format(scalar_feild_name)])),
                                   plane_orientation='z_axes',
                                   slice_index=0)

  mlab.outline()
  mlab.axes(ranges=[z_min, z_max, y_min, y_max, x_min, x_max])
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
    mlab.axes(ranges=[z_min, z_max, y_min, y_max, x_min, x_max])
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
