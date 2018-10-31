import numpy as np
import sys
import h5py
from mayavi import mlab

# list of things to plot
#scalar_data_to_plot = ['rho','pressure','temperature']
#vector_data_to_plot = ['cellWidth','velocity','velocityGradientX','velocityGradientY','velocityGradientZ','temperatureGradient']

scalar_data_to_plot = ['rho','pressure','temperature','debug_scalar']
vector_data_to_plot = ['velocity','velocityGradientY','debug_vector1', 'debug_vector2', 'debug_vector3']

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

x_slice_idx = 0
y_slice_idx = 0
z_slice_idx = 0
x_values = f['centerCoordinates'][z_slice_idx,y_slice_idx,   :      ][:,0]
y_values = f['centerCoordinates'][z_slice_idx,    :     ,x_slice_idx][:,1]
z_values = f['centerCoordinates'][    :     ,y_slice_idx,x_slice_idx][:,2]

x_min = min(x_values) - 0.5*f['cellWidth'][0,0,0][0]
y_min = min(y_values) - 0.5*f['cellWidth'][0,0,0][1]
z_min = min(z_values) - 0.5*f['cellWidth'][0,0,0][2]
x_max = max(x_values) + 0.5*f['cellWidth'][Nz-1,Ny-1,Nx-1][0]
y_max = max(y_values) + 0.5*f['cellWidth'][Nz-1,Ny-1,Nx-1][1]
z_max = max(z_values) + 0.5*f['cellWidth'][Nz-1,Ny-1,Nx-1][2]

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

#figure = mlab.figure('Pressure in Plane')
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(P),
#                                 plane_orientation='x_axes',
#                                 slice_index=0,
#                                )
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(P),
#                                 plane_orientation='y_axes',
#                                 slice_index=0,
#                                 )
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(P),
#                                 plane_orientation='z_axes',
#                                 slice_index=0,
#                                 )
#mlab.outline()
##mlab.axes(ranges=[x_min, x_max, y_min, y_max, z_min, z_max])
#mlab.xlabel('x [m]')
#mlab.ylabel('y [m]')
#mlab.zlabel('z [m]')
#mlab.title('Pressure')
#mlab.colorbar(orientation='vertical')
#
#
################################################################################
#figure = mlab.figure('Density in Plane')
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rho),
#                                 plane_orientation='x_axes',
#                                 slice_index=0,
#                                )
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rho),
#                                 plane_orientation='y_axes',
#                                 slice_index=0,
#                                 )
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rho),
#                                 plane_orientation='z_axes',
#                                 slice_index=0,
#                                 )
#mlab.outline()
##mlab.axes(ranges=[x_min, x_max, y_min, y_max, z_min, z_max])
#mlab.xlabel('x [m]')
#mlab.ylabel('y [m]')
#mlab.zlabel('z [m]')
#mlab.title('Density')
#mlab.colorbar(orientation='vertical')
#
#
################################################################################
#figure = mlab.figure('Temperature in Plane')
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(T),
#                                 plane_orientation='x_axes',
#                                 slice_index=0,
#                                )
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(T),
#                                 plane_orientation='y_axes',
#                                 slice_index=0,
#                                 )
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(T),
#                                 plane_orientation='z_axes',
#                                 slice_index=0,
#                                 )
#mlab.outline()
##mlab.axes(ranges=[x_min, x_max, y_min, y_max, z_min, z_max])
#mlab.xlabel('x [m]')
#mlab.ylabel('y [m]')
#mlab.zlabel('z [m]')
#mlab.title('Temperature')
#mlab.colorbar(orientation='vertical')
#
################################################################################
#figure = mlab.figure('u in Plane')
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocity_x),
#                                 plane_orientation='x_axes',
#                                 slice_index=0,
#                                )
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocity_x),
#                                 plane_orientation='y_axes',
#                                 slice_index=0,
#                                 )
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocity_x),
#                                 plane_orientation='z_axes',
#                                 slice_index=0,
#                                 )
#mlab.outline()
##mlab.axes(ranges=[x_min, x_max, y_min, y_max, z_min, z_max])
#mlab.xlabel('x [m]')
#mlab.ylabel('y [m]')
#mlab.zlabel('z [m]')
#mlab.title('u')
#mlab.colorbar(orientation='vertical')
#
################################################################################
#figure = mlab.figure('v in Plane')
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocity_y),
#                                 plane_orientation='x_axes',
#                                 slice_index=0,
#                                )
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocity_y),
#                                 plane_orientation='y_axes',
#                                 slice_index=0,
#                                 )
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocity_y),
#                                 plane_orientation='z_axes',
#                                 slice_index=0,
#                                 )
#mlab.outline()
##mlab.axes(ranges=[x_min, x_max, y_min, y_max, z_min, z_max])
#mlab.xlabel('x [m]')
#mlab.ylabel('y [m]')
#mlab.zlabel('z [m]')
#mlab.title('v')
#mlab.colorbar(orientation='vertical')
#
################################################################################
#figure = mlab.figure('w in Plane')
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocity_z),
#                                 plane_orientation='x_axes',
#                                 slice_index=0,
#                                )
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocity_z),
#                                 plane_orientation='y_axes',
#                                 slice_index=0,
#                                 )
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocity_z),
#                                 plane_orientation='z_axes',
#                                 slice_index=0,
#                                 )
#mlab.outline()
##mlab.axes(ranges=[x_min, x_max, y_min, y_max, z_min, z_max])
#mlab.xlabel('x [m]')
#mlab.ylabel('y [m]')
#mlab.zlabel('z [m]')
#mlab.title('w')
#mlab.colorbar(orientation='vertical')
#
################################################################################
#figure = mlab.figure('velocity in Plane')
#src = mlab.pipeline.vector_field(velocity_x, velocity_y, velocity_z)
##mlab.pipeline.vector_cut_plane(src)
#mlab.pipeline.vector_cut_plane(src, scale_factor=3)
##mlab.pipeline.vector_cut_plane(src, mask_points=1000000, scale_factor=5)
##mlab.axes(ranges=[x_min, x_max, y_min, y_max, z_min, z_max])
##mlab.axes(ranges=[x_min, x_max, y_min, y_max, z_min, z_max], extent=[1, mesh.Nx, 1, mesh.Ny, 1, mesh.Nz])
#mlab.xlabel('x [m]')
#mlab.ylabel('y [m]')
#mlab.zlabel('z [m]')
#mlab.title('velocity')
#mlab.vectorbar(orientation='vertical')
#
##################################################################################
##figure = mlab.figure('du_dx in Plane')
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocityGradientX_u),
##                                 plane_orientation='x_axes',
##                                 slice_index=0,
##                                )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocityGradientX_u),
##                                 plane_orientation='y_axes',
##                                 slice_index=0,
##                                 )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocityGradientX_u),
##                                 plane_orientation='z_axes',
##                                 slice_index=0,
##                                 )
##mlab.outline()
###mlab.axes(ranges=[x_min, x_max, y_min, y_max, z_min, z_max])
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('du_dx')
##mlab.colorbar(orientation='vertical')
##
##figure = mlab.figure('dv_dx in Plane')
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocityGradientX_v),
##                                 plane_orientation='x_axes',
##                                 slice_index=0,
##                                )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocityGradientX_v),
##                                 plane_orientation='y_axes',
##                                 slice_index=0,
##                                 )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocityGradientX_v),
##                                 plane_orientation='z_axes',
##                                 slice_index=0,
##                                 )
##mlab.outline()
###mlab.axes(ranges=[x_min, x_max, y_min, y_max, z_min, z_max])
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('dv_dx')
##mlab.colorbar(orientation='vertical')
##
##figure = mlab.figure('dw_dx in Plane')
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocityGradientX_w),
##                                 plane_orientation='x_axes',
##                                 slice_index=0,
##                                )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocityGradientX_w),
##                                 plane_orientation='y_axes',
##                                 slice_index=0,
##                                 )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocityGradientX_w),
##                                 plane_orientation='z_axes',
##                                 slice_index=0,
##                                 )
##mlab.outline()
###mlab.axes(ranges=[x_min, x_max, y_min, y_max, z_min, z_max])
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('dw_dx')
##mlab.colorbar(orientation='vertical')
##
##figure = mlab.figure('du_dy in Plane')
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocityGradientY_u),
##                                 plane_orientation='x_axes',
##                                 slice_index=0,
##                                )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocityGradientY_u),
##                                 plane_orientation='y_axes',
##                                 slice_index=0,
##                                 )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocityGradientY_u),
##                                 plane_orientation='z_axes',
##                                 slice_index=0,
##                                 )
##mlab.outline()
###mlab.axes(ranges=[x_min, x_max, y_min, y_max, z_min, z_max])
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('du_dy')
##mlab.colorbar(orientation='vertical')
##
##figure = mlab.figure('dv_dy in Plane')
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocityGradientY_v),
##                                 plane_orientation='x_axes',
##                                 slice_index=0,
##                                )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocityGradientY_v),
##                                 plane_orientation='y_axes',
##                                 slice_index=0,
##                                 )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocityGradientY_v),
##                                 plane_orientation='z_axes',
##                                 slice_index=0,
##                                 )
##mlab.outline()
###mlab.axes(ranges=[x_min, x_max, y_min, y_max, z_min, z_max])
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('dv_dy')
##mlab.colorbar(orientation='vertical')
##
##figure = mlab.figure('dw_dy in Plane')
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocityGradientY_w),
##                                 plane_orientation='x_axes',
##                                 slice_index=0,
##                                )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocityGradientY_w),
##                                 plane_orientation='y_axes',
##                                 slice_index=0,
##                                 )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocityGradientY_w),
##                                 plane_orientation='z_axes',
##                                 slice_index=0,
##                                 )
##mlab.outline()
###mlab.axes(ranges=[x_min, x_max, y_min, y_max, z_min, z_max])
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('dw_dy')
##mlab.colorbar(orientation='vertical')
##
##figure = mlab.figure('du_dz in Plane')
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocityGradientZ_u),
##                                 plane_orientation='x_axes',
##                                 slice_index=0,
##                                )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocityGradientZ_u),
##                                 plane_orientation='y_axes',
##                                 slice_index=0,
##                                 )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocityGradientZ_u),
##                                 plane_orientation='z_axes',
##                                 slice_index=0,
##                                 )
##mlab.outline()
###mlab.axes(ranges=[x_min, x_max, y_min, y_max, z_min, z_max])
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('du_dz')
##mlab.colorbar(orientation='vertical')
##
##figure = mlab.figure('dv_dz in Plane')
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocityGradientZ_v),
##                                 plane_orientation='x_axes',
##                                 slice_index=0,
##                                )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocityGradientZ_v),
##                                 plane_orientation='y_axes',
##                                 slice_index=0,
##                                 )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocityGradientZ_v),
##                                 plane_orientation='z_axes',
##                                 slice_index=0,
##                                 )
##mlab.outline()
###mlab.axes(ranges=[x_min, x_max, y_min, y_max, z_min, z_max])
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('dv_dz')
##mlab.colorbar(orientation='vertical')
##
##figure = mlab.figure('dw_dz in Plane')
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocityGradientZ_w),
##                                 plane_orientation='x_axes',
##                                 slice_index=0,
##                                )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocityGradientZ_w),
##                                 plane_orientation='y_axes',
##                                 slice_index=0,
##                                 )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(velocityGradientZ_w),
##                                 plane_orientation='z_axes',
##                                 slice_index=0,
##                                 )
##mlab.outline()
###mlab.axes(ranges=[x_min, x_max, y_min, y_max, z_min, z_max])
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('dw_dz')
##mlab.colorbar(orientation='vertical')
#
#################################################################################
##figure = mlab.figure('Du_Dx in Plane')
##src = mlab.pipeline.vector_field(velocityGradientX_u, velocityGradientX_v, velocityGradientX_w)
##mlab.pipeline.vector_cut_plane(src, scale_factor=3)
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('VelocityGradientX')
##mlab.vectorbar(orientation='vertical')
##
##figure = mlab.figure('Dv_Dy in Plane')
##src = mlab.pipeline.vector_field(velocityGradientY_u, velocityGradientY_v, velocityGradientY_w)
##mlab.pipeline.vector_cut_plane(src, scale_factor=3)
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('VelocityGradientY')
##mlab.vectorbar(orientation='vertical')
##
##
##figure = mlab.figure('Dw_Dz in Plane')
##src = mlab.pipeline.vector_field(velocityGradientZ_u, velocityGradientZ_v, velocityGradientZ_w)
##mlab.pipeline.vector_cut_plane(src, scale_factor=3)
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('VelocityGradientZ')
##mlab.vectorbar(orientation='vertical')
#
#
#################################################################################
##figure = mlab.figure('rhoVelocity in Plane')
##src = mlab.pipeline.vector_field(rhoVelocity_x, rhoVelocity_y, rhoVelocity_z)
##mlab.pipeline.vector_cut_plane(src, scale_factor=3)
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('velocity')
##mlab.vectorbar(orientation='vertical')
#
#################################################################################
##figure = mlab.figure('rhoEnergy in Plane')
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rhoEnergy),
##                                 plane_orientation='x_axes',
##                                 slice_index=0,
##                                )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rhoEnergy),
##                                 plane_orientation='y_axes',
##                                 slice_index=0,
##                                 )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rhoEnergy),
##                                 plane_orientation='z_axes',
##                                 slice_index=0,
##                                 )
##mlab.outline()
###mlab.axes(ranges=[x_min, x_max, y_min, y_max, z_min, z_max])
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('rhoEnergy')
##mlab.colorbar(orientation='vertical')
#
##################################################################################
##figure = mlab.figure('rhoVelocityFluxX in Plane')
##src = mlab.pipeline.vector_field(rhoVelocityFluxX_x, rhoVelocityFluxX_y, rhoVelocityFluxX_z)
##mlab.pipeline.vector_cut_plane(src, scale_factor=3)
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('rhoVelocityFluxX')
##mlab.vectorbar(orientation='vertical')
##
##figure = mlab.figure('rhoVelocityFluxX_x in Plane')
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rhoVelocityFluxX_x),
##                                 plane_orientation='x_axes',
##                                 slice_index=0,
##                                )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rhoVelocityFluxX_x),
##                                 plane_orientation='y_axes',
##                                 slice_index=0,
##                                 )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rhoVelocityFluxX_x),
##                                 plane_orientation='z_axes',
##                                 slice_index=0,
##                                 )
##mlab.outline()
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('rhoVelocityFluxX_x')
##mlab.colorbar(orientation='vertical')
#
#################################################################################
##figure = mlab.figure('rhoVelocityFluxY in Plane')
##src = mlab.pipeline.vector_field(rhoVelocityFluxY_x, rhoVelocityFluxY_y, rhoVelocityFluxY_z)
##mlab.pipeline.vector_cut_plane(src, scale_factor=3)
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('rhoVelocityFluxY')
##mlab.vectorbar(orientation='vertical')
##
##figure = mlab.figure('rhoVelocityFluxY_y in Plane')
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rhoVelocityFluxY_y),
##                                 plane_orientation='x_axes',
##                                 slice_index=0,
##                                )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rhoVelocityFluxY_y),
##                                 plane_orientation='y_axes',
##                                 slice_index=0,
##                                 )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rhoVelocityFluxY_y),
##                                 plane_orientation='z_axes',
##                                 slice_index=0,
##                                 )
##mlab.outline()
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('rhoVelocityFluxY_y')
##mlab.colorbar(orientation='vertical')
##
#################################################################################
##figure = mlab.figure('rhoVelocityFluxZ in Plane')
##src = mlab.pipeline.vector_field(rhoVelocityFluxZ_x, rhoVelocityFluxZ_y, rhoVelocityFluxZ_z)
##mlab.pipeline.vector_cut_plane(src, scale_factor=3)
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('rhoVelocityFluxZ')
##mlab.vectorbar(orientation='vertical')
##
##figure = mlab.figure('rhoVelocityFluxZ_z in Plane')
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rhoVelocityFluxZ_z),
##                                 plane_orientation='x_axes',
##                                 slice_index=0,
##                                )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rhoVelocityFluxZ_z),
##                                 plane_orientation='y_axes',
##                                 slice_index=0,
##                                 )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rhoVelocityFluxZ_z),
##                                 plane_orientation='z_axes',
##                                 slice_index=0,
##                                 )
##mlab.outline()
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('rhoVelocityFluxZ_z')
##mlab.colorbar(orientation='vertical')
#
###############################################################################
#figure = mlab.figure('rho_t in Plane')
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rho_t),
#                                 plane_orientation='x_axes',
#                                 slice_index=0,
#                                )
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rho_t),
#                                 plane_orientation='y_axes',
#                                 slice_index=0,
#                                 )
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rho_t),
#                                 plane_orientation='z_axes',
#                                 slice_index=0,
#                                 )
#mlab.outline()
#mlab.xlabel('x [m]')
#mlab.ylabel('y [m]')
#mlab.zlabel('z [m]')
#mlab.title('rho_t')
#mlab.colorbar(orientation='vertical')
#
################################################################################
#figure = mlab.figure('velocity_t in Plane')
#src = mlab.pipeline.vector_field(rhoVelocity_t_x, rhoVelocity_t_y, rhoVelocity_t_z)
#mlab.pipeline.vector_cut_plane(src, scale_factor=3)
#mlab.xlabel('x [m]')
#mlab.ylabel('y [m]')
#mlab.zlabel('z [m]')
#mlab.title('velocity_t')
#mlab.vectorbar(orientation='vertical')
#
#figure = mlab.figure('Velocity_t_x in Plane')
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rhoVelocity_t_x),
#                                 plane_orientation='x_axes',
#                                 slice_index=0,
#                                )
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rhoVelocity_t_x),
#                                 plane_orientation='y_axes',
#                                 slice_index=0,
#                                 )
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rhoVelocity_t_x),
#                                 plane_orientation='z_axes',
#                                 slice_index=0,
#                                 )
#mlab.outline()
#mlab.xlabel('x [m]')
#mlab.ylabel('y [m]')
#mlab.zlabel('z [m]')
#mlab.title('rhoVelocity_t_x')
#mlab.colorbar(orientation='vertical')
#
#figure = mlab.figure('Velocity_t_y in Plane')
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rhoVelocity_t_y),
#                                 plane_orientation='x_axes',
#                                 slice_index=0,
#                                )
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rhoVelocity_t_y),
#                                 plane_orientation='y_axes',
#                                 slice_index=0,
#                                 )
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rhoVelocity_t_y),
#                                 plane_orientation='z_axes',
#                                 slice_index=0,
#                                 )
#mlab.outline()
#mlab.xlabel('x [m]')
#mlab.ylabel('y [m]')
#mlab.zlabel('z [m]')
#mlab.title('rhoVelocity_t_y')
#mlab.colorbar(orientation='vertical')
#
#figure = mlab.figure('Velocity_t_z in Plane')
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rhoVelocity_t_z),
#                                 plane_orientation='x_axes',
#                                 slice_index=0,
#                                )
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rhoVelocity_t_z),
#                                 plane_orientation='y_axes',
#                                 slice_index=0,
#                                 )
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rhoVelocity_t_z),
#                                 plane_orientation='z_axes',
#                                 slice_index=0,
#                                 )
#mlab.outline()
#mlab.xlabel('x [m]')
#mlab.ylabel('y [m]')
#mlab.zlabel('z [m]')
#mlab.title('rhoVelocity_t_z')
#mlab.colorbar(orientation='vertical')
#
################################################################################
#figure = mlab.figure('rhoEnergy_t in Plane')
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rhoEnergy_t),
#                                 plane_orientation='x_axes',
#                                 slice_index=0,
#                                )
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rhoEnergy_t),
#                                 plane_orientation='y_axes',
#                                 slice_index=0,
#                                 )
#mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rhoEnergy_t),
#                                 plane_orientation='z_axes',
#                                 slice_index=0,
#                                 )
#mlab.outline()
#mlab.xlabel('x [m]')
#mlab.ylabel('y [m]')
#mlab.zlabel('z [m]')
#mlab.title('rhoEnergy_t')
#mlab.colorbar(orientation='vertical')
#
################################################################################
##figure = mlab.figure('rhoEnthalpy in Plane')
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rhoEnthalpy),
##                                 plane_orientation='x_axes',
##                                 slice_index=0,
##                                )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rhoEnthalpy),
##                                 plane_orientation='y_axes',
##                                 slice_index=0,
##                                 )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(rhoEnthalpy),
##                                 plane_orientation='z_axes',
##                                 slice_index=0,
##                                 )
##mlab.outline()
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('rhoVelocity_t_x')
##mlab.colorbar(orientation='vertical')
#
################################################################################
##figure = mlab.figure('debug_scalar_1 in Plane')
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_scalar_1),
##                                 plane_orientation='x_axes',
##                                 slice_index=0,
##                                )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_scalar_1),
##                                 plane_orientation='y_axes',
##                                 slice_index=0,
##                                 )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_scalar_1),
##                                 plane_orientation='z_axes',
##                                 slice_index=0,
##                                 )
##mlab.outline()
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('debug_scalar_1')
##mlab.colorbar(orientation='vertical')
##
################################################################################
##figure = mlab.figure('debug_scalar_2 in Plane')
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_scalar_2),
##                                 plane_orientation='x_axes',
##                                 slice_index=0,
##                                )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_scalar_2),
##                                 plane_orientation='y_axes',
##                                 slice_index=0,
##                                 )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_scalar_2),
##                                 plane_orientation='z_axes',
##                                 slice_index=0,
##                                 )
##mlab.outline()
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('debug_scalar_2')
##mlab.colorbar(orientation='vertical')
##
################################################################################
##figure = mlab.figure('debug_scalar_3 in Plane')
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_scalar_3),
##                                 plane_orientation='x_axes',
##                                 slice_index=0,
##                                )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_scalar_3),
##                                 plane_orientation='y_axes',
##                                 slice_index=0,
##                                 )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_scalar_3),
##                                 plane_orientation='z_axes',
##                                 slice_index=0,
##                                 )
##mlab.outline()
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('debug_scalar_3')
##mlab.colorbar(orientation='vertical')
#
################################################################################
##figure = mlab.figure('debug_vector_1 in Plane')
##src = mlab.pipeline.vector_field(debug_vector_1_x, debug_vector_1_y, debug_vector_1_z)
##mlab.pipeline.vector_cut_plane(src, scale_factor=3)
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('debug_vector_1')
##mlab.vectorbar(orientation='vertical')
##
##figure = mlab.figure('debug_vector_2 in Plane')
##src = mlab.pipeline.vector_field(debug_vector_2_x, debug_vector_2_y, debug_vector_2_z)
##mlab.pipeline.vector_cut_plane(src, scale_factor=3)
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('debug_vector_2')
##mlab.vectorbar(orientation='vertical')
#
##figure = mlab.figure('debug_vector_3 in Plane')
##src = mlab.pipeline.vector_field(debug_vector_3_x, debug_vector_3_y, debug_vector_3_z)
##mlab.pipeline.vector_cut_plane(src, scale_factor=3)
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('debug_vector_3')
##mlab.vectorbar(orientation='vertical')
#
#################################################################################
##figure = mlab.figure('debug_vector_1_x in Plane')
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_vector_1_x),
##                                 plane_orientation='x_axes',
##                                 slice_index=0,
##                                )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_vector_1_x),
##                                 plane_orientation='y_axes',
##                                 slice_index=0,
##                                 )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_vector_1_x),
##                                 plane_orientation='z_axes',
##                                 slice_index=0,
##                                 )
##mlab.outline()
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('debug_vector_1_x')
##mlab.colorbar(orientation='vertical')
##
##figure = mlab.figure('debug_vector_1_y in Plane')
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_vector_1_y),
##                                 plane_orientation='x_axes',
##                                 slice_index=0,
##                                )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_vector_1_y),
##                                 plane_orientation='y_axes',
##                                 slice_index=0,
##                                 )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_vector_1_y),
##                                 plane_orientation='z_axes',
##                                 slice_index=0,
##                                 )
##mlab.outline()
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('debug_vector_1_y')
##mlab.colorbar(orientation='vertical')
##
##figure = mlab.figure('debug_vector_1_z in Plane')
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_vector_1_z),
##                                 plane_orientation='x_axes',
##                                 slice_index=0,
##                                )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_vector_1_z),
##                                 plane_orientation='y_axes',
##                                 slice_index=0,
##                                 )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_vector_1_z),
##                                 plane_orientation='z_axes',
##                                 slice_index=0,
##                                 )
##mlab.outline()
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('debug_vector_1_z')
##mlab.colorbar(orientation='vertical')
##
#################################################################################
##figure = mlab.figure('debug_vector_2_x in Plane')
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_vector_2_x),
##                                 plane_orientation='x_axes',
##                                 slice_index=0,
##                                )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_vector_2_x),
##                                 plane_orientation='y_axes',
##                                 slice_index=0,
##                                 )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_vector_2_x),
##                                 plane_orientation='z_axes',
##                                 slice_index=0,
##                                 )
##mlab.outline()
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('debug_vector_2_x')
##mlab.colorbar(orientation='vertical')
##
##figure = mlab.figure('debug_vector_2_y in Plane')
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_vector_2_y),
##                                 plane_orientation='x_axes',
##                                 slice_index=0,
##                                )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_vector_2_y),
##                                 plane_orientation='y_axes',
##                                 slice_index=0,
##                                 )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_vector_2_y),
##                                 plane_orientation='z_axes',
##                                 slice_index=0,
##                                 )
##mlab.outline()
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('debug_vector_2_y')
##mlab.colorbar(orientation='vertical')
##
##figure = mlab.figure('debug_vector_2_z in Plane')
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_vector_2_z),
##                                 plane_orientation='x_axes',
##                                 slice_index=0,
##                                )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_vector_2_z),
##                                 plane_orientation='y_axes',
##                                 slice_index=0,
##                                 )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_vector_2_z),
##                                 plane_orientation='z_axes',
##                                 slice_index=0,
##                                 )
##mlab.outline()
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('debug_vector_2_z')
##mlab.colorbar(orientation='vertical')
##
#################################################################################
##figure = mlab.figure('debug_vector_3_x in Plane')
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_vector_3_x),
##                                 plane_orientation='x_axes',
##                                 slice_index=0,
##                                )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_vector_3_x),
##                                 plane_orientation='y_axes',
##                                 slice_index=0,
##                                 )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_vector_3_x),
##                                 plane_orientation='z_axes',
##                                 slice_index=0,
##                                 )
##mlab.outline()
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('debug_vector_3_x')
##mlab.colorbar(orientation='vertical')
##
##figure = mlab.figure('debug_vector_3_y in Plane')
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_vector_3_y),
##                                 plane_orientation='x_axes',
##                                 slice_index=0,
##                                )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_vector_3_y),
##                                 plane_orientation='y_axes',
##                                 slice_index=0,
##                                 )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_vector_3_y),
##                                 plane_orientation='z_axes',
##                                 slice_index=0,
##                                 )
##mlab.outline()
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('debug_vector_3_y')
##mlab.colorbar(orientation='vertical')
##
##figure = mlab.figure('debug_vector_3_z in Plane')
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_vector_3_z),
##                                 plane_orientation='x_axes',
##                                 slice_index=0,
##                                )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_vector_3_z),
##                                 plane_orientation='y_axes',
##                                 slice_index=0,
##                                 )
##mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(debug_vector_3_z),
##                                 plane_orientation='z_axes',
##                                 slice_index=0,
##                                 )
##mlab.outline()
##mlab.xlabel('x [m]')
##mlab.ylabel('y [m]')
##mlab.zlabel('z [m]')
##mlab.title('debug_vector_3_z')
##mlab.colorbar(orientation='vertical')
#
################################################################################
#mlab.show()
