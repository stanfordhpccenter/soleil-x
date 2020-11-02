import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
import bisect
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import sys
from matplotlib.patches import Rectangle

# --------------------------------------------------------------------------- #
#                                 User Input                                   #
# --------------------------------------------------------------------------- #
# list of things to save
scalar_data_to_save = ['rho','pressure','temperature']
vector_data_to_save = ['velocity']

# list of things to plot
#scalar_data_to_plot = ['rho','pressure','temperature']
scalar_data_to_plot = ['temperature']
vector_data_to_plot = ['velocity']
#vector_data_to_plot = ['cellWidth','velocity']

# Width of the domain used to obtain the spanwise averages (will be centered in simulation domain)
duct_width = 0.04 #[m]
L_spanwise_average_window = 0.075 #[m]

# --------------------------------------------------------------------------- #
#                             Utility Functions                               #
# --------------------------------------------------------------------------- #
def print_banner(message):
  print(80*'#')
  print('#{:^78s}#'.format(message))
  print(80*'#')

def print_error_message(message):
  print_banner('ERROR')
  print('{}'.format(message))

def run_command(command):
    if args.debug:
      print('Would run command:')
      print(command)
      print('')
    else:
      try:
        output=subprocess.check_output(command, shell=True)
        #print('Running command:')
        #print('{}'.format(command))
      except subprocess.CalledProcessError as e:
        print('Failed command with output:')
        print('{}'.format(command))
        print(e.output)
        sys.exit()
      else:
        print('Successfully ran command:')
        print('{}'.format(command))
        print('')
        print('With output:')
        print(output.decode('utf-8'))
        print('')

# --------------------------------------------------------------------------- #
#                            Read Command Line Input                          #
# --------------------------------------------------------------------------- #

parser = argparse.ArgumentParser()
parser.add_argument('--indir', nargs='?', const='.', default='.',
                    help='directory that contains the postprocessed soleil-x hdf5 files. basedir usually has files: fluiditer0000.hdf, fluiditer0001.hdf, particlesiter0000.hdf, particlesiter0001.hdf, ...')
parser.add_argument('--outfile', nargs='?', const='default', default='default',
                    help='directory where output is written')
parser.add_argument('-v', '--verbose',
                    action='store_true',
                    help='verbose output')
parser.add_argument('--debug',
                    action='store_true',
                    help='run in debug mode')
args = parser.parse_args()

# hack way to set default of outdir based on the basedir argument
output_file = args.outfile
if args.outfile == 'default':
  output_file = os.path.join(args.indir,'reduced_data.hdf')

# Turn pathnames into absolute paths. Doing this because later in this script os.walk() can have trouble with relative path names.
input_dir = os.path.abspath(args.indir)
output_file = os.path.abspath(output_file)

if args.verbose:
  print_banner('Input Summary')
  print('input  directory: {}'.format(input_dir))
  print('output file: {}'.format(output_file))
  print('')

# --------------------------------------------------------------------------- #
#                              Error Checking                                 #
# --------------------------------------------------------------------------- #
if not os.path.exists(input_dir):
  print_error_message('Error: the provided input directory {} does not exist'.format(base_dir))
  sys.exit()

# --------------------------------------------------------------------------- #
# Get list of fluid and particle data files in input directory
# --------------------------------------------------------------------------- #
files_in_input_dir = [f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, f))]

fluid_filenames = []
particle_filenames = []
for filename in files_in_input_dir:
  if ('fluid' in filename)  and ('.hdf' in filename):
    fluid_filenames.append(filename)
  if ('particles' in filename)  and ('.hdf' in filename):
    particle_filenames.append(filename)

fluid_filenames.sort()
particle_filenames.sort()

if not fluid_filenames:
  print_error_message("""Input directory ( {} ) did not contain any fluid data files. This script looks in the input directory for files with the substrings "fluid" and  ".hdf" in the filename. None were found.
   This script assumes a input base directory format of:
   input_dir
   |-- fluid0000000000.hdf
   |-- fluid0000000001.hdf
   |--   .
   |--   .
   |--   .
   |-- fluid##########.hdf""".format(base_dir))
  sys.exit()

if not particle_filenames:
  print_error_message("""Input directory ( {} ) did not contain any particle data files. This script looks in the input directory for files with the substrings "particle" and  ".hdf" in the filename. None were found.
   This script assumes a input base directory format of:
   input_dir
   |-- particles0000000000.hdf
   |-- particles0000000001.hdf
   |--   .
   |--   .
   |--   .
   |-- particles##########.hdf""".format(base_dir))
  sys.exit()

if args.verbose:
  print_banner('Data files in input directory')
  [print(filename) for filename in (fluid_filenames + particle_filenames)]
  print('')

## --------------------------------------------------------------------------- #
##                                                      #
## --------------------------------------------------------------------------- #

# Labels for plots
title_labels = {'rho'         : r'$\rho$',
                'pressure'    : r'$P$',
                'temperature' : r'$T$'}

unit_labels = {'rho'         : r'$\left[\frac{kg}{m^3}\right]$',
               'pressure'    : r'[Pa]',
               'temperature' : r'[K]'}

# Constant to convert meters to centimeters in the plots
m_to_cm = np.power(10,2.0)


# Open file where output will be written
if args.debug:
  print('Would create new file:')
  print('{}'.format(output_file))
else:
  if not os.path.exists(output_file):
    # hdf output filename
    hdf_out = h5py.File(output_file, 'w')
  else:
    print_error_message("""Output file exists:
    {}
    I don\'t want to clobber whatever is there.
    Delete it and re-run and this script""".format(output_file))

    sys.exit()

for filename in fluid_filenames:
  # Get time step number from the file name
  time_step = int(filename.lstrip('fluid').rstrip('.hdf'))

  # read fluid data file 
  f = h5py.File(os.path.join(input_dir,filename), 'r')
  
  # Get the data
  x_points = f['x_points']
  y_points = f['y_points']
  z_points = f['z_points']
  
  # Get number of cells in each direction
  [Nx, Ny, Nz] = np.shape(f[scalar_data_to_plot[0]])
  
  # Get length of domain in each direction
  Lx = (np.max(x_points) - np.min(x_points))
  Ly = (np.max(y_points) - np.min(y_points))
  Lz = (np.max(z_points) - np.min(z_points))

  # Get cell centers
  x_centers = (x_points[1:] + x_points[:-1] ) / 2.0
  y_centers = (y_points[1:] + y_points[:-1] ) / 2.0
  z_centers = (z_points[1:] + z_points[:-1] ) / 2.0
  
  # Get center of the domain
  x_mid = np.min(x_points) + (Lx/2.0)
  y_mid = np.min(y_points) + (Ly/2.0)
  z_mid = np.min(z_points) + (Lz/2.0)

  # Legend of plotting:
  # x = Spanwise
  # y = Spanwise and Radiation Streaming Direction
  # z = Streamwise Direction

  # Stremwise location of measurement plane used in experiments
  streamwise_location_of_measurement_plane = 0.305 # [m]
  
  x_mid_idx = bisect.bisect(x_points, x_mid)
  y_mid_idx = bisect.bisect(y_points, y_mid)
  z_mid_idx = bisect.bisect(z_points, z_mid)

  z_measure_ment_plane_idx = bisect.bisect(z_points, streamwise_location_of_measurement_plane)

  x_duct_start = x_mid - (duct_width/2.0) 
  x_duct_stop  = x_mid + (duct_width/2.0) 
  y_duct_start = y_mid - (duct_width/2.0) 
  y_duct_stop  = y_mid + (duct_width/2.0) 

  x_spanwise_average_start = x_mid - (L_spanwise_average_window/2.0) 
  x_spanwise_average_stop  = x_mid + (L_spanwise_average_window/2.0) 
  y_spanwise_average_start = y_mid - (L_spanwise_average_window/2.0) 
  y_spanwise_average_stop  = y_mid + (L_spanwise_average_window/2.0) 

  x_spanwise_average_start_idx = bisect.bisect(x_points, x_spanwise_average_start)
  x_spanwise_average_stop_idx  = bisect.bisect(x_points, x_spanwise_average_stop)
  y_spanwise_average_start_idx = bisect.bisect(y_points, y_spanwise_average_start)
  y_spanwise_average_stop_idx  = bisect.bisect(y_points, y_spanwise_average_stop)

  if args.verbose:
    print('[Nx, Ny, Nz] = [{}, {}, {}]'.format(Nx, Ny, Nz))
    print('[Lx, Ly, Lz] = [{}, {}, {}]'.format(Lx, Ly, Lz))

  # Save data to an hdf output file
  hdf_out.attrs['L_spanwise_average_window '] = L_spanwise_average_window 
  hdf_out.attrs['duct_width'] = duct_width 

  hdf_out.attrs['x_duct_start'] = x_duct_start   
  hdf_out.attrs['x_duct_stop'] = x_duct_stop   
  hdf_out.attrs['y_duct_start'] = y_duct_start  
  hdf_out.attrs['y_duct_stop'] = y_duct_stop   

  hdf_out.attrs['x_spanwise_average_start'] = x_spanwise_average_start
  hdf_out.attrs['x_spanwise_average_stop'] = x_spanwise_average_stop 
  hdf_out.attrs['y_spanwise_average_start'] = y_spanwise_average_start
  hdf_out.attrs['y_spanwise_average_stop'] = y_spanwise_average_stop 

  hdf_out.attrs['x_spanwise_average_start_idx'] = x_spanwise_average_start_idx
  hdf_out.attrs['x_spanwise_average_stop_idx'] = x_spanwise_average_stop_idx 
  hdf_out.attrs['y_spanwise_average_start_idx'] = y_spanwise_average_start_idx
  hdf_out.attrs['y_spanwise_average_stop_idx'] = y_spanwise_average_stop_idx 

  time_step_group = hdf_out.create_group("timestep_{}".format(time_step))

  for feild_name in (scalar_data_to_save + vector_data_to_save): 
    ## Streamwise Midslice `
    #hdf_out['mid_x_slice_{}'.format(feild_name)] = f['{}'.format(feild_name)][x_mid_idx,:,:]

    ## Spanwise slices Measurement plane slice
    #hdf_out['measurement_z_slice_{}'.format(feild_name)] = f['{}'.format(feild_name)][:,:,z_measure_ment_plane_idx]
    #hdf_out['inlet_z_slice_{}'.format(feild_name)] = f['{}'.format(feild_name)][:,:,0]
    #hdf_out['outlet_z_slice_{}'.format(feild_name)] = f['{}'.format(feild_name)][:,:,-1]

    #hdf_out['{}_spanwise_average'.format(feild_name)] = np.mean(f['{}'.format(feild_name)][x_spanwise_average_start_idx:x_spanwise_average_stop_idx,
    #                                                                                       y_spanwise_average_start_idx:y_spanwise_average_stop_idx,
    #                                                                                       :], axis=(0,1))
   
    # Streamwise Midslice
    time_step_group['mid_x_slice_{}'.format(feild_name)] = f['{}'.format(feild_name)][x_mid_idx,:,:]

    # Spanwise slices Measurement plane slice
    time_step_group['measurement_z_slice_{}'.format(feild_name)] = f['{}'.format(feild_name)][:,:,z_measure_ment_plane_idx]
    time_step_group['inlet_z_slice_{}'.format(feild_name)] = f['{}'.format(feild_name)][:,:,0]
    time_step_group['outlet_z_slice_{}'.format(feild_name)] = f['{}'.format(feild_name)][:,:,-1]

    time_step_group['{}_spanwise_average'.format(feild_name)] = np.mean(f['{}'.format(feild_name)][x_spanwise_average_start_idx:x_spanwise_average_stop_idx,
                                                                                           y_spanwise_average_start_idx:y_spanwise_average_stop_idx,
                                                                                           :], axis=(0,1))
  f.close()


#   # Plots
#  for scalar_feild_name in scalar_data_to_plot: 
#    plt.figure()
#    plt.plot(f['{}'.format(scalar_feild_name)][x_mid_idx, y_spanwise_average_start_idx:y_spanwise_average_stop_idx, z_measure_ment_plane_idx],
#             y_centers[y_spanwise_average_start_idx:y_spanwise_average_stop_idx],
#             '-k')
#    plt.xlabel('{} {}'.format(title_labels[scalar_feild_name], unit_labels[scalar_feild_name]), fontsize = 20)
#    plt.ylabel(r'$y$ [m]', fontsize = 20)
#    plt.title('{} on Measurement Line'.format(title_labels[scalar_feild_name]), fontsize = 20)
#    plt.savefig('{}_measurement_line.png'.format(scalar_feild_name), bbox_inches='tight')
#   
#  for scalar_feild_name in scalar_data_to_plot: 
#    plt.figure(figsize=(11,9))
#
#    #X,Y = np.meshgrid(y_points[1:-1], z_points[1:-1])
#    #plt.pcolormesh(X,
#    #               Y,
#    #               np.array(f['{}'.format(scalar_feild_name)][x_slice_idx,1:-1,1:-1]))
#  
#    #X,Y = np.meshgrid(y_points[:], z_points[:])
#    #plt.pcolormesh(X,
#    #               Y,
#    #               np.array(f['{}'.format(scalar_feild_name)][x_slice_idx,:,:]))
#  
#    #plt.pcolormesh(z_points,
#    #               y_points,
#    #               np.array(f['{}'.format(scalar_feild_name)][x_slice_idx,:,:]))
#  
#    #plt.pcolormesh(z_points*m_to_cm,
#    #               y_points[y_slice_idx_start:y_slice_idx_stop]*m_to_cm,
#    #               np.array(f['{}'.format(scalar_feild_name)][x_slice_idx,y_slice_idx_start:y_slice_idx_stop,:]))
#
#    plt.pcolormesh(z_points,
#                   y_points,
#                   time_step_group['mid_x_slice_{}'.format(scalar_feild_name)])
#
#    # Add where duct would be
#    rect = Rectangle((0.0, y_duct_start), Lz, duct_width, edgecolor='k', facecolor='none', linestyle='--')
#    ax = plt.gca()
#    ax.add_patch(rect)
#
#    ## Add where averaging area is
#    rect = Rectangle((0.0, y_spanwise_average_start), Lz, L_spanwise_average_window, edgecolor='k', facecolor='none')
#    ax = plt.gca()
#    ax.add_patch(rect)
#
#    plt.axis('scaled')
#    plt.xlim([0.0,Lz])
#    plt.ylim([0.0,Ly])
#    plt.xlabel(r'$z \ [cm]$', fontsize = 20)
#    plt.ylabel(r'$y \ [cm]$', fontsize = 20)
#    plt.title('{} {}'.format(title_labels[scalar_feild_name], unit_labels[scalar_feild_name]), fontsize = 20)
#    # create an axes on the right side of ax. The width of cax will be 5%
#    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
#    divider = make_axes_locatable(plt.gca())
#    cax = divider.append_axes("right", size="2%", pad=0.05)
#    plt.colorbar(cax=cax)
#
#    plt.savefig('{}_slice.png'.format(scalar_feild_name), dpi=300, bbox_inches='tight')
#
#  for scalar_feild_name in scalar_data_to_plot: 
#    plt.figure(figsize=(11,9))
#    #plt.pcolormesh(y_points,
#    #               x_points,
#    #               time_step_group['measurement_z_slice_{}'.format(scalar_feild_name)] )
#    plt.pcolormesh(y_points,
#                   x_points,
#                   time_step_group['outlet_z_slice_{}'.format(scalar_feild_name)] )
#
#    # Add where duct would be
#    rect = Rectangle((x_duct_start, x_duct_start), duct_width, duct_width, edgecolor='k', facecolor='none', linestyle='--')
#    ax = plt.gca()
#    ax.add_patch(rect)
#
#    # Add where averaging area is
#    rect = Rectangle((x_spanwise_average_start, x_spanwise_average_start), L_spanwise_average_window, L_spanwise_average_window, edgecolor='k', facecolor='none')
#    ax = plt.gca()
#    ax.add_patch(rect)
#
#    plt.axis('scaled')
#    plt.xlabel(r'$z \ [cm]$', fontsize = 20)
#    plt.ylabel(r'$y \ [cm]$', fontsize = 20)
#    plt.title('{} {}'.format(title_labels[scalar_feild_name], unit_labels[scalar_feild_name]), fontsize = 20)
#    # create an axes on the right side of ax. The width of cax will be 5%
#    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
#    divider = make_axes_locatable(plt.gca())
#    cax = divider.append_axes("right", size="2%", pad=0.05)
#    plt.colorbar(cax=cax)
#    plt.savefig('measurement_plane_{}_slice.png'.format(scalar_feild_name), dpi=300, bbox_inches='tight')
#
#  for scalar_feild_name in scalar_data_to_plot: 
#    plt.figure()
#    plt.plot(z_centers,
#             time_step_group['{}_spanwise_average'.format(scalar_feild_name)],
#             '-k',
#             label='Soleil-X')
#    plt.xlabel(r'$z \ [m]$', fontsize = 20)
#    plt.ylabel('{} {}'.format(title_labels[scalar_feild_name], unit_labels[scalar_feild_name]), fontsize = 20)
#    plt.title('{} Spanwise Average'.format(title_labels[scalar_feild_name]), fontsize = 20)
#    plt.savefig('{}_spanwise_average.png'.format(scalar_feild_name), bbox_inches='tight')
#    
#  
#  #for vector_feild_name in vector_data_to_plot: 
#  #  for vector_component in ['x','y','z']:
#  #    idx = {'x':0, 'y':1, 'z':2}
#  #    plt.figure()
#  #    plt.plot(f['{}'.format(vector_feild_name)][z_slice_idx,:,x_slice_idx][:,idx[vector_component]], y_values, 'ok', label='Soleil-X')
#  #
#  #    x_label = ''
#  #    if vector_feild_name == 'velocityGradientX':
#  #      if vector_component == 'x':
#  #        x_label = r'$\frac{{\partial {}}}{{ \partial {} }}$'.format('u','x')
#  #      elif vector_component == 'y':
#  #        x_label = r'$\frac{{\partial {}}}{{ \partial {} }}$'.format('v','x')
#  #      elif vector_component == 'z':
#  #        x_label = r'$\frac{{\partial {}}}{{ \partial {} }}$'.format('w','x')
#  #    elif vector_feild_name == 'velocityGradientY':
#  #      if vector_component == 'x':
#  #        x_label = r'$\frac{{\partial {}}}{{ \partial {} }}$'.format('u','y')
#  #      elif vector_component == 'y':
#  #        x_label = r'$\frac{{\partial {}}}{{ \partial {} }}$'.format('v','y')
#  #      elif vector_component == 'z':
#  #        x_label = r'$\frac{{\partial {}}}{{ \partial {} }}$'.format('w','y')
#  #    elif vector_feild_name == 'velocityGradientZ':
#  #      if vector_component == 'x':
#  #        x_label = r'$\frac{{\partial {}}}{{ \partial {} }}$'.format('u','Z')
#  #      elif vector_component == 'y':
#  #        x_label = r'$\frac{{\partial {}}}{{ \partial {} }}$'.format('v','Z')
#  #      elif vector_component == 'z':
#  #        x_label = r'$\frac{{\partial {}}}{{ \partial {} }}$'.format('w','Z')
#  #    else:
#  #        x_label = '{}_{}'.format(vector_feild_name, vector_component)
#  #
#  #    plt.title(x_label)
#  #    plt.xlabel(x_label, fontsize = 20)
#  #    plt.ylabel(r'$y \ [m]$', fontsize = 20)
#  #    #plt.legend()
#  #    #plt.savefig('{}.pdf'.format(vector_feild_name), bbox_inches='tight')
#
#plt.show()

hdf_out.close()
