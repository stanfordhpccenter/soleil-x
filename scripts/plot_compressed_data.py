import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import bisect
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Rectangle
#import matplotlib.animation as animation
import subprocess
import shutil

# --------------------------------------------------------------------------- #
#                                User Input                                   #
# --------------------------------------------------------------------------- #
dt = 9.816128316174614e-08
nt_between_restarts = 10000

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
                    help='directory that contains the postprocessed and compressed soleil-x hdf5 files. indir has files: ')
parser.add_argument('--outdir', nargs='?', const='.', default='.',
                    help='directory where output is written')
parser.add_argument('-v', '--verbose',
                    action='store_true',
                    help='verbose output')
parser.add_argument('--debug',
                    action='store_true',
                    help='run in debug mode')
args = parser.parse_args()

# Turn pathnames into absolute paths. Doing this because later in this script os.walk() can have trouble with relative path names.
input_dir = os.path.abspath(args.indir)
output_dir = os.path.abspath(args.outdir)

if args.verbose:
  print_banner('Input Summary')
  print('input  directory: {}'.format(input_dir))
  print('output file: {}'.format(output_dir))
  print('')

# --------------------------------------------------------------------------- #
#                              Error Checking                                 #
# --------------------------------------------------------------------------- #
if not os.path.exists(input_dir):
  print_error_message('Error: the provided input directory {} does not exist'.format(input_dir))
  sys.exit()

if not os.path.exists(output_dir):
  print_error_message('Error: the provided outputt directory {} does not exist'.format(output_dir))
  sys.exit()

# --------------------------------------------------------------------------- #
# Get list of fluid and particle data files in input directory
# --------------------------------------------------------------------------- #
files_in_input_dir = [f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, f))]

hdf_filenames = []
for filename in files_in_input_dir:
  if '.hdf' in filename:
    hdf_filenames.append(filename)

hdf_filenames.sort()

if not hdf_filenames:
  print_error_message("""Input directory ( {} ) did not contain any hdf data files. This script looks in the input directory for files with the substring ".hdf" in the filename. None were found.
   This script assumes a input base directory format of:
   input_dir
   |-- job0.hdf
   |-- job1.hdf
   |--   .
   |--   .
   |--   .
   |-- jobN.hdf""".format(base_dir))
  sys.exit()

if args.verbose:
  print_banner('Data files in input directory')
  [print(filename) for filename in (hdf_filenames)]
  print('')


# --------------------------------------------------------------------------- #
# Pull the Data Out of the HDF Files                                          #
# --------------------------------------------------------------------------- #
# list of things to save
scalar_data_saved = ['rho','pressure','temperature']
vector_data_saved = ['velocity']

#scalar_data_to_plot = ['rho','pressure','temperature']
scalar_data_to_plot = ['temperature']
vector_data_to_plot = ['velocity']

# Labels for plots
title_labels = {'rho'         : r'$\rho$',
                'pressure'    : r'$P$',
                'temperature' : r'$T$'}

unit_labels = {'rho'         : r'$\left[\frac{kg}{m^3}\right]$',
               'pressure'    : r'[Pa]',
               'temperature' : r'[K]'}


# Read the mesh info from one of the hdf files (assume this is the same of all files)
filename = hdf_filenames[0]
h5_file  = h5py.File(os.path.join(input_dir,filename), "r")
L_spanwise_average_window  = h5_file.attrs['L_spanwise_average_window '] 
duct_width = h5_file.attrs['duct_width'] 

x_duct_start = h5_file.attrs['x_duct_start'] 
x_duct_stop = h5_file.attrs['x_duct_stop'] 
y_duct_start = h5_file.attrs['y_duct_start'] 
y_duct_stop = h5_file.attrs['y_duct_stop'] 

x_spanwise_average_start = h5_file.attrs['x_spanwise_average_start'] 
x_spanwise_average_stop = h5_file.attrs['x_spanwise_average_stop'] 
y_spanwise_average_start = h5_file.attrs['y_spanwise_average_start'] 
y_spanwise_average_stop = h5_file.attrs['y_spanwise_average_stop'] 

x_spanwise_average_start_idx = h5_file.attrs['x_spanwise_average_start_idx'] 
x_spanwise_average_stop_idx = h5_file.attrs['x_spanwise_average_stop_idx'] 
y_spanwise_average_start_idx = h5_file.attrs['y_spanwise_average_start_idx'] 
y_spanwise_average_stop_idx = h5_file.attrs['y_spanwise_average_stop_idx'] 

x_points = np.array(h5_file['timestep_0']['x_points'])
y_points = np.array(h5_file['timestep_0']['y_points'])
z_points = np.array(h5_file['timestep_0']['z_points'])

# Do some math based on that data
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

# Index of center of domain
x_mid_idx = bisect.bisect(x_points, x_mid)
y_mid_idx = bisect.bisect(y_points, y_mid)
z_mid_idx = bisect.bisect(z_points, z_mid)

h5_file.close()


###############################################################################
# Find time averages
###############################################################################

# Setup arrays for the time averages
filename = hdf_filenames[0]
h5_file  = h5py.File(os.path.join(input_dir,filename), "r")
time_averaged_data = {}
for feild_name in (scalar_data_saved + vector_data_saved): 
  time_averaged_data['mid_x_slice_{}'.format(feild_name)]         = np.zeros(np.shape(np.array(h5_file['timestep_0']['mid_x_slice_{}'.format(feild_name)])))
  time_averaged_data['mid_y_slice_{}'.format(feild_name)]         = np.zeros(np.shape(np.array(h5_file['timestep_0']['mid_y_slice_{}'.format(feild_name)])))
  time_averaged_data['measurement_z_slice_{}'.format(feild_name)] = np.zeros(np.shape(np.array(h5_file['timestep_0']['measurement_z_slice_{}'.format(feild_name)])))
  time_averaged_data['inlet_z_slice_{}'.format(feild_name)]       = np.zeros(np.shape(np.array(h5_file['timestep_0']['inlet_z_slice_{}'.format(feild_name)])))
  time_averaged_data['outlet_z_slice_{}'.format(feild_name)]      = np.zeros(np.shape(np.array(h5_file['timestep_0']['outlet_z_slice_{}'.format(feild_name)])))
  time_averaged_data['{}_spanwise_average'.format(feild_name)]    = np.zeros(np.shape(np.array(h5_file['timestep_0']['{}_spanwise_average'.format(feild_name)])))
  time_averaged_data['{}_average_0'.format(feild_name)]           = np.zeros(np.shape(np.array(h5_file['timestep_0']['{}_average_0'.format(feild_name)])))
  time_averaged_data['{}_average_1'.format(feild_name)]           = np.zeros(np.shape(np.array(h5_file['timestep_0']['{}_average_1'.format(feild_name)])))
  #time_averaged_data['{}_average_2'.format(feild_name)]          = np.zeros(np.shape(np.array(h5_file['timestep_0']['{}_average_2'.format(feild_name)])))
  time_averaged_data['measurement_line_{}'.format(feild_name)]    = np.zeros(np.shape(np.array(h5_file['timestep_0']['measurement_z_slice_{}'.format(feild_name)][x_mid_idx, :])))
  time_averaged_data['center_line_{}'.format(feild_name)]         = np.zeros(np.shape(np.array(h5_file['timestep_0']['mid_x_slice_{}'.format(feild_name)][x_mid_idx, :])))
h5_file.close()

# Accumulate over all time steps
n_timesteps_total = 0
for filename in hdf_filenames:
  h5_file  = h5py.File(os.path.join(input_dir,filename), "r")
  # Get the timesteps in the file
  for timestep_group in h5_file.keys():
    # skip timestep zero since it is the same as the first timestep of the next restart file
    if timestep_group != 'timestep_0':
      n_timesteps_total += 1 
      for feild_name in (scalar_data_saved + vector_data_saved): 
        time_averaged_data['mid_x_slice_{}'.format(feild_name)]         += np.array(h5_file[timestep_group]['mid_x_slice_{}'.format(feild_name)])
        time_averaged_data['mid_y_slice_{}'.format(feild_name)]         += np.array(h5_file[timestep_group]['mid_y_slice_{}'.format(feild_name)])
        time_averaged_data['measurement_z_slice_{}'.format(feild_name)] += np.array(h5_file[timestep_group]['measurement_z_slice_{}'.format(feild_name)])
        time_averaged_data['inlet_z_slice_{}'.format(feild_name)]       += np.array(h5_file[timestep_group]['inlet_z_slice_{}'.format(feild_name)])
        time_averaged_data['outlet_z_slice_{}'.format(feild_name)]      += np.array(h5_file[timestep_group]['outlet_z_slice_{}'.format(feild_name)])
        time_averaged_data['{}_spanwise_average'.format(feild_name)]    += np.array(h5_file[timestep_group]['{}_spanwise_average'.format(feild_name)])
        time_averaged_data['{}_average_0'.format(feild_name)]           += np.array(h5_file[timestep_group]['{}_average_0'.format(feild_name)])
        time_averaged_data['{}_average_1'.format(feild_name)]           += np.array(h5_file[timestep_group]['{}_average_1'.format(feild_name)])
        #time_averaged_data['{}_average_2'.format(feild_name)]           += np.array(h5_file[timestep_group]['{}_average_2'.format(feild_name)])
        time_averaged_data['measurement_line_{}'.format(feild_name)]    += np.array(h5_file[timestep_group]['measurement_z_slice_{}'.format(feild_name)][x_mid_idx, :])
        time_averaged_data['center_line_{}'.format(feild_name)]         += np.array(h5_file[timestep_group]['mid_x_slice_{}'.format(feild_name)][y_mid_idx, :])
  h5_file.close() 

# Average
for feild_name in (scalar_data_saved + vector_data_saved): 
  time_averaged_data['mid_x_slice_{}'.format(feild_name)]         /= n_timesteps_total 
  time_averaged_data['mid_y_slice_{}'.format(feild_name)]         /= n_timesteps_total 
  time_averaged_data['measurement_z_slice_{}'.format(feild_name)] /= n_timesteps_total 
  time_averaged_data['inlet_z_slice_{}'.format(feild_name)]       /= n_timesteps_total 
  time_averaged_data['outlet_z_slice_{}'.format(feild_name)]      /= n_timesteps_total 
  time_averaged_data['{}_spanwise_average'.format(feild_name)]    /= n_timesteps_total 
  time_averaged_data['{}_average_0'.format(feild_name)]           /= n_timesteps_total
  time_averaged_data['{}_average_1'.format(feild_name)]           /= n_timesteps_total
  #time_averaged_data['{}_average_2'.format(feild_name)]           /= n_timesteps_total
  time_averaged_data['measurement_line_{}'.format(feild_name)]    /= n_timesteps_total 
  time_averaged_data['center_line_{}'.format(feild_name)]         /= n_timesteps_total 

###############################################################################
# plot the time averages
###############################################################################

# Mid X Slice
for scalar_feild_name in scalar_data_to_plot: 
  plt.figure(figsize=(11,9))
  plt.pcolormesh(np.array(z_points),
                 np.array(y_points),
                 time_averaged_data['mid_x_slice_{}'.format(scalar_feild_name)])
  # Add where duct would be
  rect = Rectangle((0.0, y_duct_start), Lz, duct_width, edgecolor='k', facecolor='none', linestyle='--')
  ax = plt.gca()
  ax.add_patch(rect)
  ## Add where averaging area is
  rect = Rectangle((0.0, y_spanwise_average_start), Lz, L_spanwise_average_window, edgecolor='k', facecolor='none')
  ax = plt.gca()
  ax.add_patch(rect)
  # Labels
  plt.axis('scaled')
  plt.xlim([0.0,Lz])
  plt.ylim([0.0,Ly])
  plt.xlabel(r'$x \ [m]$', fontsize = 20)
  plt.ylabel(r'$y \ [m]$', fontsize = 20)
  plt.title('{} {}'.format(title_labels[scalar_feild_name], unit_labels[scalar_feild_name]), fontsize = 20)
  # create an axes on the right side of ax. The width of cax will be 5%
  # of ax and the padding between cax and ax will be fixed at 0.05 inch.
  divider = make_axes_locatable(plt.gca())
  cax = divider.append_axes("right", size="2%", pad=0.05)
  plt.colorbar(cax=cax)
  plt.savefig(os.path.join(output_dir,'{}_mid_x_slice_time_average.png'.format(scalar_feild_name)), dpi=300, bbox_inches='tight')

# Mid Y Slice
for scalar_feild_name in scalar_data_to_plot: 
  plt.figure(figsize=(11,9))
  plt.pcolormesh(np.array(z_points),
                 np.array(x_points),
                 time_averaged_data['mid_y_slice_{}'.format(scalar_feild_name)])
  # Add where duct would be
  rect = Rectangle((0.0, y_duct_start), Lz, duct_width, edgecolor='k', facecolor='none', linestyle='--')
  ax = plt.gca()
  ax.add_patch(rect)
  ## Add where averaging area is
  rect = Rectangle((0.0, y_spanwise_average_start), Lz, L_spanwise_average_window, edgecolor='k', facecolor='none')
  ax = plt.gca()
  ax.add_patch(rect)
  # Labels
  plt.axis('scaled')
  plt.xlim([0.0,Lz])
  plt.ylim([0.0,Ly])
  plt.xlabel(r'$x \ [m]$', fontsize = 20)
  plt.ylabel(r'$z \ [m]$', fontsize = 20)
  plt.title('{} {}'.format(title_labels[scalar_feild_name], unit_labels[scalar_feild_name]), fontsize = 20)
  # create an axes on the right side of ax. The width of cax will be 5%
  # of ax and the padding between cax and ax will be fixed at 0.05 inch.
  divider = make_axes_locatable(plt.gca())
  cax = divider.append_axes("right", size="2%", pad=0.05)
  plt.colorbar(cax=cax)
  plt.savefig(os.path.join(output_dir,'{}_mid_y_slice_time_average.png'.format(scalar_feild_name)), dpi=300, bbox_inches='tight')

# Average axis 0 (spanwise )
for scalar_feild_name in scalar_data_to_plot: 
  plt.figure(figsize=(11,9))
  plt.pcolormesh(np.array(z_points),
                 np.array(y_points),
                 time_averaged_data['{}_average_0'.format(scalar_feild_name)])
  # Add where duct would be
  rect = Rectangle((0.0, y_duct_start), Lz, duct_width, edgecolor='k', facecolor='none', linestyle='--')
  ax = plt.gca()
  ax.add_patch(rect)
  ## Add where averaging area is
  rect = Rectangle((0.0, y_spanwise_average_start), Lz, L_spanwise_average_window, edgecolor='k', facecolor='none')
  ax = plt.gca()
  ax.add_patch(rect)
  # Labels
  plt.axis('scaled')
  plt.xlim([0.0,Lz])
  plt.ylim([0.0,Ly])
  plt.xlabel(r'$x \ [m]$', fontsize = 20)
  plt.ylabel(r'$y \ [m]$', fontsize = 20)
  plt.title('{} {}'.format(title_labels[scalar_feild_name], unit_labels[scalar_feild_name]), fontsize = 20)
  # create an axes on the right side of ax. The width of cax will be 5%
  # of ax and the padding between cax and ax will be fixed at 0.05 inch.
  divider = make_axes_locatable(plt.gca())
  cax = divider.append_axes("right", size="2%", pad=0.05)
  plt.colorbar(cax=cax)
  plt.savefig(os.path.join(output_dir,'{}_average_0_time_average.png'.format(scalar_feild_name)), dpi=300, bbox_inches='tight')

# Average axis 1 (spanwise and direction of radiation propigation)
for scalar_feild_name in scalar_data_to_plot: 
  plt.figure(figsize=(11,9))
  plt.pcolormesh(np.array(z_points),
                 np.array(x_points),
                 time_averaged_data['{}_average_1'.format(scalar_feild_name)])
  # Add where duct would be
  rect = Rectangle((0.0, y_duct_start), Lz, duct_width, edgecolor='k', facecolor='none', linestyle='--')
  ax = plt.gca()
  ax.add_patch(rect)
  ## Add where averaging area is
  rect = Rectangle((0.0, y_spanwise_average_start), Lz, L_spanwise_average_window, edgecolor='k', facecolor='none')
  ax = plt.gca()
  ax.add_patch(rect)
  # Labels
  plt.axis('scaled')
  plt.xlim([0.0,Lz])
  plt.ylim([0.0,Ly])
  plt.xlabel(r'$x \ [m]$', fontsize = 20)
  plt.ylabel(r'$z \ [m]$', fontsize = 20)
  plt.title('{} {}'.format(title_labels[scalar_feild_name], unit_labels[scalar_feild_name]), fontsize = 20)
  # create an axes on the right side of ax. The width of cax will be 5%
  # of ax and the padding between cax and ax will be fixed at 0.05 inch.
  divider = make_axes_locatable(plt.gca())
  cax = divider.append_axes("right", size="2%", pad=0.05)
  plt.colorbar(cax=cax)
  plt.savefig(os.path.join(output_dir,'{}_average_1_time_average.png'.format(scalar_feild_name)), dpi=300, bbox_inches='tight')

  #plt.figure()
  #plt.pcolormesh(np.array(x_points),
  #               np.array(y_points),
  #               time_averaged_data['{}_average_2'.format(scalar_feild_name)])
 

# Spanwise average
for scalar_feild_name in scalar_data_to_plot: 
  plt.figure()
  plt.plot(z_centers,
           time_averaged_data['{}_spanwise_average'.format(scalar_feild_name)],
           '-k',
           label='Soleil-X')
  plt.xlabel(r'$x \ [m]$', fontsize = 20)
  plt.ylabel('{} {}'.format(title_labels[scalar_feild_name], unit_labels[scalar_feild_name]), fontsize = 20)
  plt.title('{} Spanwise Average'.format(title_labels[scalar_feild_name]), fontsize = 20)
  plt.savefig(os.path.join(output_dir,'{}_spanwise_average_time_average.png'.format(scalar_feild_name)), bbox_inches='tight')

# Measurement Plane
for scalar_feild_name in scalar_data_to_plot: 
  plt.figure(figsize=(11,9))
  plt.pcolormesh(y_points,
                 x_points,
                 time_averaged_data['measurement_z_slice_{}'.format(scalar_feild_name)] )

  # Add where duct would be
  rect = Rectangle((x_duct_start, x_duct_start), duct_width, duct_width, edgecolor='k', facecolor='none', linestyle='--')
  ax = plt.gca()
  ax.add_patch(rect)

  # Add where averaging area is
  rect = Rectangle((x_spanwise_average_start, x_spanwise_average_start), L_spanwise_average_window, L_spanwise_average_window, edgecolor='k', facecolor='none')
  ax = plt.gca()
  ax.add_patch(rect)

  plt.axis('scaled')
  plt.xlabel(r'$y \ [cm]$', fontsize = 20)
  plt.ylabel(r'$z \ [cm]$', fontsize = 20)
  plt.title('{} {}'.format(title_labels[scalar_feild_name], unit_labels[scalar_feild_name]), fontsize = 20)
  # create an axes on the right side of ax. The width of cax will be 5%
  # of ax and the padding between cax and ax will be fixed at 0.05 inch.
  divider = make_axes_locatable(plt.gca())
  cax = divider.append_axes("right", size="2%", pad=0.05)
  plt.colorbar(cax=cax)
  plt.savefig(os.path.join(output_dir,'measurement_plane_{}_slice_time_average.png'.format(scalar_feild_name)), dpi=300, bbox_inches='tight')

# Measurement Line
for scalar_feild_name in scalar_data_to_plot: 
  plt.figure()
  plt.plot(time_averaged_data['measurement_line_{}'.format(scalar_feild_name)],
           y_centers,
           '-k')
  plt.xlabel('{} {}'.format(title_labels[scalar_feild_name], unit_labels[scalar_feild_name]), fontsize = 20)
  plt.ylabel(r'$y$ [m]', fontsize = 20)
  plt.title('{} on Measurement Line'.format(title_labels[scalar_feild_name]), fontsize = 20)
  plt.savefig(os.path.join(output_dir,'measurement_line_{}_time_average.png'.format(scalar_feild_name)), dpi=300, bbox_inches='tight')

# Center Line
for scalar_feild_name in scalar_data_to_plot: 
  plt.figure()
  plt.plot(z_centers,
           time_averaged_data['center_line_{}'.format(scalar_feild_name)],
           '-k')
  plt.xlabel(r'$x$ [m]', fontsize = 20)
  plt.ylabel('{} {}'.format(title_labels[scalar_feild_name], unit_labels[scalar_feild_name]), fontsize = 20)
  plt.title('{} on Center Line'.format(title_labels[scalar_feild_name]), fontsize = 20)
  plt.savefig(os.path.join(output_dir,'center_line_{}_time_average.png'.format(scalar_feild_name)), dpi=300, bbox_inches='tight')


# Inlet Plane

# Outlet Plane

print('n_timesteps_total = {}'.format(n_timesteps_total))


###############################################################################
# find max and min over the time range 
###############################################################################
stuff_for_plots = {}
for feild_name in (scalar_data_saved + vector_data_saved): 
  stuff_for_plots['mid_x_slice_{}'.format(feild_name)]         = {}
  stuff_for_plots['mid_y_slice_{}'.format(feild_name)]         = {}
  stuff_for_plots['measurement_z_slice_{}'.format(feild_name)] = {}
  stuff_for_plots['inlet_z_slice_{}'.format(feild_name)]       = {}
  stuff_for_plots['outlet_z_slice_{}'.format(feild_name)]      = {}
  stuff_for_plots['{}_spanwise_average'.format(feild_name)]    = {}
  stuff_for_plots['{}_center_line'.format(feild_name)]         = {}
  stuff_for_plots['measurement_line_{}'.format(feild_name)]    = {}

  stuff_for_plots['mid_x_slice_{}'.format(feild_name)]['min']         = np.inf
  stuff_for_plots['mid_y_slice_{}'.format(feild_name)]['min']         = np.inf
  stuff_for_plots['measurement_z_slice_{}'.format(feild_name)]['min'] = np.inf
  stuff_for_plots['inlet_z_slice_{}'.format(feild_name)]['min']       = np.inf
  stuff_for_plots['outlet_z_slice_{}'.format(feild_name)]['min']      = np.inf
  stuff_for_plots['{}_spanwise_average'.format(feild_name)]['min']    = np.inf
  stuff_for_plots['{}_center_line'.format(feild_name)]['min']         = np.inf
  stuff_for_plots['measurement_line_{}'.format(feild_name)]['min']    = np.inf

  stuff_for_plots['mid_x_slice_{}'.format(feild_name)]['max']         = -np.inf
  stuff_for_plots['mid_y_slice_{}'.format(feild_name)]['max']         = -np.inf
  stuff_for_plots['measurement_z_slice_{}'.format(feild_name)]['max'] = -np.inf
  stuff_for_plots['inlet_z_slice_{}'.format(feild_name)]['max']       = -np.inf
  stuff_for_plots['outlet_z_slice_{}'.format(feild_name)]['max']      = -np.inf
  stuff_for_plots['{}_spanwise_average'.format(feild_name)]['max']    = -np.inf
  stuff_for_plots['{}_center_line'.format(feild_name)]['max']         = -np.inf
  stuff_for_plots['measurement_line_{}'.format(feild_name)]['max']    = -np.inf

for filename in hdf_filenames:
  h5_file  = h5py.File(os.path.join(input_dir,filename), "r")
  # Get the timesteps in the file
  for timestep_group in h5_file.keys():
    # skip timestep zero since it is the same as the first timestep of the next restart file
    if timestep_group != 'timestep_0':
      #for feild_name in (scalar_data_saved + vector_data_saved): 
      for feild_name in (scalar_data_to_plot): 
        #for data_name in ['mid_x_slice_', 'measurement_z_slice_', 'inlet_z_slice_', 'outlet_z_slice_', '_spanwise_average']
        if stuff_for_plots['mid_x_slice_{}'.format(feild_name)]['min']         > np.min(np.array(h5_file[timestep_group]['mid_x_slice_{}'.format(feild_name)])):
          stuff_for_plots['mid_x_slice_{}'.format(feild_name)]['min']          = np.min(np.array(h5_file[timestep_group]['mid_x_slice_{}'.format(feild_name)]))
        if stuff_for_plots['mid_y_slice_{}'.format(feild_name)]['min']         > np.min(np.array(h5_file[timestep_group]['mid_y_slice_{}'.format(feild_name)])):
          stuff_for_plots['mid_y_slice_{}'.format(feild_name)]['min']          = np.min(np.array(h5_file[timestep_group]['mid_y_slice_{}'.format(feild_name)]))
        if stuff_for_plots['measurement_z_slice_{}'.format(feild_name)]['min'] > np.min(np.array(h5_file[timestep_group]['measurement_z_slice_{}'.format(feild_name)])):
          stuff_for_plots['measurement_z_slice_{}'.format(feild_name)]['min']  = np.min(np.array(h5_file[timestep_group]['measurement_z_slice_{}'.format(feild_name)]))
        if stuff_for_plots['inlet_z_slice_{}'.format(feild_name)]['min']       > np.min(np.array(h5_file[timestep_group]['inlet_z_slice_{}'.format(feild_name)])):
          stuff_for_plots['inlet_z_slice_{}'.format(feild_name)]['min']        = np.min(np.array(h5_file[timestep_group]['inlet_z_slice_{}'.format(feild_name)]))
        if stuff_for_plots['outlet_z_slice_{}'.format(feild_name)]['min']      > np.min(np.array(h5_file[timestep_group]['outlet_z_slice_{}'.format(feild_name)])):
          stuff_for_plots['outlet_z_slice_{}'.format(feild_name)]['min']       = np.min(np.array(h5_file[timestep_group]['outlet_z_slice_{}'.format(feild_name)]))
        if stuff_for_plots['{}_spanwise_average'.format(feild_name)]['min']    > np.min(np.array(h5_file[timestep_group]['{}_spanwise_average'.format(feild_name)])):
          stuff_for_plots['{}_spanwise_average'.format(feild_name)]['min']     = np.min(np.array(h5_file[timestep_group]['{}_spanwise_average'.format(feild_name)]))
        if stuff_for_plots['{}_center_line'.format(feild_name)]['min']         > np.min(np.array(h5_file[timestep_group]['mid_x_slice_{}'.format(feild_name)][y_mid_idx,:])):
          stuff_for_plots['{}_center_line'.format(feild_name)]['min']          = np.min(np.array(h5_file[timestep_group]['mid_x_slice_{}'.format(feild_name)][y_mid_idx,:]))
        if stuff_for_plots['measurement_line_{}'.format(feild_name)]['min']         > np.min(np.array(h5_file[timestep_group]['measurement_z_slice_{}'.format(feild_name)][x_mid_idx, :])):
          stuff_for_plots['measurement_line_{}'.format(feild_name)]['min']          = np.min(np.array(h5_file[timestep_group]['measurement_z_slice_{}'.format(feild_name)][x_mid_idx, :]))

        if stuff_for_plots['mid_x_slice_{}'.format(feild_name)]['max']         < np.max(np.array(h5_file[timestep_group]['mid_x_slice_{}'.format(feild_name)])):
           stuff_for_plots['mid_x_slice_{}'.format(feild_name)]['max']         = np.max(np.array(h5_file[timestep_group]['mid_x_slice_{}'.format(feild_name)]))
        if stuff_for_plots['mid_y_slice_{}'.format(feild_name)]['max']         < np.max(np.array(h5_file[timestep_group]['mid_y_slice_{}'.format(feild_name)])):
           stuff_for_plots['mid_y_slice_{}'.format(feild_name)]['max']         = np.max(np.array(h5_file[timestep_group]['mid_y_slice_{}'.format(feild_name)]))
        if stuff_for_plots['measurement_z_slice_{}'.format(feild_name)]['max'] < np.max(np.array(h5_file[timestep_group]['measurement_z_slice_{}'.format(feild_name)])):
           stuff_for_plots['measurement_z_slice_{}'.format(feild_name)]['max'] = np.max(np.array(h5_file[timestep_group]['measurement_z_slice_{}'.format(feild_name)]))
        if stuff_for_plots['inlet_z_slice_{}'.format(feild_name)]['max']       < np.max(np.array(h5_file[timestep_group]['inlet_z_slice_{}'.format(feild_name)])):
           stuff_for_plots['inlet_z_slice_{}'.format(feild_name)]['max']       = np.max(np.array(h5_file[timestep_group]['inlet_z_slice_{}'.format(feild_name)]))
        if stuff_for_plots['outlet_z_slice_{}'.format(feild_name)]['max']      < np.max(np.array(h5_file[timestep_group]['outlet_z_slice_{}'.format(feild_name)])):
           stuff_for_plots['outlet_z_slice_{}'.format(feild_name)]['max']      = np.max(np.array(h5_file[timestep_group]['outlet_z_slice_{}'.format(feild_name)]))
        if stuff_for_plots['{}_spanwise_average'.format(feild_name)]['max']    < np.max(np.array(h5_file[timestep_group]['{}_spanwise_average'.format(feild_name)])):
           stuff_for_plots['{}_spanwise_average'.format(feild_name)]['max']    = np.max(np.array(h5_file[timestep_group]['{}_spanwise_average'.format(feild_name)]))
        if stuff_for_plots['{}_center_line'.format(feild_name)]['max']         < np.max(np.array(h5_file[timestep_group]['mid_x_slice_{}'.format(feild_name)][y_mid_idx,:])):
           stuff_for_plots['{}_center_line'.format(feild_name)]['max']         = np.max(np.array(h5_file[timestep_group]['mid_x_slice_{}'.format(feild_name)][y_mid_idx,:]))
        if stuff_for_plots['measurement_line_{}'.format(feild_name)]['max']    < np.max(np.array(h5_file[timestep_group]['measurement_z_slice_{}'.format(feild_name)][x_mid_idx, :])):
          stuff_for_plots['measurement_line_{}'.format(feild_name)]['max']     = np.max(np.array(h5_file[timestep_group]['measurement_z_slice_{}'.format(feild_name)][x_mid_idx, :]))

  h5_file.close() 

###############################################################################
# save animation frames to a folder
###############################################################################

# make folder for movie frames (delete it if it does exist already)
frames_dir = os.path.join(output_dir,'frames')
if not os.path.exists(frames_dir):
  os.makedirs(frames_dir)
  if args.verbose:
    print('Created new directory for frames:')
    print('{}'.format(frames_dir))
else:
  shutil.rmtree(frames_dir)
  os.makedirs(frames_dir)
  if args.verbose:
    print('Clobbered old directory and created new directory for visualization ready data:')
    print('{}'.format(frames_dir))


restart_number = 0
for filename in hdf_filenames:
  h5_file  = h5py.File(os.path.join(input_dir,filename), "r")
  # Get the timesteps in the file
  for timestep_group in h5_file.keys():
    # skip timestep zero since it is the same as the first timestep of the next restart file
    if timestep_group != 'timestep_0':
      restart_number += 1
      t = restart_number*nt_between_restarts*dt
      time_string = r'$t={:7.3f} \ [\mathrm{{s}}]$'.format(t)

      # Mid X Slice
      for scalar_feild_name in scalar_data_to_plot: 
        plt.figure(figsize=(11,9))
        #plt.pcolormesh(np.array(z_points),
        #               np.array(y_points),
        #               h5_file[timestep_group]['mid_x_slice_{}'.format(scalar_feild_name)])
        plt.pcolormesh(np.array(z_points),
                       np.array(y_points),
                       h5_file[timestep_group]['mid_x_slice_{}'.format(scalar_feild_name)],
                       vmin =  stuff_for_plots['mid_x_slice_{}'.format(scalar_feild_name)]['min'],
                       vmax =  stuff_for_plots['mid_x_slice_{}'.format(scalar_feild_name)]['max'])
        ## Add where duct would be
        #rect = Rectangle((0.0, y_duct_start), Lz, duct_width, edgecolor='k', facecolor='none', linestyle='--')
        #ax = plt.gca()
        #ax.add_patch(rect)
        ### Add where averaging area is
        #rect = Rectangle((0.0, y_spanwise_average_start), Lz, L_spanwise_average_window, edgecolor='k', facecolor='none')
        #ax = plt.gca()
        #ax.add_patch(rect)
        # Labels
        plt.axis('scaled')
        plt.xlim([0.0,Lz])
        plt.ylim([0.0,Ly])
        plt.xlabel(r'$x \ [m]$', fontsize = 20)
        plt.ylabel(r'$y \ [m]$', fontsize = 20)
        plt.title('{} {}'.format(title_labels[scalar_feild_name], unit_labels[scalar_feild_name]), fontsize = 20)
        # create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        divider = make_axes_locatable(plt.gca())
        cax = divider.append_axes("right", size="2%", pad=0.05)
        plt.colorbar(cax=cax)
        plt.savefig(os.path.join(frames_dir,'mid_x_slice_{}_{:0>4d}.png'.format(scalar_feild_name,restart_number)), bbox_inches='tight')
        plt.close()

      # Measurement Plane
      for scalar_feild_name in scalar_data_to_plot: 
        plt.figure(figsize=(11,9))
        plt.pcolormesh(y_points,
                       x_points,
                       h5_file[timestep_group]['measurement_z_slice_{}'.format(scalar_feild_name)],
                       vmin =  stuff_for_plots['measurement_z_slice_{}'.format(scalar_feild_name)]['min'],
                       vmax =  stuff_for_plots['measurement_z_slice_{}'.format(scalar_feild_name)]['max'])
      
        # Add where duct would be
        rect = Rectangle((x_duct_start, x_duct_start), duct_width, duct_width, edgecolor='k', facecolor='none', linestyle='--')
        ax = plt.gca()
        ax.add_patch(rect)
      
        # Add where averaging area is
        rect = Rectangle((x_spanwise_average_start, x_spanwise_average_start), L_spanwise_average_window, L_spanwise_average_window, edgecolor='k', facecolor='none')
        ax = plt.gca()
        ax.add_patch(rect)

        plt.axis('scaled')
        plt.xlabel(r'$x \ [cm]$', fontsize = 20)
        plt.ylabel(r'$y \ [cm]$', fontsize = 20)
        plt.title('{} {}'.format(title_labels[scalar_feild_name], unit_labels[scalar_feild_name]), fontsize = 20)
        # create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        divider = make_axes_locatable(plt.gca())
        cax = divider.append_axes("right", size="2%", pad=0.05)
        plt.colorbar(cax=cax)
        plt.savefig(os.path.join(frames_dir,'measurement_plane_{}_{:04d}.png'.format(scalar_feild_name, restart_number)), dpi=300, bbox_inches='tight')
        plt.close()

      # Spanwise average
      for scalar_feild_name in scalar_data_to_plot: 
        plt.figure()
        plt.plot(z_centers,
                 h5_file[timestep_group]['{}_spanwise_average'.format(scalar_feild_name)],
                 '-k',
                 label=time_string )
        plt.plot(z_centers,
                 time_averaged_data['{}_spanwise_average'.format(scalar_feild_name)],
                 '--k',
                 label='time average')
        plt.xlabel(r'$x \ [m]$', fontsize = 20)
        plt.ylabel('{} {}'.format(title_labels[scalar_feild_name], unit_labels[scalar_feild_name]), fontsize = 20)
        plt.ylim([stuff_for_plots['{}_spanwise_average'.format(scalar_feild_name)]['min'], stuff_for_plots['{}_spanwise_average'.format(scalar_feild_name)]['max']])
        plt.title('{} Spanwise Average'.format(title_labels[scalar_feild_name]), fontsize = 20)
        plt.legend(loc='upper left')
        plt.savefig(os.path.join(frames_dir,'spanwise_average_{}_{:0>4d}.png'.format(scalar_feild_name,restart_number)), bbox_inches='tight')
        plt.close()

      # Measurement Line
      for scalar_feild_name in scalar_data_to_plot: 
        plt.figure()
        #plt.plot(h5_file[timestep_group]['measurement_z_slice_{}'.format(scalar_feild_name)][x_mid_idx, y_spanwise_average_start_idx:y_spanwise_average_stop_idx],
        #         y_centers[y_spanwise_average_start_idx:y_spanwise_average_stop_idx],
        #         '-k')
        plt.plot(h5_file[timestep_group]['measurement_z_slice_{}'.format(scalar_feild_name)][x_mid_idx, :],
                 y_centers,
                 '-k',
                 label=time_string)
        plt.plot(time_averaged_data['measurement_line_{}'.format(scalar_feild_name)],
                 y_centers,
                 '--k',
                 label='time average')
        plt.xlim([stuff_for_plots['measurement_line_{}'.format(scalar_feild_name)]['min'],stuff_for_plots['measurement_line_{}'.format(scalar_feild_name)]['max']])
        plt.xlabel('{} {}'.format(title_labels[scalar_feild_name], unit_labels[scalar_feild_name]), fontsize = 20)
        plt.ylabel(r'$y$ [m]', fontsize = 20)
        plt.title('{} on Measurement Line'.format(title_labels[scalar_feild_name]), fontsize = 20)
        plt.savefig(os.path.join(frames_dir,'measurement_line_{}_{:04d}.png'.format(scalar_feild_name, restart_number)), dpi=300, bbox_inches='tight')
        plt.close()

      # Center Line
      for scalar_feild_name in scalar_data_to_plot: 
        plt.figure()
        plt.plot(z_centers,
                 h5_file[timestep_group]['mid_x_slice_{}'.format(scalar_feild_name)][y_mid_idx,:],
                 '-k',
                 label=time_string)
        plt.plot(z_centers,
                 time_averaged_data['center_line_{}'.format(scalar_feild_name)],
                 '--k',
                 label='time average')
        plt.xlabel(r'$x \ [m]$', fontsize = 20)
        plt.ylabel('{} {}'.format(title_labels[scalar_feild_name], unit_labels[scalar_feild_name]), fontsize = 20)
        plt.ylim([stuff_for_plots['measurement_line_{}'.format(scalar_feild_name)]['min'], stuff_for_plots['measurement_line_{}'.format(scalar_feild_name)]['max']])
        plt.title('{} Center Line'.format(title_labels[scalar_feild_name]), fontsize = 20)
        plt.legend(loc='upper left')
        plt.savefig(os.path.join(frames_dir,'center_line_{}_{:0>4d}.png'.format(scalar_feild_name,restart_number)), bbox_inches='tight')
        plt.close()
  
  h5_file.close() 


framerate = 8
for location in ['center_line','measurement_line','spanwise_average', 'measurement_plane', 'mid_x_slice']:
  for scalar_feild_name in scalar_data_to_plot: 
    movie_name = location + '_' + scalar_feild_name + '.mp4'
    frame_basename = os.path.join(frames_dir, location + '_' + scalar_feild_name + '_%04d.png')
    run_command('ffmpeg -r {FRAMERATE} -i {FRAME_BASENAME} -vcodec libx264 -crf 25 {MOVIE_NAME}'.format(FRAMERATE = framerate, FRAME_BASENAME=frame_basename, MOVIE_NAME=movie_name))

# example command to turn images into a movie
#ffmpeg -r 8  -i measurement_line_temperature_%04d.png -vcodec libx264 -crf 25 measurement_line_temperature.mp4
#ffmpeg -r 8  -i measurement_plane_temperature_slice_%04d.png -vcodec libx264 -crf 25 measurement_plane_temperature.mp4
#ffmpeg -r 10 -i temperature_mid_x_slice_%04d.png -vcodec libx264 -crf 25 temperature_mid_x_slice.mp4
#ffmpeg -r 10 -i temperature_spanwise_average_%04d.png -vcodec libx264 -crf 25 temperature_spanwise_average.mp4




################################################################################
## make animation
################################################################################
#
#animation_frames = {}
#animation_figures = {}
#for feild_name in (scalar_data_to_plot + vector_data_to_plot): 
#  animation_frames['mid_x_slice_{}'.format(feild_name)]         = []
#  animation_frames['measurement_z_slice_{}'.format(feild_name)] = []
#  animation_frames['inlet_z_slice_{}'.format(feild_name)]       = []
#  animation_frames['outlet_z_slice_{}'.format(feild_name)]      = []
#  animation_frames['{}_spanwise_average'.format(feild_name)]    = []
#
#  animation_figures['mid_x_slice_{}'.format(feild_name)]         = plt.figure(figsize=(11,9))
#  animation_figures['measurement_z_slice_{}'.format(feild_name)] = plt.figure(figsize=(11,9))
#  animation_figures['inlet_z_slice_{}'.format(feild_name)]       = plt.figure(figsize=(11,9))
#  animation_figures['outlet_z_slice_{}'.format(feild_name)]      = plt.figure(figsize=(11,9))
#
#  animation_figures['{}_spanwise_average'.format(feild_name)]    = plt.figure()
#  plt.xlabel(r'$z \ [m]$', fontsize = 20)
#  plt.ylabel('{} {}'.format(title_labels[scalar_feild_name], unit_labels[scalar_feild_name]), fontsize = 20)
#  plt.title('{} Spanwise Average'.format(title_labels[scalar_feild_name]), fontsize = 20)
#  plt.legend()
#
#timestep = 0
#for filename in hdf_filenames:
#  h5_file  = h5py.File(os.path.join(input_dir,filename), "r")
#  # Get the timesteps in the file
#  for timestep_group in h5_file.keys():
#    # skip timestep zero since it is the same as the first timestep of the next restart file
#    if timestep_group != 'timestep_0':
#      timestep += 1 
#
#      for scalar_feild_name in scalar_data_to_plot: 
#        #plt.figure()
#        #plt.plot(z_centers,
#        #         time_averaged_data['{}_spanwise_average'.format(scalar_feild_name)],
#        #         '-k',
#        #         label='t={}'.format(timestep))
#        #plt.xlabel(r'$z \ [m]$', fontsize = 20)
#        #plt.ylabel('{} {}'.format(title_labels[scalar_feild_name], unit_labels[scalar_feild_name]), fontsize = 20)
#        #plt.title('{} Spanwise Average'.format(title_labels[scalar_feild_name]), fontsize = 20)
#        ##plt.savefig('{:0>4d}_spanwise_average.png'.format(scalar_feild_name), bbox_inches='tight')
#
#        animation_frames['{}_spanwise_average'.format(feild_name)].append(plt.plot(z_centers,
#                               h5_file[timestep_group]['{}_spanwise_average'.format(scalar_feild_name)],
#                               '-k',
#                               label='t={}'.format(timestep)))
#        #plt.savefig('{}_spanwise_average_{:0>4d}.png'.format(scalar_feild_name,timestep), bbox_inches='tight')
#  
#  h5_file.close() 
#
#for scalar_feild_name in scalar_data_to_plot: 
#  animation = animation.ArtistAnimation(animation_figures['{}_spanwise_average'.format(feild_name)],
#                                        animation_frames['{}_spanwise_average'.format(feild_name)], interval=300, repeat_delay=10000, repeat=False, blit=False)
#    

plt.show()
