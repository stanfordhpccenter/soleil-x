#!/usr/bin/env python3

import os
import shutil
import sys
import subprocess
import argparse
import glob

################################################################################
#                               Parse Input                                    #
################################################################################
parser = argparse.ArgumentParser()
parser.add_argument('--basedir', nargs='?', const='.', default='.',
                    help='directory that contains the sample directories. basedir usually has subdirectories: sample0, sample1, ...')
parser.add_argument('--outdir', nargs='?', const='default', default='default',
                    help='directory where output is written')
parser.add_argument('-v', '--verbose',
                    action='store_true',
                    help='verbose output')
parser.add_argument('--zip',
                    action='store_true',
                    help='create zip file')
parser.add_argument('--debug',
                    action='store_true',
                    help='run in debug mode')
args = parser.parse_args()

# hack way to set default of outdir based on the basedir argument
out_base_dir = args.outdir
if args.outdir == 'default':
  out_base_dir = os.path.join(args.basedir,'viz_ready_data')

# Turn pathnames into absolute paths
# doing this because later os.walk() can have trouble with relative path names.
base_dir = os.path.abspath(args.basedir)
out_base_dir = os.path.abspath(out_base_dir)

if args.verbose:
  print('################################################################################')
  print('#                              Input Summary                                   #')
  print('################################################################################')
  print('input  base directory: {}'.format(base_dir))
  print('output base directory: {}'.format(out_base_dir))
  print('')

################################################################################
#                              Error Checking                                  #
################################################################################

if not os.path.exists(base_dir):
  print('################################################################################')
  print('#                                 ERROR                                        #')
  print('################################################################################')
  print('Error: the provided sample base directory {} does not exist'.format(base_dir))
  sys.exit()

if base_dir == out_base_dir:
  print('################################################################################')
  print('#                                 ERROR                                        #')
  print('################################################################################')
  print('Error: arguments basedir outdir cannot be the same')
  sys.exit()

if 'SOLEIL_DIR' not in os.environ:
  print('################################################################################')
  print('#                                 ERROR                                        #')
  print('################################################################################')
  print('Environment variable SOLEIL_DIR not set')
  print('This should point to the directory where you installed Soliel-X')
  sys.exit()

merge_data_script = os.path.expandvars('$SOLEIL_DIR/scripts/merge_data.py')
if not os.path.isfile(merge_data_script):
  print('################################################################################')
  print('#                                 ERROR                                        #')
  print('################################################################################')
  print('file: {} does not exist'.format(merge_data_script))
  sys.exit()

viz_script = os.path.expandvars('$SOLEIL_DIR/scripts/viz_sample.py')
if not os.path.isfile(viz_script):
  print('################################################################################')
  print('#                                 ERROR                                        #')
  print('################################################################################')
  print('file: {} does not exist'.format())
  sys.exit()

################################################################################
# Get list of sample direcories
################################################################################
sample_dirs = [os.path.join(base_dir,directory) for directory in os.listdir(base_dir) if 'sample' in directory]
sample_dirs.sort()

if not sample_dirs:
  print('################################################################################')
  print('#                                 ERROR                                        #')
  print('################################################################################')
  print('Input base directory ( {} ) did not contain any sample directories. This script looks in the input base directory for sub-directores with the substring "sample" in the name. None were found.'.format(base_dir))
  print('This script assumes a input base directory format of:')
  print('input_base_dir')
  print('|-- sample0')
  print('|-- sample1')
  print('|--   .')
  print('|--   .')
  print('|--   .')
  print('|-- sampleN')
  sys.exit()

if args.verbose:
  print('##############################################################################')
  print('#                Find sample directories in input base directory             #')
  print('##############################################################################')
  print('Found sample directories:')
  [print(sample_dir) for sample_dir in sample_dirs]
  print('')

################################################################################
#                        Create Directory for Output                           #
################################################################################
if args.verbose:
  print('##############################################################################')
  print('#                Make directory for output paraview friendly output          #')
  print('##############################################################################')

if args.debug:
  print('Would create new directory for visualization ready data:')
  print('{}'.format(out_base_dir))
else:
  if not os.path.exists(out_base_dir):
    os.makedirs(out_base_dir)
    if args.verbose:
      print('Created new directory for visualization ready data:')
      print('{}'.format(out_base_dir))
  else:
    print('################################################################################')
    print('#                                 ERROR                                        #')
    print('################################################################################')
    print('Directory for visualization ready data already exists:')
    print('{}'.format(out_base_dir))
    print('I don\'t want to clobber whatever is there.')
    print('Delete it and re-run and this script')
    print('')
    sys.exit()
print('')

################################################################################
# Process data in sample directories
################################################################################
if args.verbose:
  print('##############################################################################')
  print('#                             Process Data                                   #')
  print('##############################################################################')

for sample_dir in sample_dirs:

  path, sample_dir_basename = os.path.split(sample_dir)
  out_dir = os.path.join(out_base_dir,sample_dir_basename)

  sample_dir_to_viz = sample_dir

  fluid_dirs = [os.path.join(sample_dir,fluid_directory) for fluid_directory in os.listdir(sample_dir) if 'fluid' in fluid_directory and os.path.isdir(os.path.join(sample_dir,fluid_directory))]
  fluid_dirs.sort()

  particle_dirs = [os.path.join(sample_dir,particle_directory) for particle_directory in os.listdir(sample_dir) if 'particles' in particle_directory and os.path.isdir(os.path.join(sample_dir,particle_directory))]
  particle_dirs.sort()

  if args.verbose:
    print('In sample directory: {}'.format(sample_dir))
    print('Found fluid directories:'.format(sample_dir))
    [print(fluid_dir) for fluid_dir in fluid_dirs]
    print('Found particle directories:'.format(sample_dir))
    [print(particle_dir) for particle_dir in particle_dirs]
    print('')

  # Determine if the data needs to be merged 
  merge_data = False
  for fluid_dir in fluid_dirs:
    fluid_restarts = [os.path.join(fluid_dir,restart_file) for restart_file in os.listdir(fluid_dir) if os.path.isfile(os.path.join(fluid_dir,restart_file))]
    # If more than one file in the restart directories we will need to merge them
    if len(fluid_restarts) > 1:
      merge_data = True
      break

  # Merge the data into one hdf5 file if needed
  if merge_data == True:
    if args.verbose:
      print('################################################################################')
      print('#                               Merge Data                                     #')
      print('################################################################################')
      print('merging data in sample dir: {}'.format(sample_dir))

    sample_dir_merged_data = os.path.join(sample_dir,'merged_data')
    command = 'python {} --sampledir {} --outdir {}'.format(
                                                     merge_data_script,
                                                     sample_dir,
                                                     sample_dir_merged_data)
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

    # Change input dir for viz script to point to merged data
    sample_dir_to_viz = sample_dir_merged_data 

  # Convert the hdf5 files to a more paraview friendly format
  if args.verbose:
    print('################################################################################')
    print('#                               Convert Data                                   #')
    print('################################################################################')
    print('converting data in sample dir: {}'.format(sample_dir))
  # Convert data to output
  command = 'python {} --sampledir {} --outdir {}'.format(
                                                   viz_script,
                                                   sample_dir_to_viz,
                                                   out_dir)
  if args.verbose:
    command = command + ' --verbose'
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


if args.zip:
  if args.verbose:
    print('##############################################################################')
    print('#                        Zip Visualization Ready Data                        #')
    print('##############################################################################')

  zip_filename = out_base_dir + '.tar.gz'
  # Convert data to output
  command = 'tar -czvf {} {}'.format(zip_filename,
                                                  out_dir)

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
