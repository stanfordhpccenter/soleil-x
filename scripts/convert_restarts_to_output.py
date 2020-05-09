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
                    help='directory with all the simulation output')
parser.add_argument('--outdir', nargs='?', const='default', default='default',
                    help='directory where the simulation output were go')
parser.add_argument('-v', '--verbose',
                    action='store_true',
                    help='verbose output')
parser.add_argument('--debug',
                    action='store_true',
                    help='run in debug mode')
args = parser.parse_args()


# Turn pathnames into absolute paths
sample_base_dir = os.path.abspath(args.basedir)

# hack way to set default of outdir based on the sampledir argument
out_base_dir = args.outdir
if args.outdir == 'default':
  out_base_dir = os.path.join(sample_base_dir,'viz_ready_data')
out_base_dir = os.path.abspath(out_base_dir)

if args.verbose:
  print('################################################################################')
  print('#                              Input Summary                                   #')
  print('################################################################################')
  print('sample base directory: {}'.format(sample_base_dir))
  print('output base directory: {}'.format(out_base_dir))
  print('')

################################################################################
#                              Error Checking                                  #
################################################################################

if not os.path.exists(sample_base_dir):
  print('################################################################################')
  print('#                                 ERROR                                        #')
  print('################################################################################')
  print('Error: the provided sample base directory {} does not exist'.format(sample_base_dir))
  sys.exit()

#sample0_base_dir = os.path.join(sample_base_dir,'sample0')
#if not os.path.exists(sample0_base_dir):
#  print('################################################################################')
#  print('#                                 ERROR                                        #')
#  print('################################################################################')
#  print('Error: the provided sample0 base directory {} does not exist'.format(sample0_base_dir))
#  sys.exit()
#
#sample1_base_dir = os.path.join(sample_base_dir,'sample1')
#if not os.path.exists(sample1_base_dir):
#  print('################################################################################')
#  print('#                                 ERROR                                        #')
#  print('################################################################################')
#  print('Error: the provided sample1 base directory {} does not exist'.format(sample1_base_dir))
#  sys.exit()

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

viz_script = os.path.expandvars('$SOLEIL_DIR/scripts/convert_output_for_viz.py')
if not os.path.isfile(viz_script):
  print('################################################################################')
  print('#                                 ERROR                                        #')
  print('################################################################################')
  print('file: {} does not exist'.format())
  sys.exit()

print('##############################################################################')
print('#                Set up directory for visualization ready data files         #')
print('##############################################################################')

if args.debug:
    print('Would create new directory for visualization ready data:')
    print('{}'.format(out_base_dir))
else:
  if not os.path.exists(out_base_dir):
    os.makedirs(out_base_dir)
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
# Get list of sample direcories
################################################################################
#for root, dirs, files in os.walk(sample_base_dir):
#  #print('root = {}'.format(root))
#  #print('dirs = {}'.format(dirs))
#  #print('files = {}'.format(files))
#  #print('')
#
#  path, dirname = os.path.split(root)
#  if 'sample' in dirname:
#    print(dirname)

#sample_dirs = [directory for directory in os.listdir(sample_base_dir) if 'sample' in directory]
#print(sample_dirs)

# Get list of sample directories
sample_dirs = [os.path.join(sample_base_dir,directory) for directory in os.listdir(sample_base_dir) if 'sample' in directory]
sample_dirs.sort()
#print(sample_dirs)


# See if need to merge data
for sample_dir in sample_dirs:
  merge_data = False

  path, sample_dir_basename = os.path.split(sample_dir)
  out_dir = os.path.join(out_base_dir,sample_dir_basename)

  fluid_dirs = [os.path.join(sample_dir,directory) for directory in os.listdir(sample_dir) if 'fluid' in directory and os.path.isdir(os.path.join(sample_dir,directory))]
  fluid_dirs.sort()
  #print(fluid_dirs)

  for fluid_dir in fluid_dirs:
    fluid_restarts = [os.path.join(fluid_dir,restart_file) for restart_file in os.listdir(fluid_dir) if os.path.isfile(os.path.join(fluid_dir,restart_file))]
    #print(fluid_restarts)

    if len(fluid_restarts) > 1:
      #print(len(fluid_restarts))
      merge_data = True
      break

  if merge_data == True:
    sample_dir_merged_data = os.path.join(sample_dir,'merged_data')
    command = 'python {} --sampledir {} --outdir {}'.format(
                                                     merge_data_script,
                                                     sample_dir,
                                                     sample_dir_merged_data)
    # Merge Data
    if args.debug:
      print('Would run command:')
      print(command)
    else:
      try:
        subprocess.check_output(command, shell=True)
        print('Running command:')
        print('{}'.format(command))
      except subprocess.CalledProcessError as e:
        print('Failed command with output:')
        print('{}'.format(command))
        print(e.output)
        sys.exit()
      else:
        print('Successfully ran command:')
        print('{}'.format(command))
    print('')


    # Convert for viz
    command = 'python {} --sampledir {} --outdir {}'.format(
                                                     viz_script,
                                                     sample_dir_merged_data,
                                                     out_dir)
    if args.debug:
      print('Would run command:')
      print(command)
    else:
      try:
        subprocess.check_output(command, shell=True)
        print('Running command:')
        print('{}'.format(command))
      except subprocess.CalledProcessError as e:
        print('Failed command with output:')
        print('{}'.format(command))
        print(e.output)
        sys.exit()
      else:
        print('Successfully ran command:')
        print('{}'.format(command))
    print('')



  else : 
    command = 'python {} --sampledir {} --outdir {}'.format(
                                                     viz_script,
                                                     sample_dir,
                                                     out_dir)
    if args.debug:
      print('Would run command:')
      print(command)
    else:
      try:
        subprocess.check_output(command, shell=True)
        print('Running command:')
        print('{}'.format(command))
      except subprocess.CalledProcessError as e:
        print('Failed command with output:')
        print('{}'.format(command))
        print(e.output)
        sys.exit()
      else:
        print('Successfully ran command:')
        print('{}'.format(command))
    print('')



  #particle_dirs = [os.path.join(sample_dir,directory) for directory in os.listdir(sample_dir) if 'particles' in directory and os.path.isdir(os.path.join(sample_dir,directory))]
  #particle_dirs.sort()
  #print(particle_dirs)
  #print('')

#################################################################################
## Create top level file 
#################################################################################
#
#sample0_out_dir = os.path.join(out_base_dir, 'sample0')
#sample1_out_dir = os.path.join(out_base_dir, 'sample1')
#
#for sample_dir, out_dir in zip([sample0_base_dir, sample1_base_dir], [sample0_out_dir, sample1_out_dir]):
#  
#  command = 'python {} --sampledir {} --outdir {}'.format(
#                                                   viz_script,
#                                                   sample_dir,
#                                                   out_dir)
#  if args.debug:
#    print('Would run command:')
#    print(command)
#  else:
#    try:
#      subprocess.check_output(command, shell=True)
#      print('Running command:')
#      print('{}'.format(command))
#    except subprocess.CalledProcessError as e:
#      print('Failed command with output:')
#      print('{}'.format(command))
#      print(e.output)
#      sys.exit()
#    else:
#      print('Successfully ran command:')
#      print('{}'.format(command))
#  print('')
