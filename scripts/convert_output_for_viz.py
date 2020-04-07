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
#parser.add_argument('json_file',
#                    help='original simulation configuration file')
#parser.add_argument('-s', '--section', choices=['1','2'], 
#                    help='which section to visualize (if multi-section sim)')
parser.add_argument('--sampledir', nargs='?', const='.', default='.',
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
sample_dir = os.path.abspath(args.sampledir)

# hack way to set default of outdir based on the sampledir argument
out_dir = args.outdir
if args.outdir == 'default':
  out_dir = os.path.join(sample_dir,'vis_ready_data')
out_dir = os.path.abspath(out_dir)

if args.verbose:
  print('################################################################################')
  print('#                              Input Summary                                   #')
  print('################################################################################')
  #print('case file: {}'.format(args.json_file))
  #print('config in case file (section): {}'.format(args.section))
  print('sample directory: {}'.format(sample_dir))
  print('output directory: {}'.format(out_dir))
  print('')

if not os.path.exists(sample_dir):
  print('################################################################################')
  print('#                                 ERROR                                        #')
  print('################################################################################')
  print('Error: the provided sample directory {} does not exist'.format(sample_dir))
  sys.exit()

if 'SOLEIL_DIR' not in os.environ:
  print('################################################################################')
  print('#                                 ERROR                                        #')
  print('################################################################################')
  print('Environment variable SOLEIL_DIR not set')
  print('This should point to the directory where you installed Soliel-X')
  sys.exit()

fluid_viz_script = os.path.expandvars('$SOLEIL_DIR/scripts/viz_fluid_nonuniform.py')
if not os.path.isfile(fluid_viz_script):
  print('################################################################################')
  print('#                                 ERROR                                        #')
  print('################################################################################')
  print('file: {} does not exist'.format())
  sys.exit()

particles_viz_script = os.path.expandvars('$SOLEIL_DIR/scripts/viz_particles.py')
if not os.path.isfile(particles_viz_script):
  print('################################################################################')
  print('#                                 ERROR                                        #')
  print('################################################################################')
  print('file: {} does not exist'.format())
  sys.exit()


################################################################################
# Create top level file 
################################################################################
print('##############################################################################')
print('#                Set up directory for visualization ready data files         #')
print('##############################################################################')
if args.debug:
    print('Would create new directory for visualization ready data:')
    print('{}'.format(out_dir))
else:
  if not os.path.exists(out_dir):
    os.makedirs(out_dir)
    print('Created new directory for visualization ready data:')
    print('{}'.format(out_dir))
    os.chdir(out_dir)
  else:
    print('################################################################################')
    print('#                                 ERROR                                        #')
    print('################################################################################')
    print('Directory for visualization ready data already exists:')
    print('{}'.format(out_dir))
    print('Delete it and rerun this script')
    print('')
    sys.exit()
print('')


print('##############################################################################')
print('#                          Generate fluid viz files                          #')
print('##############################################################################')

fluid_hdf_files = sorted(glob.glob(os.path.join(sample_dir,'fluid_iter*/*.hdf')))
if args.verbose:
  print('Found fluid files:')
  for filename in fluid_hdf_files:
    print(filename)
  print('')

#viz_fluid_command = 'python {} -s {} {} {}'.format(fluid_viz_script,
#                                                   args.section,
#                                                   args.json_file,
#                                                   ' '.join(fluid_hdf_files))

viz_fluid_command = 'python {} {}'.format(fluid_viz_script,
                                          ' '.join(fluid_hdf_files))
if args.debug:
  print('Would run command for fluid:')
  print(viz_fluid_command)
else:
  try:
    subprocess.check_output(viz_fluid_command, shell=True)
  except subprocess.CalledProcessError as e:
    print('Failed command:')
    print('{}'.format(viz_fluid_command))
    print(e.output)
    sys.exit()
  else:
    print('Successfully ran command:')
    print('{}'.format(viz_fluid_command))
print('')

print('##############################################################################')
print('#                       Generate particle viz files                          #')
print('##############################################################################')

particle_hdf_files = sorted(glob.glob(os.path.join(sample_dir,'particles_iter*/*.hdf')))
print('Found particle files:')
for filename in particle_hdf_files:
  print(filename)
print('')

viz_particles_command = 'python {} {}'.format(particles_viz_script,
                                              ' '.join(particle_hdf_files))
if args.debug:
  print('Would run command for particles:')
  print(viz_particles_command)
else:
  try:
    subprocess.check_output(viz_particles_command, shell=True)
  except subprocess.CalledProcessError as e:
    print('Failed command:')
    print('{}'.format(viz_particles_command))
    print(e.output)
    sys.exit()
  else:
    print('Successfully ran command:')
    print('{}'.format(viz_particles_command))
