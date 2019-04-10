import os
import shutil
import sys
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--section', choices=['1','2'],
                    help='which section to visualize (if multi-section sim)')
parser.add_argument('json_file',
                    help='original simulation configuration file')
parser.add_argument('--sampledir', nargs='?', const='.', default='.',
                    help='directory with all the simulation output')
args = parser.parse_args()

sample_dir = args.sampledir
out_dir    = os.path.join(sample_dir,'viz_ready_data')

print('##############################################################################')
print('                Set up directory for visualization ready data files')
print('##############################################################################')

if os.path.exists(sample_dir):
  if not os.path.exists(out_dir):
    os.makedirs(out_dir)
    print("Created new directory for visualization ready data: {}".format(out_dir))
  else:
    print("Directory for visualization ready data already exists: {}".format(out_dir))
else:
  print('Error: the provided sample directory {} does not exist'.format(sample_dir))
  sys.exit()
print('')

print('##############################################################################')
print('                          Generate fluid viz files ')
print('##############################################################################')

viz_fluid_command = 'python {} -s {} {} {}'.format(os.path.expandvars('$SOLEIL_DIR/scripts/viz_fluid.py'), 
                                                   args.section,
                                                   args.json_file,
                                                   os.path.join(sample_dir,'fluid*/*.hdf'))
try:
  subprocess.check_output(viz_fluid_command, shell=True)
except subprocess.CalledProcessError as e:
  print("Failed command: {}".format(viz_fluid_command))
  print(e.output)
  sys.exit()
else:  
  print("Successfully ran command: {}".format(viz_fluid_command))
  #print("Successfully moved merged fluid hdf file to: {}".format(merged_hdf_data_path))

# Move the generated files to the output directory
try:  
    subprocess.check_output('mv out_fluid.xmf out_fluid*.hdf {}'.format(out_dir), shell=True)
except subprocess.CalledProcessError as e:
  #print("Failed command: {}".format(viz_fluid_command))
  print("Failed to move fluid xmf and hdf files to: {}".format(out_dir))
  print(e.output)
  sys.exit()
else:  
  print("Successfully moved fluid xmf and hdf files to: {}".format(out_dir))

print('')

print('##############################################################################')
print('                       Generate particle viz files ')
print('##############################################################################')

viz_particles_command = 'python {} {}'.format(os.path.expandvars('$SOLEIL_DIR/scripts/viz_particles.py'),
                                              os.path.join(sample_dir,'particles*/*.hdf'))
try:
  subprocess.check_output(viz_particles_command, shell=True)
except subprocess.CalledProcessError as e:
  print("Failed command: {}".format(viz_particles_command))
  print(e.output)
  sys.exit()
else:  
  print("Successfully ran command: {}".format(viz_particles_command))
  #print("Successfully moved merged fluid hdf file to: {}".format(merged_hdf_data_path))

# Move the generated files to the output directory
try:  
    subprocess.check_output('mv out_particles.xmf out_particles*.hdf {}'.format(out_dir), shell=True)
except subprocess.CalledProcessError as e:
  #print("Failed command: {}".format(mv_particle_command))
  print("Failed to move particles xmf and hdf files to {}".format(out_dir))
  print(e.output)
  sys.exit()
else:  
  print("Successfully moved particles xmf and hdf files to {}".format(out_dir))

