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
print('                     Set up directory for viz ready data files')
print('##############################################################################')
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

print('##############################################################################')
print('                          Generate fluid viz files ')
print('##############################################################################')

viz_fluid_command = 'python {} -s {} {} {}'.format(os.path.expandvars('$SOLEIL_DIR/scripts/viz_fluid.py'), 
                                                   args.section,
                                                   args.json_file,
                                                   os.path.join(sample_dir,'fluid_iter*/*.hdf'))
try:
  subprocess.call(viz_fluid_command, shell=True)
except OSError:
  print("Failed command: {}".format(viz_fluid_command))
  sys.exit()
else: 
  print("Successfully ran command: {}".format(viz_fluid_command))

# Move the generated files to the output directory
try:  
    subprocess.call('mv out_fluid.xmf out_fluid*.hdf {}'.format(out_dir), shell=True)
except OSError:  
    print("Failed to move fluid xmf and hdf files to: {}".format(out_dir))
    sys.exit()
else:  
    print("Successfully moved fluid xmf and hdf files to: {}".format(out_dir))

print('')

print('##############################################################################')
print('                       Generate particle viz files ')
print('##############################################################################')

viz_particles_command = 'python {} {}'.format(os.path.expandvars('$SOLEIL_DIR/scripts/viz_particles.py'),
                                              os.path.join(sample_dir,'particles_iter*/*.hdf'))
try:
  subprocess.call(viz_particles_command, shell=True)
except OSError:
  print("Failed command: {}".format(viz_particles_command))
  sys.exit()
else: 
  print("Successfully ran command: {}".format(viz_particles_command))

# Move the generated files to the output directory
try:  
    subprocess.call('mv out_particles.xmf out_particles*.hdf {}'.format(out_dir), shell=True)
except OSError:  
    print("Failed to move particles xmf and hdf files to {}".format(out_dir))
    sys.exit()
else:  
    print("Successfully moved particles xmf and hdf files to {}".format(out_dir))

