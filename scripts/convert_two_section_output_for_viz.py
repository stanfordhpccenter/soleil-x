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
  #print('case file: {}'.format(args.json_file))
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

if 'SOLEIL_DIR' not in os.environ:
  print('################################################################################')
  print('#                                 ERROR                                        #')
  print('################################################################################')
  print('Environment variable SOLEIL_DIR not set')
  print('This should point to the directory where you installed Soliel-X')
  sys.exit()

viz_script = os.path.expandvars('$SOLEIL_DIR/scripts/convert_output_for_viz.py')
if not os.path.isfile(viz_script):
  print('################################################################################')
  print('#                                 ERROR                                        #')
  print('################################################################################')
  print('file: {} does not exist'.format())
  sys.exit()

sample0_base_dir = os.path.join(sample_base_dir,'sample0')
if not os.path.exists(sample0_base_dir):
  print('################################################################################')
  print('#                                 ERROR                                        #')
  print('################################################################################')
  print('Error: the provided sample0 base directory {} does not exist'.format(sample0_base_dir))
  sys.exit()

sample1_base_dir = os.path.join(sample_base_dir,'sample1')
if not os.path.exists(sample1_base_dir):
  print('################################################################################')
  print('#                                 ERROR                                        #')
  print('################################################################################')
  print('Error: the provided sample1 base directory {} does not exist'.format(sample1_base_dir))
  sys.exit()

################################################################################
# Create top level file 
################################################################################
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
    print('Delete it and rerun this script')
    print('')
    sys.exit()
print('')

sample0_out_dir = os.path.join(out_base_dir, 'sample0')
sample1_out_dir = os.path.join(out_base_dir, 'sample1')

#for section_number, (sample_dir, out_dir) in enumerate(zip([sample0_base_dir, sample1_base_dir], [sample0_out_dir, sample0_out_dir])):
for sample_dir, out_dir in zip([sample0_base_dir, sample1_base_dir], [sample0_out_dir, sample1_out_dir]):
  
  #print('##############################################################################')
  #print('#                          Sample {}                                         #'.format(section_number))
  #print('##############################################################################')
  
  #command = 'python {} --section {} --sampledir {} --outdir {} {}'.format(
  #                                                   viz_script,
  #                                                   section_number+1,
  #                                                   sample_dir,
  #                                                   out_dir,
  #                                                   args.json_file)

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
    except subprocess.CalledProcessError as e:
      print('Failed command:')
      print('{}'.format(command))
      print(e.output)
      sys.exit()
    else:
      print('Successfully ran command:')
      print('{}'.format(command))
  print('')
