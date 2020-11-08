import argparse
import os
import subprocess
import sys

parser = argparse.ArgumentParser()
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

sample_dir = args.sampledir

# hack way to set default of outdir based on the sampledir argument
out_dir = args.outdir
if args.outdir == 'default':
  outdir = os.path.join(sampledir,'merged_data')
merged_dir = os.path.abspath(out_dir)


print('##############################################################################')
print('                     Set up directory for merged data files')
print('##############################################################################')
if os.path.exists(sample_dir):
  if not os.path.exists(merged_dir):
    os.makedirs(merged_dir)
    print("Created new directory for merged data: {}".format(merged_dir))
  else:
    print("Directory for merged data already exists: {}".format(merged_dir))
else:
  print('Error: the provided sample directory {} does not exist'.format(sample_dir))
  sys.exit()
print('')

print('##############################################################################')
print('                           Generate merged data files ')
print('##############################################################################')

for content in os.listdir(sample_dir):
  
  if os.path.isdir(os.path.join(sample_dir, content)):

    origional_hdf_data_path = os.path.join(sample_dir, content)
    merged_hdf_data_path = os.path.join(merged_dir, content)

    merge_command = ''
    if "fluid" in os.path.basename(os.path.normpath(origional_hdf_data_path)):
      if not os.path.exists(merged_hdf_data_path):
        os.makedirs(merged_hdf_data_path)

      merge_command = 'python {} --output_filename {} {}'.format(os.path.expandvars('$SOLEIL_DIR/scripts/merge_fluid.py'),
                                                                 os.path.join(merged_hdf_data_path,'merged_fluid.hdf'),
                                                                 os.path.join(origional_hdf_data_path,'*.hdf'))
    elif "particle" in os.path.basename(os.path.normpath(origional_hdf_data_path)):
      if not os.path.exists(merged_hdf_data_path):
        os.makedirs(merged_hdf_data_path)

      merge_command = 'python {} --output_filename {} {}'.format(os.path.expandvars('$SOLEIL_DIR/scripts/merge_particles.py'),
                                                                 os.path.join(merged_hdf_data_path,'merged_particles.hdf'),
                                                                 os.path.join(origional_hdf_data_path,'*.hdf'))


    if "fluid" in os.path.basename(os.path.normpath(origional_hdf_data_path)) or "particle" in os.path.basename(os.path.normpath(origional_hdf_data_path)):
      try:
        subprocess.check_output(merge_command, shell=True)
      except subprocess.CalledProcessError as e:
        print("Failed command: {}".format(merge_command))
        print(e.output)
        sys.exit()
      else: 
        print("Successfully merged hdf data in: {}".format(origional_hdf_data_path))
    else:
       print('Skipping merging hdf files in the directory "{}" because it does not have "fluid" or "particle" in the basename'.format(origional_hdf_data_path))
    
    print('')
