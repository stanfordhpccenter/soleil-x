import argparse
import os
import subprocess
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--sampledir', nargs='?', const='.', default='.',
                    help='directory with all the simulation output')
args = parser.parse_args()

sample_dir = args.sampledir
merged_dir = os.path.join(sample_dir,'merged_data')

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
print('                        Generate and move merged data files ')
print('##############################################################################')

for content in os.listdir(sample_dir):
  
  if os.path.isdir(os.path.join(sample_dir, content)):

    origional_hdf_data_path = os.path.join(sample_dir, content)
    merged_hdf_data_path = os.path.join(merged_dir, content)

    merge_command = ''
    mv_command = ''
    if "fluid" in os.path.basename(os.path.normpath(origional_hdf_data_path)):
      merge_command = 'python {} --output_filename {} {}'.format(os.path.expandvars('$SOLEIL_DIR/scripts/merge_fluid.py'),
                                                                 'merged_fluid.hdf',
                                                                 os.path.join(origional_hdf_data_path,'*.hdf'))
      if not os.path.exists(merged_hdf_data_path):
        os.makedirs(merged_hdf_data_path)

      mv_command = 'mv ./merged_fluid.hdf {}'.format(os.path.join(merged_hdf_data_path,'merged_fluid.hdf'))
    elif "particle" in os.path.basename(os.path.normpath(origional_hdf_data_path)):
      merge_command = 'python {} --output_filename {} {}'.format(os.path.expandvars('$SOLEIL_DIR/scripts/merge_particles.py'),
                                                                 'merged_particles.hdf',
                                                                 os.path.join(origional_hdf_data_path,'*.hdf'))

      if not os.path.exists(merged_hdf_data_path):
        os.makedirs(merged_hdf_data_path)

      mv_command = 'mv ./merged_particles.hdf {}'.format(os.path.join(merged_hdf_data_path,'merged_particles.hdf'))


    if "fluid" in os.path.basename(os.path.normpath(origional_hdf_data_path)) or "particle" in os.path.basename(os.path.normpath(origional_hdf_data_path)):
      try:
        subprocess.check_output(merge_command, shell=True)
      except subprocess.CalledProcessError as e:
        print("Failed command: {}".format(merge_command))
        print(e.output)
        sys.exit()
      else: 
        print("Successfully merged hdf data in: {}".format(origional_hdf_data_path))

      # Move the generated files to the output directory
      try:  
        subprocess.check_output(mv_command, shell=True)
      except subprocess.CalledProcessError as e:
        print("Failed command: {}".format(mv_command))
        print(e.output)
        sys.exit()
      else:  
        print("Successfully moved merged fluid hdf file to: {}".format(merged_hdf_data_path))
    else:
       print('Skipping merging hdf files in the directory "{}" because it does not have "fluid" or "particle" in the basename'.format(origional_hdf_data_path))
    
    print('')


