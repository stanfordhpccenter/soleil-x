import os
import shutil
import sys
import subprocess

rootdir = '.'
outdir  = os.path.join(rootdir,'xmf_output')


# if there is already an output directory then delete it
if os.path.isdir(outdir):
  try:  
      #os.rmdir(outdir)
      shutil.rmtree(outdir)
  except OSError:  
      print("Failed to delete directory: {}".format(outdir))
      sys.exit()
  else:  
      print("Successfully deleted directory: {}".format(outdir))

# Create the direcotry for the output
try:  
    os.mkdir(outdir)
except OSError:  
    print("Failed to create directory: {}".format(outdir))
    sys.exit()
else:  
    print("Successfully created directory: {}".format(outdir))


dirs = os.listdir(rootdir)
dirs.sort()

#break_ = False
#for type_of_data in ['fluid', 'particles']:
for type_of_data in ['particles']:
  dump_number=0
  for dir_name in dirs:
    if dir_name.startswith(type_of_data):

      sub_dir_name = os.path.join(rootdir, dir_name)

      #print(sub_dir_name)

      # assumes one file name that I want is the only file in the directory
      original_hdf_file_name = os.path.join(sub_dir_name, os.listdir(sub_dir_name)[0])
      python_script = os.path.expandvars('$SOLEIL_DIR/scripts/viz_{}.py'.format(type_of_data))
      command = 'python {} {}'.format(python_script, original_hdf_file_name)

      try:
        subprocess.call(command, shell=True)
      except OSError:
        print("command failed: {}".format(command))
        sys.exit()
      else: 
        print("Successfully ran command: {}".format(command))

      
      output_file_basename = dir_name

      # rename and move the xmf output of the viz_particles.py script
      viz_output_xmf_filename = os.path.join(rootdir,'out.xmf')
      new_xmf_filename = os.path.join(outdir,'{}_timestep_{:04}.xmf'.format(type_of_data,dump_number))
      try:  
        os.rename(viz_output_xmf_filename, new_xmf_filename)
      except OSError:  
        print("Failed to move file: {} -> {}".format(viz_output_xmf_filename,new_xmf_filename))
        sys.exit()
      else:  
        print("Successfully moved file: {} -> {}".format(viz_output_xmf_filename,new_xmf_filename))


      # rename and move the hdf output of the viz_particles.py script
      viz_output_hdf_filename = os.path.join(rootdir,'out.hdf')
      new_hdf_filename = os.path.join(outdir,output_file_basename+'.hdf')
      try:
        os.rename(viz_output_hdf_filename , new_hdf_filename)
      except OSError:  
        print("Failed to move file: {} -> {}".format(viz_output_hdf_filename,new_hdf_filename))
        sys.exit()
      else:  
        print("Successfully moved file: {} -> {}".format(viz_output_hdf_filename,new_hdf_filename))

      # Change the xmf file contents to point to the correct hdf file
      try: 
        with open(new_xmf_filename, 'r') as myfile:
          xmf_text=myfile.read().replace('out.hdf', output_file_basename+'.hdf')
        with open(new_xmf_filename, 'w') as myfile:
          myfile.write(xmf_text)
      except:
        print("Failed to edit file: {}".format(new_xmf_filename))
        sys.exit()
      else:
        print("Successfully edited file: {}".format(new_xmf_filename))

      dump_number+=1

#      if break_:
#        break
#      else:
#        break_ = True


