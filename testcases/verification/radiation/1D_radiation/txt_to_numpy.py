import sys 
import numpy as np

# --------------------------------------------------------------------------- #
#                            Read Command Line Input                          #
# --------------------------------------------------------------------------- #
if len(sys.argv) == 2:
    # Directory containing the out data of mcrt.rg
    data_dir = sys.argv[1]

    # DOM results numpy file (will be the input to postprocessing file)
    # make it the same name as original input file
    #output_filename = data_filename[:-4]  
    # Note: this assumes filename.xxx format for input and strips off the
    #       the dot and three character extension. 
#elif len(sys.argv) == 3:
#    # MCRT data to postprocess
#    data_filename = sys.argv[1]
#
#    # MCRT data to postprocess
#    output_filename = sys.argv[2]
else :
    print('*************************************')
    print('Error: Not enough command line arguments')
    print('$ python {} data_DIR'.format(sys.argv[0]))
    #print('$ python {} data_filename output_filename'.format(sys.argv[0]))
    print('*************************************')
    sys.exit(2)

# MCRT results text file (output of simulation)
volume_data_filename = data_dir+'/volume_solution.txt'

###############################################################################
#                                   Volume                                    #
###############################################################################

print("Reading {}".format(volume_data_filename) )

# Open file
f = open(volume_data_filename, 'rb')

# Read the number of cells in each direction
c = [int(a) for a in f.readline().split()]
Nx, Ny, Nz = c
print('Nx = {}'.format(Nx))
print('Ny = {}'.format(Ny))
print('Nz = {}'.format(Nz))

# Set up the numpy arrays for the solution
x   = np.zeros([Nx, Ny, Nz])
y   = np.zeros([Nx, Ny, Nz])
z   = np.zeros([Nx, Ny, Nz])
G   = np.zeros([Nx, Ny, Nz])

# Skip the header line
f.readline()

# Read the data
for line in f:
    # Turn the line into an array of floats
    c = [float(a) for a in line.split()]

    # Get the index of the cell
    i = int(c[0])
    j = int(c[1])
    k = int(c[2])

    # Fill in the data
    x[i,j,k] = c[3]
    y[i,j,k] = c[4]
    z[i,j,k] = c[5]
    G[i,j,k] = c[6]

# Close the file
f.close()

# Save the output to a compressed numpy file
output_filename = volume_data_filename[:-4]
np.savez(output_filename,
         x=x, y=y, z=z,
         G=G)

