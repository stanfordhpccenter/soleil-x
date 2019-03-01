# generated using paraview version 5.6.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

import argparse
import os
import shutil

parser = argparse.ArgumentParser()

# required
parser.add_argument('viz_ready_data_dir',
                    help='directory with the viz ready data')
# optional
parser.add_argument('--fluid_xmf_filename', nargs='?', const='out_fluid.xmf', default='out_fluid.xmf', type=str,
                    help='The fluid xmf file')
parser.add_argument('--particles_xmf_filename', nargs='?', const='out_particles.xmf', default='out_particles.xmf', type=str,
                    help='The partilces xmf file')
parser.add_argument('--image_size_x', nargs='?', const=1280, default=1280, type=int,
                    help='Sets the number of horizontal pixels in the produced images')
parser.add_argument('--image_size_y', nargs='?', const=720, default=720, type=int,
                    help='Sets the number of horizontal pixels in the produced images')
parser.add_argument('--out_dir', nargs='?', const=None, default=None, type=str,
                    help='Where to save the output of this script')
args = parser.parse_args()


###############################################################################
# Parse the input
###############################################################################
fluid_xmf_filename     = os.path.join(args.viz_ready_data_dir, args.fluid_xmf_filename)
particles_xmf_filename = os.path.join(args.viz_ready_data_dir, args.particles_xmf_filename)

if args.out_dir is not None:
  out_dir = args.out_dir
else:
  # set a default place to send the output if not specified on the command line
  out_dir = os.path.join(args.viz_ready_data_dir,'viz_output')

view_size_x = args.image_size_x
view_size_y = args.image_size_y

###############################################################################
# Create directory for the output if it is not already there
###############################################################################
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

###############################################################################
# Read in the data
###############################################################################
## Fluid data
out_fluidxmf = XDMFReader(FileNames=[fluid_xmf_filename])
out_fluidxmf.CellArrayStatus = ['pressure', 'rho', 'temperature', 'velocity']
out_fluidxmf.GridStatus = out_fluidxmf.GetPropertyValue('GridInfo').GetData()

# Partilces data
out_particlesxmf = XDMFReader(FileNames=[particles_xmf_filename])
out_particlesxmf.PointArrayStatus = ['diameter', 'temperature', 'velocity']
out_particlesxmf.GridStatus = out_particlesxmf.GetPropertyValue('GridInfo').GetData()

###############################################################################
# cosmetic stuff for line plots
###############################################################################

# To be used to set up how the plots look
series_line_style = ['pressure (stats)', '1',
                     'rho (stats)', '1',
                     'temperature (stats)', '1',
                     'velocity (0) (stats)', '1',
                     'velocity (1) (stats)', '1',
                     'velocity (2) (stats)', '1',
                     'velocity (Magnitude) (stats)', '1',
                     'vtkOriginalCellIds (stats)', '1',
                     'N (stats)', '1',
                     'Time (stats)', '1',
                     'vtkValidPointMask (stats)', '1']

series_marker_style = ['pressure (stats)', '4',
                       'rho (stats)', '4',
                       'temperature (stats)', '4',
                       'velocity (0) (stats)', '4',
                       'velocity (1) (stats)', '4',
                       'velocity (2) (stats)', '4',
                       'velocity (Magnitude) (stats)', '4',
                       'vtkOriginalCellIds (stats)', '4',
                       'N (stats)', '4',
                       'Time (stats)', '4',
                       'vtkValidPointMask (stats)', '4']

fluid_series_label = ['pressure (stats)',     'fluid pressure (stats)',
                      'rho (stats)',          'fluid rho (stats)',
                      'temperature (stats)',  'fluid temperature (stats)',
                      'velocity (0) (stats)', 'fluid velocity (0) (stats)',
                      'velocity (1) (stats)', 'fluid velocity (1) (stats)',
                      'velocity (2) (stats)', 'fluid velocity (2) (stats)',
                      'velocity (Magnitude) (stats)', 'fluid velocity (Magnitude) (stats)', 
                      'vtkOriginalCellIds (stats)', 'vtkOriginalCellIds (stats)', 'N (stats)', 'N (stats)', 'Time (stats)', 'Time (stats)', 'vtkValidPointMask (stats)', 'vtkValidPointMask (stats)']

particles_series_label = ['diameter (stats)',     'particle diameter (stats)',
                          'temperature (stats)',  'partilce temperature (stats)',
                          'velocity (0) (stats)', 'partilce velocity (0) (stats)',
                          'velocity (1) (stats)', 'partilce velocity (1) (stats)',
                          'velocity (2) (stats)', 'partilce velocity (2) (stats)',
                          'velocity (Magnitude) (stats)', 'particle velocity (Magnitude) (stats)',
                          'vtkOriginalPointIds (stats)', 'vtkOriginalPointIds (stats)',
                          'X (stats)', 'X (stats)',
                          'Y (stats)', 'Y (stats)',
                          'Z (stats)', 'Z (stats)',
                          'N (stats)', 'N (stats)',
                          'Time (stats)', 'Time (stats)',
                          'vtkValidPointMask (stats)', 'vtkValidPointMask (stats)']

###############################################################################
# Plot data max and min vs time
###############################################################################

# To be used in the loop
fluid_variable_names    = ['temperature', 'velocity (Magnitude)', 'velocity (0)','velocity (1)','velocity (2)', 'rho', 'pressure']
particle_variable_names = ['temperature', 'velocity (Magnitude)', 'velocity (0)','velocity (1)','velocity (2)',]

query_operations = ['max', 'min']

# key value pairs where (key = variable name) and (value = stuff to be used in loops below)
query_variable_names = {'temperature'          : 'temperature',
                        'velocity (Magnitude)' : 'mag(velocity)',
                        'velocity (0)'         : 'velocity[:,0]',
                        'velocity (1)'         : 'velocity[:,1]',
                        'velocity (2)'         : 'velocity[:,2]',
                        'rho'                  : 'rho',
                        'pressure'             : 'pressure'}

label_names = {'temperature'          : 'Temperature',
               'velocity (Magnitude)' : 'Velocity_Magnitude',
               'velocity (0)'         : 'u',
               'velocity (1)'         : 'v',
               'velocity (2)'         : 'w',
               'rho'                  : 'rho',
               'pressure'             : 'Pressure'}
     

for variable_name in fluid_variable_names:
  for query_operation in query_operations:

    query_string = '{} == {}({})'.format(query_variable_names[variable_name], query_operation, query_variable_names[variable_name])

    # Get the fluid selection
    fluid_selection = SelectionQuerySource()
    fluid_selection.QueryString = query_string 

    # Create a new 'Quartile Chart View'
    quartileChartView1 = CreateView('QuartileChartView')
    quartileChartView1.ViewSize = [view_size_x,view_size_y]
    
    # Plot the fluid selection over time
    plotSelectionOverTime1 = PlotSelectionOverTime(Input=out_fluidxmf,
                                                   Selection=fluid_selection)
    # show data in view and set up how it looks
    plotSelectionOverTime1Display = Show(plotSelectionOverTime1, quartileChartView1)
    plotSelectionOverTime1Display.AttributeType = 'Row Data'
    plotSelectionOverTime1Display.UseIndexForXAxis = 0
    plotSelectionOverTime1Display.XArrayName = 'Time'
    plotSelectionOverTime1Display.SeriesVisibility = [variable_name + ' (stats)']
    plotSelectionOverTime1Display.SeriesLabel = fluid_series_label 
    plotSelectionOverTime1Display.SeriesLineStyle = series_line_style 
    plotSelectionOverTime1Display.SeriesMarkerStyle = series_marker_style 

    # Get the partilce selection if needed
    if variable_name in particle_variable_names: 
      particle_selection = SelectionQuerySource()
      particle_selection.QueryString = query_string 
      particle_selection.FieldType = 'POINT' 

      plotSelectionOverTime2 = PlotSelectionOverTime(Input=out_particlesxmf,
                                                     Selection=particle_selection)
      # show data in view and set up how it looks
      plotSelectionOverTime2Display = Show(plotSelectionOverTime2, quartileChartView1)
      plotSelectionOverTime2Display.AttributeType = 'Row Data'
      plotSelectionOverTime2Display.UseIndexForXAxis = 0
      plotSelectionOverTime2Display.XArrayName = 'Time'
      plotSelectionOverTime2Display.SeriesVisibility = [variable_name + ' (stats)']
      plotSelectionOverTime2Display.SeriesLabel = particles_series_label 
      plotSelectionOverTime2Display.SeriesLineStyle = series_line_style
      plotSelectionOverTime2Display.SeriesMarkerStyle = series_marker_style

    # Properties modified on quartileChartView1
    quartileChartView1.ChartTitle = '{} {} vs Time'.format(query_operation, label_names[variable_name])
    
    # Properties modified on quartileChartView1
    quartileChartView1.LeftAxisTitle = '{} {}'.format(query_operation, label_names[variable_name])
    
    # Properties modified on quartileChartView1
    quartileChartView1.BottomAxisTitle = 'Restart File Number'
    
    # Properties modified on quartileChartView1
    quartileChartView1.ShowLegend = 1
    
    # Properties modified on quartileChartView1
    quartileChartView1.LegendLocation = 'TopLeft'
    
    # update the view to ensure updated data information
    quartileChartView1.Update()

    # save screenshot
    screenshot_name = os.path.join(out_dir,'{}_{}_vs_time.png'.format(label_names[variable_name],query_operation))
    SaveScreenshot(screenshot_name,
                   quartileChartView1,
                   ImageResolution=[view_size_x, view_size_y],
                   CompressionLevel='0')
    print('Saved file: {}'.format(screenshot_name))
   

