# generated using paraview version 5.6.0

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
parser.add_argument('--fluid_xmf_filename', nargs='?', const='fluid.xmf', default='fluid.xmf', type=str,
                    help='The fluid xmf file')
parser.add_argument('--particles_xmf_filename', nargs='?', const='particles.xmf', default='particles.xmf', type=str,
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
fluid = XDMFReader(FileNames=[fluid_xmf_filename])

# Partilces data
particles = XDMFReader(FileNames=[particles_xmf_filename])

###############################################################################
# Plot the data vs. time
###############################################################################

def plot_stats_vs_time(input_source, selection, variable_name, paraview_variable_name):


  # Plot the selection applied to the input source over time
  plotSelectionOverTime1 = PlotSelectionOverTime(Input=input_source,
                                                 Selection=selection)
  # Create a new 'Quartile Chart View'
  quartileChartView1 = CreateView('QuartileChartView')

  # show data in view and set up how it looks
  plotSelectionOverTime1Display = Show(plotSelectionOverTime1, quartileChartView1)
  plotSelectionOverTime1Display.AttributeType = 'Row Data'
  plotSelectionOverTime1Display.UseIndexForXAxis = 0
  plotSelectionOverTime1Display.XArrayName = 'Time'
  plotSelectionOverTime1Display.SeriesVisibility = [paraview_variable_name + ' (stats)']
  plotSelectionOverTime1Display.SeriesMarkerStyle  = [paraview_variable_name + ' (stats)', '4']

  # Properties modified on quartileChartView1
  quartileChartView1.ChartTitle = '{} vs Time in Entire Domain'.format(variable_name)

  # Properties modified on quartileChartView1
  quartileChartView1.LeftAxisTitle = '{}'.format(variable_name)

  # Properties modified on quartileChartView1
  quartileChartView1.BottomAxisTitle = 'Restart File Number'

  # Properties modified on quartileChartView1
  quartileChartView1.ShowLegend = 1

  # Properties modified on quartileChartView1
  quartileChartView1.LegendLocation = 'TopLeft'

  # update the view to ensure updated data information
  quartileChartView1.Update()

  return quartileChartView1

# List of which plots to make
fluid_variable_names = ['Temperature',
                        'Velocity_Magnitude',
                        'u',
                        'v',
                        'w',
                        'rho',
                        'Pressure']

particle_variable_names = ['Temperature',
                           'Velocity_Magnitude',
                           'u',
                           'v',
                           'w',
                           'Number_Particles']

# translates to what the variables are called in paraview
paraview_fluid_variable_names = {'Temperature'        : 'temperature',
                                 'Velocity_Magnitude' : 'velocity (Magnitude)',
                                 'u'                  : 'velocity (0)',
                                 'v'                  : 'velocity (1)',
                                 'w'                  : 'velocity (2)',
                                 'rho'                : 'rho',
                                 'Pressure'           : 'pressure'}

paraview_particle_variable_names = {'Temperature'        : 'temperature',
                                    'Velocity_Magnitude' : 'velocity (Magnitude)',
                                    'u'                  : 'velocity (0)',
                                    'v'                  : 'velocity (1)',
                                    'w'                  : 'velocity (2)',
                                    'Number_Particles'   : 'N'}


# Get the selections set up
# hack query to select all everything in the source
query_string = 'id >= -1'
# select the points for the particles
points_selection = SelectionQuerySource()
points_selection.QueryString = query_string
points_selection.FieldType = 'POINT'
# select the cells for the cells
cells_selection = SelectionQuerySource()
cells_selection.QueryString = query_string
cells_selection.FieldType = 'CELL'


for variable_name in fluid_variable_names:

   chartView = plot_stats_vs_time(fluid,
                                  cells_selection,
                                  variable_name,
                                  paraview_fluid_variable_names[variable_name])

   # make the view the same size as the screen shot that you want to take
   chartView.ViewSize = [view_size_x,view_size_y]
   chartView.Update()

   # save screenshot
   screenshot_filename = os.path.join(out_dir,'fluid_{}_vs_time.png'.format(variable_name))
   SaveScreenshot(screenshot_filename,
                  chartView,
                  ImageResolution=[view_size_x, view_size_y],
                  CompressionLevel='0')
   print('Saved file: {}'.format(screenshot_filename))


for variable_name in particle_variable_names:

   chartView = plot_stats_vs_time(particles,
                                  points_selection,
                                  variable_name,
                                  paraview_particle_variable_names[variable_name])

   # make the view the same size as the screen shot that you want to take
   chartView.ViewSize = [view_size_x,view_size_y]
   chartView.Update()

   # save screenshot
   screenshot_filename = os.path.join(out_dir,'particle_{}_vs_time.png'.format(variable_name))
   SaveScreenshot(screenshot_filename,
                  chartView,
                  ImageResolution=[view_size_x, view_size_y],
                  CompressionLevel='0')
   print('Saved file: {}'.format(screenshot_filename))
