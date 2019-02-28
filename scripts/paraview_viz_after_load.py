# trace generated using paraview version 5.6.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

###############################################################################
# Variables that are assumed
###############################################################################
# Get from extents in the future:

fluid_xmf_filename     = '/home/lalo/Desktop/debug-crash-data/lf32/job0/sample1/viz_ready_output/out_fluid.xmf'
particles_xmf_filename = '/home/lalo/Desktop/debug-crash-data/lf32/job0/sample1/viz_ready_output/out_particles.xmf'

channel_center_point = [0.08, 0.02, 0.02]
z_hat = [0.0, 0.0, 1.0]

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

## Find and name the sources if the files were read before this scrit was run
#out_fluidxmf     = FindSource('out_fluid.xmf')
#out_particlesxmf = FindSource('out_particles.xmf')


RenameSource('fluid', out_fluidxmf)
RenameSource('particles', out_particlesxmf)

###############################################################################
# Set up render view
###############################################################################

# Get or create render view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1612, 673]


###############################################################################
# Set up layout
###############################################################################
layout1 = GetLayout()
RenameLayout('Main Layout', layout1)


###############################################################################
# Set up animation 
###############################################################################

# get animation scene
animationScene1 = GetAnimationScene()
# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()
# Properties modified on animationScene1
animationScene1.PlayMode = 'Real Time'
# Properties modified on animationScene1
animationScene1.Duration = 5
# Make the animation loop
animationScene1.Loop = 1


###############################################################################
# Set up view of raw fluid data  
###############################################################################

## show data in view
out_fluidxmfDisplay = Show(out_fluidxmf, renderView1)
out_fluidxmfDisplay.Representation = 'Outline'

# reset view to fit data
renderView1.ResetCamera()

## update the view to ensure updated data information
renderView1.Update()


###############################################################################
# Set up view of raw particle data
###############################################################################

# Use the last dump. Because the first one may not have particles in it
# and also because the last time step will likely be the hottest. So good for
# using to scale the color bar.
animationScene1.GoToLast()

# show data in view
out_particlesxmfDisplay = Show(out_particlesxmf, renderView1)
out_particlesxmfDisplay.Representation = 'Points'
out_particlesxmfDisplay.ColorArrayName = ['temperature']
# rescale color and/or opacity maps used to include current data range
out_particlesxmfDisplay.RescaleTransferFunctionToDataRange(True, False)
# show color bar/color legend
out_particlesxmfDisplay.SetScalarBarVisibility(renderView1, True)

# hide the particles from the view
Hide(out_particlesxmf, renderView1)


################################################################################
## Create fluid z plane slice though mid-point of simulation
################################################################################
# create a new 'Slice'
z_mid_slice = Slice(Input=out_fluidxmf)
z_mid_slice.SliceType = 'Plane'
z_mid_slice.SliceOffsetValues = [0.0]
z_mid_slice.SliceType.Origin = channel_center_point
z_mid_slice.SliceType.Normal = z_hat

# show data in view
z_mid_slice_Display = Show(z_mid_slice, renderView1)

# get color transfer function/color map for 'temperature'
temperatureLUT = GetColorTransferFunction('temperature')
temperatureLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941,
                            0.5, 0.865003, 0.865003, 0.865003, 
                            1.0, 0.705882, 0.0156863, 0.14902]
temperatureLUT.ScalarRangeInitialized = 1.0
# rescale color and/or opacity maps used to include current data range (of the last time step)
# would be better to rescale over the whole time range but for now you will just have to do that via the GUI
animationScene1.GoToLast()
z_mid_slice_Display.RescaleTransferFunctionToDataRange(True, False)

# display properties.
z_mid_slice_Display.Representation = 'Surface'
z_mid_slice_Display.ColorArrayName = ['CELLS', 'temperature']
z_mid_slice_Display.LookupTable = temperatureLUT

# show color bar/color legend
z_mid_slice_Display.SetScalarBarVisibility(renderView1, True)

# rename source object
RenameSource('z_mid_slice', z_mid_slice)

###############################################################################
# Create fluid y plane slice though mid-point of simulation
###############################################################################
# create a new 'Slice'
y_mid_slice = Slice(Input=out_fluidxmf)
y_mid_slice.SliceType = 'Plane'
y_mid_slice.SliceOffsetValues = [0.0]
y_mid_slice.SliceType.Origin = channel_center_point
y_mid_slice.SliceType.Normal = [0,1,0]

# show data in view
y_mid_slice_Display = Show(y_mid_slice, renderView1)

# get color transfer function/color map for 'temperature'
temperatureLUT = GetColorTransferFunction('temperature')
temperatureLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941,
                            0.5, 0.865003, 0.865003, 0.865003, 
                            1.0, 0.705882, 0.0156863, 0.14902]
temperatureLUT.ScalarRangeInitialized = 1.0
# rescale color and/or opacity maps used to include current data range (of the last time step)
# would be better to rescale over the whole time range but for now you will just have to do that via the GUI
animationScene1.GoToLast()
y_mid_slice_Display.RescaleTransferFunctionToDataRange(True, False)

# display properties.
y_mid_slice_Display.Representation = 'Surface'
y_mid_slice_Display.ColorArrayName = ['CELLS', 'temperature']
y_mid_slice_Display.LookupTable = temperatureLUT

# show color bar/color legend
y_mid_slice_Display.SetScalarBarVisibility(renderView1, True)

## rename source object
RenameSource('y_mid_slice', y_mid_slice)

Hide(y_mid_slice, renderView1)

###############################################################################
# Inlet Slice
###############################################################################

# create a new 'Slice'
inlet_slice = Slice(Input=out_fluidxmf)
inlet_slice.SliceType = 'Plane'
inlet_slice.SliceOffsetValues = [0.0]

# Properties modified on inlet_slice.SliceType
inlet_slice.SliceType.Origin = [0.0, 0.02, 0.02]

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
renderView1.Update()

# show data in view
inlet_slice_Display = Show(inlet_slice, renderView1)
inlet_slice_Display.Representation = 'Surface'
inlet_slice_Display.ColorArrayName = ['CELLS', 'temperature']
inlet_slice_Display.LookupTable = temperatureLUT

# set scalar coloring
ColorBy(inlet_slice_Display, ('CELLS', 'temperature'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(temperatureLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
inlet_slice_Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
inlet_slice_Display.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(inlet_slice, renderView1)

# rename source object
RenameSource('inlet_slice', inlet_slice)

###############################################################################
# Outlet Slice
###############################################################################

# create a new 'Slice'
outlet_slice = Slice(Input=out_fluidxmf)
outlet_slice.SliceType = 'Plane'
outlet_slice.SliceOffsetValues = [0.0]

# Properties modified on outlet_slice.SliceType
outlet_slice.SliceType.Origin = [0.16, 0.02, 0.02]

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
renderView1.Update()

# show data in view
outlet_slice_Display = Show(outlet_slice, renderView1)
outlet_slice_Display.Representation = 'Surface'
outlet_slice_Display.ColorArrayName = ['CELLS', 'temperature']
outlet_slice_Display.LookupTable = temperatureLUT

# set scalar coloring
ColorBy(outlet_slice_Display, ('CELLS', 'temperature'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(temperatureLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
outlet_slice_Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
outlet_slice_Display.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(outlet_slice, renderView1)

# rename source object
RenameSource('outlet_slice', outlet_slice)



###############################################################################
# Create the particle clip
###############################################################################

# create a new 'Clip'
clip1 = Clip(Input=out_particlesxmf)
clip1.ClipType = 'Plane'
clip1.Scalars = ['POINTS', 'temperature']

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = channel_center_point

# Properties modified on clip1.ClipType
clip1.ClipType.Normal = z_hat

# show data in view
clip1Display = Show(clip1, renderView1)

## trace defaults for the display properties.
clip1Display.Representation = 'Surface'
clip1Display.ColorArrayName = ['POINTS', 'temperature']
clip1Display.LookupTable = temperatureLUT

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()


###############################################################################
# Fluid temperature threshold
###############################################################################

# create a new 'Threshold'
threshold1 = Threshold(Input=out_fluidxmf)

# Get the bounds of the temperature on the last time step
animationScene1.GoToLast()
temperature_min = out_fluidxmf.CellData['temperature'].GetRange()[0]
temperature_max = out_fluidxmf.CellData['temperature'].GetRange()[1]

# Properties modified on threshold1
threshold1.Scalars = ['CELLS', 'temperature']
threshold1.ThresholdRange = [0.8*temperature_max, temperature_max]

# show data in view
threshold1Display = Show(threshold1, renderView1)

## trace defaults for the display properties.
threshold1Display.Representation = 'Surface'
threshold1Display.ColorArrayName = ['CELLS', 'temperature']
threshold1Display.LookupTable = temperatureLUT

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# hide data in view
Hide(threshold1, renderView1)

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(temperatureLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
threshold1Display.RescaleTransferFunctionToDataRange(True, False)



###############################################################################
# Stats in spread sheet on layout 2
###############################################################################

layout2 = CreateLayout('Layout #2')
RenameLayout('Discriptive Statistics Layout', layout2)

###################
# Fluid Stats
###################
# create a new 'Descriptive Statistics'
descriptiveStatistics1 = DescriptiveStatistics(Input=out_fluidxmf,
                                               ModelInput=None)
descriptiveStatistics1.AttributeMode = 'Cell Data'
descriptiveStatistics1.VariablesofInterest = ['pressure', 'rho', 'temperature', 'velocity']

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.ColumnToSort = ''
spreadSheetView1.BlockSize = 1024L
# uncomment following to set a specific view size
# spreadSheetView1.ViewSize = [400, 400]

# place view in the layout
layout2.AssignView(0, spreadSheetView1)

# show data in view
descriptiveStatistics1Display = Show(descriptiveStatistics1, spreadSheetView1)

# trace defaults for the display properties.
descriptiveStatistics1Display.CompositeDataSetIndex = [1]

# update the view to ensure updated data information
spreadSheetView1.Update()

###################
# Particle Stats
###################

# split cell
layout2.SplitVertical(0, 0.5)

# Create a new 'SpreadSheet View'
spreadSheetView2 = CreateView('SpreadSheetView')
spreadSheetView2.ColumnToSort = ''
spreadSheetView2.BlockSize = 1024L
# uncomment following to set a specific view size
# spreadSheetView2.ViewSize = [400, 400]

# place view in the layout
layout2.AssignView(2, spreadSheetView2)

# create a new 'Descriptive Statistics'
descriptiveStatistics2 = DescriptiveStatistics(Input=out_particlesxmf,
                                               ModelInput=None)
descriptiveStatistics2.VariablesofInterest = ['diameter', 'temperature', 'velocity']

# show data in view
descriptiveStatistics2Display = Show(descriptiveStatistics2, spreadSheetView2)

# trace defaults for the display properties.
descriptiveStatistics2Display.CompositeDataSetIndex = [1]

# update the view to ensure updated data information
spreadSheetView2.Update()

###############################################################################
# Fluid Temporal Statistics
###############################################################################

# create a new 'Temporal Statistics'
temporalStatistics1 = TemporalStatistics(Input=out_fluidxmf)

# show data in view
temporalStatistics1Display = Show(temporalStatistics1, renderView1)

# trace defaults for the display properties.
temporalStatistics1Display.Representation = 'Surface'
temporalStatistics1Display.ColorArrayName = [None, '']
# set scalar coloring
ColorBy(temporalStatistics1Display, ('CELLS', 'temperature_maximum'))

# rescale color and/or opacity maps used to include current data range
temporalStatistics1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
temporalStatistics1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'temperature_maximum'
temperature_maximumLUT = GetColorTransferFunction('temperature_maximum')

# get opacity transfer function/opacity map for 'temperature_maximum'
temperature_maximumPWF = GetOpacityTransferFunction('temperature_maximum')

# hide data in view
Hide(temporalStatistics1, renderView1)

# update the view to ensure updated data information
renderView1.Update()


###############################################################################
# Fluid max temperature vs. time
###############################################################################
# create a new 'Plot Selection Over Time'

SetActiveSource(out_fluidxmf)
SetActiveView(renderView1)

max_fluid_temperature_selection = SelectionQuerySource()
max_fluid_temperature_selection.QueryString = 'temperature  == max(temperature)'

plotSelectionOverTime1 = PlotSelectionOverTime(Input=out_fluidxmf,
                                               Selection=max_fluid_temperature_selection)

# Create a new 'Quartile Chart View'
quartileChartView1 = CreateView('QuartileChartView')
#quartileChartView1.ViewSize = [733, 224]
#quartileChartView1.ChartTitleFontFile = ''
#quartileChartView1.LeftAxisTitleFontFile = ''
#quartileChartView1.LeftAxisRangeMaximum = 6.66
#quartileChartView1.LeftAxisLabelFontFile = ''
#quartileChartView1.BottomAxisTitleFontFile = ''
#quartileChartView1.BottomAxisRangeMaximum = 6.66
#quartileChartView1.BottomAxisLabelFontFile = ''
#quartileChartView1.RightAxisRangeMaximum = 6.66
#quartileChartView1.RightAxisLabelFontFile = ''
#quartileChartView1.TopAxisTitleFontFile = ''
#quartileChartView1.TopAxisRangeMaximum = 6.66
#quartileChartView1.TopAxisLabelFontFile = ''

# place view in the layout
# split cell
#layout1.SplitHorizontal(0, 0.5)
layout1.AssignView(2, quartileChartView1)

# show data in view
plotSelectionOverTime1Display = Show(plotSelectionOverTime1, quartileChartView1)

# trace defaults for the display properties.
plotSelectionOverTime1Display.AttributeType = 'Row Data'
plotSelectionOverTime1Display.UseIndexForXAxis = 0
plotSelectionOverTime1Display.XArrayName = 'Time'
plotSelectionOverTime1Display.SeriesVisibility = ['pressure (stats)', 'rho (stats)', 'temperature (stats)', 'velocity (Magnitude) (stats)']
plotSelectionOverTime1Display.SeriesLabel = ['pressure (stats)',     'fluid pressure (stats)',
                                             'rho (stats)',          'fluid rho (stats)',
                                             'temperature (stats)',  'fluid temperature (stats)',
                                             'velocity (0) (stats)', 'fluid velocity (0) (stats)',
                                             'velocity (1) (stats)', 'fluid velocity (1) (stats)',
                                             'velocity (2) (stats)', 'fluid velocity (2) (stats)',
                                             'velocity (Magnitude) (stats)', 'fluid velocity (Magnitude) (stats)', 'vtkOriginalCellIds (stats)', 'vtkOriginalCellIds (stats)', 'N (stats)', 'N (stats)', 'Time (stats)', 'Time (stats)', 'vtkValidPointMask (stats)', 'vtkValidPointMask (stats)']
plotSelectionOverTime1Display.SeriesColor = ['pressure (stats)', '0', '0', '0', 'rho (stats)', '0.89', '0.1', '0.11', 'temperature (stats)', '0.22', '0.49', '0.72', 'velocity (0) (stats)', '0.3', '0.69', '0.29', 'velocity (1) (stats)', '0.6', '0.31', '0.64', 'velocity (2) (stats)', '1', '0.5', '0', 'velocity (Magnitude) (stats)', '0.65', '0.34', '0.16', 'vtkOriginalCellIds (stats)', '0', '0', '0', 'N (stats)', '0.89', '0.1', '0.11', 'Time (stats)', '0.22', '0.49', '0.72', 'vtkValidPointMask (stats)', '0.3', '0.69', '0.29']
plotSelectionOverTime1Display.SeriesPlotCorner = ['pressure (stats)', '0', 'rho (stats)', '0', 'temperature (stats)', '0', 'velocity (0) (stats)', '0', 'velocity (1) (stats)', '0', 'velocity (2) (stats)', '0', 'velocity (Magnitude) (stats)', '0', 'vtkOriginalCellIds (stats)', '0', 'N (stats)', '0', 'Time (stats)', '0', 'vtkValidPointMask (stats)', '0']
plotSelectionOverTime1Display.SeriesLabelPrefix = ''
plotSelectionOverTime1Display.SeriesLineStyle = ['pressure (stats)', '1', 'rho (stats)', '1', 'temperature (stats)', '1', 'velocity (0) (stats)', '1', 'velocity (1) (stats)', '1', 'velocity (2) (stats)', '1', 'velocity (Magnitude) (stats)', '1', 'vtkOriginalCellIds (stats)', '1', 'N (stats)', '1', 'Time (stats)', '1', 'vtkValidPointMask (stats)', '1']
plotSelectionOverTime1Display.SeriesLineThickness = ['pressure (stats)', '2', 'rho (stats)', '2', 'temperature (stats)', '2', 'velocity (0) (stats)', '2', 'velocity (1) (stats)', '2', 'velocity (2) (stats)', '2', 'velocity (Magnitude) (stats)', '2', 'vtkOriginalCellIds (stats)', '2', 'N (stats)', '2', 'Time (stats)', '2', 'vtkValidPointMask (stats)', '2']
plotSelectionOverTime1Display.SeriesMarkerStyle = ['pressure (stats)', '0', 'rho (stats)', '0', 'temperature (stats)', '0', 'velocity (0) (stats)', '0', 'velocity (1) (stats)', '0', 'velocity (2) (stats)', '0', 'velocity (Magnitude) (stats)', '0', 'vtkOriginalCellIds (stats)', '0', 'N (stats)', '0', 'Time (stats)', '0', 'vtkValidPointMask (stats)', '0']


## get color transfer function/color map for 'temperature'
#temperatureLUT = GetColorTransferFunction('temperature')
#
## Rescale transfer function
#temperatureLUT.RescaleTransferFunction(288.949059081, 1676.09476296)
#
## get opacity transfer function/opacity map for 'temperature'
#temperaturePWF = GetOpacityTransferFunction('temperature')
#
## Rescale transfer function
#temperaturePWF.RescaleTransferFunction(288.949059081, 1676.09476296)

# Properties modified on plotSelectionOverTime1Display
plotSelectionOverTime1Display.SeriesVisibility = ['rho (stats)', 'temperature (stats)', 'velocity (Magnitude) (stats)']
plotSelectionOverTime1Display.SeriesColor = ['pressure (stats)', '0', '0', '0', 'rho (stats)', '0.889998', '0.100008', '0.110002', 'temperature (stats)', '0.220005', '0.489998', '0.719997', 'velocity (0) (stats)', '0.300008', '0.689998', '0.289998', 'velocity (1) (stats)', '0.6', '0.310002', '0.639994', 'velocity (2) (stats)', '1', '0.500008', '0', 'velocity (Magnitude) (stats)', '0.650004', '0.340002', '0.160006', 'vtkOriginalCellIds (stats)', '0', '0', '0', 'N (stats)', '0.889998', '0.100008', '0.110002', 'Time (stats)', '0.220005', '0.489998', '0.719997', 'vtkValidPointMask (stats)', '0.300008', '0.689998', '0.289998']
plotSelectionOverTime1Display.SeriesPlotCorner = ['N (stats)', '0', 'Time (stats)', '0', 'pressure (stats)', '0', 'rho (stats)', '0', 'temperature (stats)', '0', 'velocity (0) (stats)', '0', 'velocity (1) (stats)', '0', 'velocity (2) (stats)', '0', 'velocity (Magnitude) (stats)', '0', 'vtkOriginalCellIds (stats)', '0', 'vtkValidPointMask (stats)', '0']
plotSelectionOverTime1Display.SeriesLineStyle = ['N (stats)', '1', 'Time (stats)', '1', 'pressure (stats)', '1', 'rho (stats)', '1', 'temperature (stats)', '1', 'velocity (0) (stats)', '1', 'velocity (1) (stats)', '1', 'velocity (2) (stats)', '1', 'velocity (Magnitude) (stats)', '1', 'vtkOriginalCellIds (stats)', '1', 'vtkValidPointMask (stats)', '1']
plotSelectionOverTime1Display.SeriesLineThickness = ['N (stats)', '2', 'Time (stats)', '2', 'pressure (stats)', '2', 'rho (stats)', '2', 'temperature (stats)', '2', 'velocity (0) (stats)', '2', 'velocity (1) (stats)', '2', 'velocity (2) (stats)', '2', 'velocity (Magnitude) (stats)', '2', 'vtkOriginalCellIds (stats)', '2', 'vtkValidPointMask (stats)', '2']
plotSelectionOverTime1Display.SeriesMarkerStyle = ['N (stats)', '0', 'Time (stats)', '0', 'pressure (stats)', '0', 'rho (stats)', '0', 'temperature (stats)', '4', 'velocity (0) (stats)', '0', 'velocity (1) (stats)', '0', 'velocity (2) (stats)', '0', 'velocity (Magnitude) (stats)', '0', 'vtkOriginalCellIds (stats)', '0', 'vtkValidPointMask (stats)', '0']

# Properties modified on plotSelectionOverTime1Display
plotSelectionOverTime1Display.SeriesVisibility = ['temperature (stats)']

# Properties modified on quartileChartView1
quartileChartView1.ChartTitle = 'Max Temperature vs Time'

# Properties modified on quartileChartView1
#quartileChartView1.LeftAxisTitle = 'Temperature'

# Properties modified on quartileChartView1
quartileChartView1.BottomAxisTitle = 'Restart File Number'

# Properties modified on quartileChartView1
quartileChartView1.ShowLegend = 1

# Properties modified on quartileChartView1
quartileChartView1.LegendLocation = 'TopLeft'

# update the view to ensure updated data information
quartileChartView1.Update()

###############################################################################
# particle max temperature vs. time
###############################################################################
# You may get some complaints if there is no particle data in the first dump.

max_particle_temperature_selection = SelectionQuerySource()
max_particle_temperature_selection.QueryString = 'temperature  == max(temperature)'
max_particle_temperature_selection.FieldType = 'POINT' 

plotSelectionOverTime2 = PlotSelectionOverTime(Input=out_particlesxmf,
                                               Selection=max_particle_temperature_selection)
# show data in view
plotSelectionOverTime2Display = Show(plotSelectionOverTime2, quartileChartView1)


# trace defaults for the display properties.
plotSelectionOverTime2Display.AttributeType = 'Row Data'
plotSelectionOverTime2Display.UseIndexForXAxis = 0
plotSelectionOverTime2Display.XArrayName = 'Time'
plotSelectionOverTime2Display.SeriesVisibility = ['diameter (stats)', 'temperature (stats)', 'velocity (Magnitude) (stats)']
# Properties modified2on plotSelectionOverTime1Display
plotSelectionOverTime2Display.SeriesLabel = ['diameter (stats)',     'diameter (stats)',
                                             'temperature (stats)',  'partilce temperature (stats)',
                                             'velocity (0) (stats)', 'partilce velocity (0) (stats)',
                                             'velocity (1) (stats)', 'partilce velocity (1) (stats)',
                                             'velocity (2) (stats)', 'partilce velocity (2) (stats)',
                                             'velocity (Magnitude) (stats)', 'velocity (Magnitude) (stats)',
                                             'vtkOriginalPointIds (stats)', 'vtkOriginalPointIds (stats)', 'X (stats)', 'X (stats)', 'Y (stats)', 'Y (stats)', 'Z (stats)', 'Z (stats)', 'N (stats)', 'N (stats)', 'Time (stats)', 'Time (stats)', 'vtkValidPointMask (stats)', 'vtkValidPointMask (stats)']
plotSelectionOverTime2Display.SeriesColor = ['diameter (stats)', '0', '0', '0', 'temperature (stats)', '0.89', '0.1', '0.11', 'velocity (0) (stats)', '0.22', '0.49', '0.72', 'velocity (1) (stats)', '0.3', '0.69', '0.29', 'velocity (2) (stats)', '0.6', '0.31', '0.64', 'velocity (Magnitude) (stats)', '1', '0.5', '0', 'vtkOriginalPointIds (stats)', '0.65', '0.34', '0.16', 'X (stats)', '0', '0', '0', 'Y (stats)', '0.89', '0.1', '0.11', 'Z (stats)', '0.22', '0.49', '0.72', 'N (stats)', '0.3', '0.69', '0.29', 'Time (stats)', '0.6', '0.31', '0.64', 'vtkValidPointMask (stats)', '1', '0.5', '0']
plotSelectionOverTime2Display.SeriesPlotCorner = ['diameter (stats)', '0', 'temperature (stats)', '0', 'velocity (0) (stats)', '0', 'velocity (1) (stats)', '0', 'velocity (2) (stats)', '0', 'velocity (Magnitude) (stats)', '0', 'vtkOriginalPointIds (stats)', '0', 'X (stats)', '0', 'Y (stats)', '0', 'Z (stats)', '0', 'N (stats)', '0', 'Time (stats)', '0', 'vtkValidPointMask (stats)', '0']
plotSelectionOverTime2Display.SeriesLabelPrefix = ''
plotSelectionOverTime2Display.SeriesLineStyle = ['diameter (stats)',     '1',
                                                 'temperature (stats)',  '1',
                                                 'velocity (0) (stats)', '1', 
                                                 'velocity (1) (stats)', '1',
                                                 'velocity (2) (stats)', '1',
                                                 'velocity (Magnitude) (stats)', '1', 'vtkOriginalPointIds (stats)', '1', 'X (stats)', '1', 'Y (stats)', '1', 'Z (stats)', '1', 'N (stats)', '1', 'Time (stats)', '1', 'vtkValidPointMask (stats)', '1']
plotSelectionOverTime2Display.SeriesLineThickness = ['diameter (stats)', '2', 'temperature (stats)', '2', 'velocity (0) (stats)', '2', 'velocity (1) (stats)', '2', 'velocity (2) (stats)', '2', 'velocity (Magnitude) (stats)', '2', 'vtkOriginalPointIds (stats)', '2', 'X (stats)', '2', 'Y (stats)', '2', 'Z (stats)', '2', 'N (stats)', '2', 'Time (stats)', '2', 'vtkValidPointMask (stats)', '2']
plotSelectionOverTime2Display.SeriesMarkerStyle = ['diameter (stats)', '0', 'temperature (stats)', '4', 'velocity (0) (stats)', '0', 'velocity (1) (stats)', '0', 'velocity (2) (stats)', '0', 'velocity (Magnitude) (stats)', '0', 'vtkOriginalPointIds (stats)', '0', 'X (stats)', '0', 'Y (stats)', '0', 'Z (stats)', '0', 'N (stats)', '0', 'Time (stats)', '0', 'vtkValidPointMask (stats)', '0']

plotSelectionOverTime2Display.SeriesVisibility = ['temperature (stats)']

## Rescale transfer function
#temperatureLUT.RescaleTransferFunction(288.949059081, 1676.09476296)
#
## Rescale transfer function
#temperaturePWF.RescaleTransferFunction(288.949059081, 1676.09476296)

# update the view to ensure updated data information
quartileChartView1.Update()

# save screenshot
SaveScreenshot('/home/lalo/Desktop/max_temperature_vs_time.png',
               quartileChartView1,
               ImageResolution=[731, 646],
               CompressionLevel='0')


################################################################################
## Plot over line
################################################################################
## set active source
#SetActiveSource(out_fluidxmf)
#
## create a new 'Plot Over Line'
#plotOverLine1 = PlotOverLine(Input=out_fluidxmf,
#                             Source='High Resolution Line Source')
#
## init the 'High Resolution Line Source' selected for 'Source'
#plotOverLine1.Source.Point1 = [-0.00125, 0.02, 0.02]
#plotOverLine1.Source.Point2 = [ 0.16125, 0.02, 0.02]
#
## show data in view
#plotOverLine1Display = Show(plotOverLine1, renderView1)
#
## trace defaults for the display properties.
#plotOverLine1Display.Representation = 'Surface'
#
## Create a new 'Line Chart View'
#lineChartView1 = CreateView('XYChartView')
#lineChartView1.ViewSize = [758, 321]
#
## place view in the layout
##layout1.AssignView(4, lineChartView1)
#
## show data in view
#plotOverLine1Display_1 = Show(plotOverLine1, lineChartView1)
#
## trace defaults for the display properties.
#plotOverLine1Display_1.CompositeDataSetIndex = [0]
#plotOverLine1Display_1.UseIndexForXAxis = 0
#plotOverLine1Display_1.XArrayName = 'arc_length'
##plotOverLine1Display_1.SeriesVisibility = ['pressure', 'rho', 'temperature', 'velocity_Magnitude']
#plotOverLine1Display_1.SeriesVisibility = ['temperature']
#plotOverLine1Display_1.SeriesLabel = ['arc_length', 'arc_length', 'pressure', 'pressure', 'rho', 'rho', 'temperature', 'temperature', 'velocity_X', 'velocity_X', 'velocity_Y', 'velocity_Y', 'velocity_Z', 'velocity_Z', 'velocity_Magnitude', 'velocity_Magnitude', 'vtkValidPointMask', 'vtkValidPointMask', 'Points_X', 'Points_X', 'Points_Y', 'Points_Y', 'Points_Z', 'Points_Z', 'Points_Magnitude', 'Points_Magnitude']
#plotOverLine1Display_1.SeriesColor = ['arc_length', '0', '0', '0', 'pressure', '0.89', '0.1', '0.11', 'rho', '0.22', '0.49', '0.72', 'temperature', '0.3', '0.69', '0.29', 'velocity_X', '0.6', '0.31', '0.64', 'velocity_Y', '1', '0.5', '0', 'velocity_Z', '0.65', '0.34', '0.16', 'velocity_Magnitude', '0', '0', '0', 'vtkValidPointMask', '0.89', '0.1', '0.11', 'Points_X', '0.22', '0.49', '0.72', 'Points_Y', '0.3', '0.69', '0.29', 'Points_Z', '0.6', '0.31', '0.64', 'Points_Magnitude', '1', '0.5', '0']
#plotOverLine1Display_1.SeriesPlotCorner = ['arc_length', '0', 'pressure', '0', 'rho', '0', 'temperature', '0', 'velocity_X', '0', 'velocity_Y', '0', 'velocity_Z', '0', 'velocity_Magnitude', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
#plotOverLine1Display_1.SeriesLabelPrefix = ''
#plotOverLine1Display_1.SeriesLineStyle = ['arc_length', '1', 'pressure', '1', 'rho', '1', 'temperature', '1', 'velocity_X', '1', 'velocity_Y', '1', 'velocity_Z', '1', 'velocity_Magnitude', '1', 'vtkValidPointMask', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'Points_Magnitude', '1']
#plotOverLine1Display_1.SeriesLineThickness = ['arc_length', '2', 'pressure', '2', 'rho', '2', 'temperature', '2', 'velocity_X', '2', 'velocity_Y', '2', 'velocity_Z', '2', 'velocity_Magnitude', '2', 'vtkValidPointMask', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'Points_Magnitude', '2']
#plotOverLine1Display_1.SeriesMarkerStyle = ['arc_length', '0', 'pressure', '0', 'rho', '0', 'temperature', '0', 'velocity_X', '0', 'velocity_Y', '0', 'velocity_Z', '0', 'velocity_Magnitude', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
#
## update the view to ensure updated data information
#lineChartView1.Update()
#
### hide data in view
##Hide(plotOverLine1, lineChartView1)


###############################################################################
# Display and/or take screenshots
###############################################################################

# Display the render view
SetActiveView(renderView1)

#### uncomment the following to render all views
RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).


