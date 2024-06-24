# trace generated using paraview version 5.11.2
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11
import os

#### import the simple module from the paraview
from paraview.simple import *

pagerank_vtk_file_name =os.environ['pagerank_vtk_file_name']

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'Legacy VTK Reader'
pagerankvtk = LegacyVTKReader(registrationName='pagerank', FileNames=[pagerank_vtk_file_name])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
pagerankvtkDisplay = Show(pagerankvtk, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'pagerank'
pagerankLUT = GetColorTransferFunction('pagerank')

# trace defaults for the display properties.
pagerankvtkDisplay.Representation = 'Surface'
pagerankvtkDisplay.ColorArrayName = ['POINTS', 'pagerank']
pagerankvtkDisplay.LookupTable = pagerankLUT
pagerankvtkDisplay.SelectTCoordArray = 'None'
pagerankvtkDisplay.SelectNormalArray = 'None'
pagerankvtkDisplay.SelectTangentArray = 'None'
pagerankvtkDisplay.OSPRayScaleArray = 'pagerank'
pagerankvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
pagerankvtkDisplay.SelectOrientationVectors = 'None'
pagerankvtkDisplay.ScaleFactor = 5.945670783519745
pagerankvtkDisplay.SelectScaleArray = 'pagerank'
pagerankvtkDisplay.GlyphType = 'Arrow'
pagerankvtkDisplay.GlyphTableIndexArray = 'pagerank'
pagerankvtkDisplay.GaussianRadius = 0.29728353917598727
pagerankvtkDisplay.SetScaleArray = ['POINTS', 'pagerank']
pagerankvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
pagerankvtkDisplay.OpacityArray = ['POINTS', 'pagerank']
pagerankvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
pagerankvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
pagerankvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
pagerankvtkDisplay.SelectInputVectors = [None, '']
pagerankvtkDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
pagerankvtkDisplay.ScaleTransferFunction.Points = [306.0, 0.0, 0.5, 0.0, 393.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
pagerankvtkDisplay.OpacityTransferFunction.Points = [306.0, 0.0, 0.5, 0.0, 393.0, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera(False)

# get the material library
materialLibrary1 = GetMaterialLibrary()

# show color bar/color legend
pagerankvtkDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get opacity transfer function/opacity map for 'pagerank'
pagerankPWF = GetOpacityTransferFunction('pagerank')

# get 2D transfer function for 'pagerank'
pagerankTF2D = GetTransferFunction2D('pagerank')

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
pagerankLUT.ApplyPreset('Cold and Hot', True)


# rescale color and/or opacity maps used to exactly fit the current data range
pagerankvtkDisplay.RescaleTransferFunctionToDataRange(False, True)

pagerankLUT = GetColorTransferFunction('pagerank')
pagerankLUT.ShowDataHistogram = 1


# Properties modified on pagerankLUT
pagerankLUT.AutomaticDataHistogramComputation = 1

# Properties modified on pagerankLUT
pagerankLUT.DataHistogramNumberOfBins = 50

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

# get layout
layout1 = GetLayout()

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(1103, 622)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
# renderView1.CameraPosition = [22.392648462206125, 30.10434490442276, 211.45100651253964]
# renderView1.CameraFocalPoint = [22.392648462206125, 30.10434490442276, 29.176661547273397]
# renderView1.CameraParallelScale = 47.17607191059766

#--------------------------------------------
# uncomment the following to render all views
RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).