# trace generated using paraview version 5.11.2
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11
import sys
import os
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()



method_names = ['SFC_morton','parMETIS','fastPart','ptscotch']


file_names = [os.environ[method_name] for method_name in method_names]

# file_names = [sys.argv[1:]]
# print(file_names)
assert len(file_names) ==  len(method_names), f'{len(method_names)} should be provided' 

readers = [LegacyVTKReader(registrationName=method_names[i], FileNames=[file_names[i]]) for i in range(len(method_names))]

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# get layout
layout1 = GetLayout()

# split cell
layout1.SplitHorizontal(0, 0.5)


# split cell
layout1.SplitVertical(2, 0.5)


# split cell
layout1.SplitVertical(1, 0.5)


# Create a new 'Render View'
renderView2 = CreateView('RenderView')


# assign view to a particular cell in the layout
AssignViewToLayout(view=renderView2, layout=layout1, hint=1)


# Create a new 'Render View'
renderView3 = CreateView('RenderView')


# assign view to a particular cell in the layout
AssignViewToLayout(view=renderView3, layout=layout1, hint=2)

# Create a new 'Render View'
renderView4 = CreateView('RenderView')

# assign view to a particular cell in the layout
AssignViewToLayout(view=renderView4, layout=layout1, hint=3)


renderViews = [renderView1,renderView2,renderView3,renderView4]

renderTextLabels = [Text(Text=method_names[i],registrationName=f'Text{i}') for i in range(4)]

colorsLUT = GetColorTransferFunction('Colors')


for i in range(4):
    # show data in view
    display = Show(readers[i], renderViews[i], 'GeometryRepresentation')

    # trace defaults for the display properties.
    display.Representation = 'Surface'
    display.ColorArrayName = ['POINTS', 'Colors']
    display.LookupTable = colorsLUT
    display.SelectTCoordArray = 'None'
    display.SelectNormalArray = 'None'
    display.SelectTangentArray = 'None'
    display.OSPRayScaleArray = 'Colors'
    display.OSPRayScaleFunction = 'PiecewiseFunction'
    display.SelectOrientationVectors = 'Colors'
    display.SelectScaleArray = 'Colors'
    display.GlyphType = 'Arrow'
    display.GlyphTableIndexArray = 'Colors'
    display.SetScaleArray = ['POINTS', 'Colors']
    display.ScaleTransferFunction = 'PiecewiseFunction'
    display.OpacityArray = ['POINTS', 'Colors']
    display.OpacityTransferFunction = 'PiecewiseFunction'
    display.DataAxesGrid = 'GridAxesRepresentation'
    display.PolarAxes = 'PolarAxesRepresentation'
    display.SelectInputVectors = ['POINTS', 'Colors']
    display.WriteLog = ''

    display.MapScalars = 0


    # show color bar/color legend
    display.SetScalarBarVisibility(renderViews[i], False)

    Show(renderTextLabels[i], renderViews[i], 'TextSourceRepresentation')



# reset view to fit data
renderView1.ResetCamera(False)          # TODO: do not remove



# link cameras in two views
AddCameraLink(renderView1, renderView2, 'CameraLink0')



# link cameras in two views
AddCameraLink(renderView4, renderView2, 'CameraLink1')



# link cameras in two views
AddCameraLink(renderView3, renderView2, 'CameraLink2')



#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(1396, 734)


#--------------------------------------------
# uncomment the following to render all views
RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).