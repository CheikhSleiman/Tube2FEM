# trace generated using paraview version 5.10.0-RC1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'PVD Reader'
velocitypvd = PVDReader(registrationName='velocity.pvd', FileName=velPath)
velocitypvd.PointArrays = ['f_29']

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
velocitypvdDisplay = Show(velocitypvd, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
velocitypvdDisplay.Representation = 'Surface'
velocitypvdDisplay.ColorArrayName = [None, '']
velocitypvdDisplay.SelectTCoordArray = 'None'
velocitypvdDisplay.SelectNormalArray = 'None'
velocitypvdDisplay.SelectTangentArray = 'None'
velocitypvdDisplay.OSPRayScaleArray = 'f_29'
velocitypvdDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
velocitypvdDisplay.SelectOrientationVectors = 'f_29'
velocitypvdDisplay.ScaleFactor = 0.00793351722098945
velocitypvdDisplay.SelectScaleArray = 'None'
velocitypvdDisplay.GlyphType = 'Arrow'
velocitypvdDisplay.GlyphTableIndexArray = 'None'
velocitypvdDisplay.GaussianRadius = 0.00039667586104947253
velocitypvdDisplay.SetScaleArray = ['POINTS', 'f_29']
velocitypvdDisplay.ScaleTransferFunction = 'PiecewiseFunction'
velocitypvdDisplay.OpacityArray = ['POINTS', 'f_29']
velocitypvdDisplay.OpacityTransferFunction = 'PiecewiseFunction'
velocitypvdDisplay.DataAxesGrid = 'GridAxesRepresentation'
velocitypvdDisplay.PolarAxes = 'PolarAxesRepresentation'
velocitypvdDisplay.ScalarOpacityUnitDistance = 0.0010970428225556563
velocitypvdDisplay.OpacityArrayName = ['POINTS', 'f_29']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
velocitypvdDisplay.ScaleTransferFunction.Points = [-0.001194082437291646, 0.0, 0.5, 0.0, -7.743731586527737e-05, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
velocitypvdDisplay.OpacityTransferFunction.Points = [-0.001194082437291646, 0.0, 0.5, 0.0, -7.743731586527737e-05, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera(False)

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(velocitypvdDisplay, ('POINTS', 'f_29', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
velocitypvdDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
velocitypvdDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'f_29'
f_29LUT = GetColorTransferFunction('f_29')

# get opacity transfer function/opacity map for 'f_29'
f_29PWF = GetOpacityTransferFunction('f_29')

animationScene1.GoToLast()

# rescale color and/or opacity maps used to exactly fit the current data range
velocitypvdDisplay.RescaleTransferFunctionToDataRange(False, True)

# create a new 'Glyph'
glyph1 = Glyph(registrationName='Glyph1', Input=velocitypvd,
    GlyphType='Arrow')
glyph1.OrientationArray = ['POINTS', 'f_29']
glyph1.ScaleArray = ['POINTS', 'No scale array']
glyph1.ScaleFactor = 0.00793351722098945
glyph1.GlyphTransform = 'Transform2'

# show data in view
glyph1Display = Show(glyph1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.ColorArrayName = ['POINTS', 'f_29']
glyph1Display.LookupTable = f_29LUT
glyph1Display.SelectTCoordArray = 'None'
glyph1Display.SelectNormalArray = 'None'
glyph1Display.SelectTangentArray = 'None'
glyph1Display.OSPRayScaleArray = 'f_29'
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.SelectOrientationVectors = 'f_29'
glyph1Display.ScaleFactor = 0.00917029562406242
glyph1Display.SelectScaleArray = 'None'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GlyphTableIndexArray = 'None'
glyph1Display.GaussianRadius = 0.00045851478120312094
glyph1Display.SetScaleArray = ['POINTS', 'f_29']
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityArray = ['POINTS', 'f_29']
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
glyph1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
glyph1Display.ScaleTransferFunction.Points = [-0.10119822714411862, 0.0, 0.5, 0.0, 5.668172509006031e-05, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
glyph1Display.OpacityTransferFunction.Points = [-0.10119822714411862, 0.0, 0.5, 0.0, 5.668172509006031e-05, 1.0, 0.5, 0.0]

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on glyph1
glyph1.ScaleArray = ['POINTS', 'f_29']

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(velocitypvd)

# Properties modified on velocitypvdDisplay
velocitypvdDisplay.Opacity = 0.9

# Properties modified on velocitypvdDisplay
velocitypvdDisplay.Opacity = 0.8

# Properties modified on velocitypvdDisplay
velocitypvdDisplay.Opacity = 0.7

# Properties modified on velocitypvdDisplay
velocitypvdDisplay.Opacity = 0.6

# Properties modified on velocitypvdDisplay
velocitypvdDisplay.Opacity = 0.5

# Properties modified on velocitypvdDisplay
velocitypvdDisplay.Opacity = 0.4

# Properties modified on velocitypvdDisplay
velocitypvdDisplay.Opacity = 0.3

# Properties modified on velocitypvdDisplay
velocitypvdDisplay.Opacity = 0.2

# set active source
SetActiveSource(glyph1)

# Properties modified on glyph1
glyph1.ScaleArray = ['POINTS', 'No scale array']

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on glyph1
glyph1.ScaleArray = ['POINTS', 'f_29']

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on glyph1
glyph1.ScaleFactor = 0.0793351722098945

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on glyph1
glyph1.ScaleFactor = 0.1

# update the view to ensure updated data information
renderView1.Update()

animationScene1.Play()

# get layout
layout1 = GetLayout()

# layout/tab size in pixels
layout1.SetSize(1284, 540)

# current camera placement for renderView1
renderView1.CameraPosition = [-0.1451826147392004, -0.13583695211484575, 0.04862766592187073]
renderView1.CameraFocalPoint = [-0.03999973846740475, -8.364999999999762e-06, 0.0055552]
renderView1.CameraViewUp = [0.1220154933095204, 0.2129254685452332, 0.9694199112032037]
renderView1.CameraParallelScale = 0.04583953442809416

# save animation
SaveAnimation(postPath, renderView1, ImageResolution=[1284, 540],
    FrameWindow=[0, 9])

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(1304, 540)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.CameraPosition = [-0.1451826147392004, -0.13583695211484575, 0.04862766592187073]
renderView1.CameraFocalPoint = [-0.03999973846740475, -8.364999999999762e-06, 0.0055552]
renderView1.CameraViewUp = [0.1220154933095204, 0.2129254685452332, 0.9694199112032037]
renderView1.CameraParallelScale = 0.04583953442809416

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).