# This is currently the primary python script for animating a
# bh merger with polar gravitational wave meshes. 

import sys, csv, os
sys.path.append("/usr/local/3.3.1/linux-x86_64/lib/site-packages") 
import visit
visit.Launch()
import visit
v = visit

parent_directory = os.path.dirname(os.path.dirname(__file__))

# In order for the render to record, 'record_render' must be True
record_render = False

bh_database = "/home/guest/data/horizon_data_obj/h.t277632.ah2_state0.obj"
bh_data = "/home/guest/Documents/bh_vis/scripts/sorted.csv"
gw_database = "/home/guest/Documents/BH_Vis_local/data/mesh/gw_test8_polar_zeroR_r=200/state*.vtu database"

# Parameters
green_screen = False
show_axis_annotations = False

def default_atts():
    a = v.AnnotationAttributes()
    if show_axis_annotations == False:
        a.axes2D.visible = 0
        a.axes3D.visible = 0
        a.axes3D.triadFlag = 0
        a.axes3D.bboxFlag = 0
        a.axes3D.xAxis.grid = 0
        a.axes3D.yAxis.grid = 0
        a.axes3D.zAxis.grid = 0
        a.axes3D.triadLineWidth = 0
    a.userInfoFlag = 0
    a.databaseInfoFlag = 0
    a.legendInfoFlag = 0
    a.axesArray.lineWidth = 0
    a.axesArray.axes.grid = 0

    if green_screen:
        a.backgroundColor = (0, 255, 0, 255)
        a.foregroundColor = (255, 255, 255, 255)
        a.backgroundMode = a.Solid  # Solid, Gradient, Image, ImageSphere
    else:
        a.gradientBackgroundStyle = a.Image
        a.gradientColor1 = (80, 80, 80, 255)
        a.gradientColor2 = (70, 70, 70, 255)
        #a.backgroundImage = background_image
        a.backgroundMode = a.Gradient  # Solid, Gradient, Image, ImageSphere
        a.imageRepeatX = 1
        a.imageRepeatY = 1
    v.SetAnnotationAttributes(a)
    
    View3DAtts = v.View3DAttributes()
    View3DAtts.viewNormal = (0, 0, 1)
    View3DAtts.focus = (0, 0, 0)
    View3DAtts.viewUp = (0, 1, 0)
    View3DAtts.nearPlane = -28.4429
    View3DAtts.farPlane = 28.4429
    View3DAtts.windowValid = 1
    v.SetView3D(View3DAtts)

    light = v.LightAttributes()
    light.enabledFlag = 1
    light.type = light.Object 
    light.direction = (0.666, -0.666, -0.666)
    v.SetLight(0, light)
        
    RenderingAtts = v.RenderingAttributes()
    RenderingAtts.specularFlag = 1
    RenderingAtts.specularCoeff = 0.13
    RenderingAtts.specularPower = 3.2
    RenderingAtts.specularColor = (255, 255, 255, 255)
    v.SetRenderingAttributes(RenderingAtts)

    s = v.SaveWindowAttributes()
    s.outputToCurrentDirectory = 0
    #s.outputDirectory = movie_output_destination
    #s.fileName = frame_name
    s.format = s.PNG
    s.progressive = 1
    s.width = 772
    s.height = 702
    s.screenCapture = 1
    v.SetSaveWindowAttributes(s)

def create_spheres():
    visit.OpenDatabase(bh_database)
    v.PseudocolorAtts = v.PseudocolorAttributes()
    v.PseudocolorAtts.minFlag = 1
    v.PseudocolorAtts.min = -0.1
    v.PseudocolorAtts.useBelowMinColor = 1
    v.PseudocolorAtts.belowMinColor = (31, 31, 31, 255)
    v.PseudocolorAtts.maxFlag = 1
    v.PseudocolorAtts.max = 0.1
    v.PseudocolorAtts.useAboveMaxColor = 1
    v.PseudocolorAtts.aboveMaxColor = (0, 255, 0, 255)
    v.AddPlot("Pseudocolor", "mesh_quality/warpage", 1, 1)
    v.AddOperator("Transform")
    v.SetPlotOptions(v.PseudocolorAtts)
    v.AddPlot("Pseudocolor", "mesh_quality/warpage", 1, 1)
    v.AddOperator("Transform")
    v.SetPlotOptions(v.PseudocolorAtts)

def create_gw():
    visit.OpenDatabase(gw_database)
    v.AddPlot("Pseudocolor", "Strain", 1, 1)
    PseudocolorAtts = v.PseudocolorAttributes()
    PseudocolorAtts.scaling = PseudocolorAtts.Linear  # Linear, Log, Skew
    PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData  # OriginalData, ActualData
    PseudocolorAtts.min = -10
    PseudocolorAtts.maxFlag = 1
    PseudocolorAtts.max = 10
    PseudocolorAtts.centering = PseudocolorAtts.Nodal  # Natural, Nodal, Zonal
    PseudocolorAtts.colorTableName = "viridis_light" # imola-seq and viridis_light are nice
    PseudocolorAtts.invertColorTable = 0
    PseudocolorAtts.opacityType = PseudocolorAtts.FullyOpaque  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
    PseudocolorAtts.pointType = PseudocolorAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
    PseudocolorAtts.pointSizeVarEnabled = 0
    PseudocolorAtts.pointSizeVar = "default"
    PseudocolorAtts.pointSizePixels = 2
    PseudocolorAtts.lineType = PseudocolorAtts.Line  # Line, Tube, Ribbon
    PseudocolorAtts.lineWidth = 0
    PseudocolorAtts.renderSurfaces = 1
    PseudocolorAtts.renderWireframe = 0
    PseudocolorAtts.renderPoints = 0
    PseudocolorAtts.smoothingLevel = 0
    PseudocolorAtts.legendFlag = 1
    PseudocolorAtts.lightingFlag = 1
    v.SetPlotOptions(PseudocolorAtts)
    v.DrawPlots()

def set_coords(objNum, x, y, z):
    v.SetActivePlots(objNum)
    trasnformAtts = v.TransformAttributes()
    trasnformAtts.doTranslate = 1
    trasnformAtts.translateX = x
    trasnformAtts.translateY = y
    trasnformAtts.translateZ = z
    v.SetOperatorOptions(trasnformAtts, 0, 0)

default_atts()
create_spheres()
#create_gw()

with open(bh_data, 'r') as file:
    reader = csv.reader(file)
    header = next(reader)  # skip the header row
    v.DrawPlots()
    for row in reader:
        # moves to the next gw vtk file
        v.TimeSliderNextState()
        t = float(row[0])
        bh1_x = float(row[1])
        bh1_y = float(row[3])
        bh1_z = 2 #float(row[3])
        bh2_x = float(row[2])
        bh2_y = float(row[4])
        bh2_z = 2 #float(row[9])
        set_coords(0, bh1_x, bh1_y, bh1_z)
        set_coords(1, bh2_x, bh2_y, bh2_z)
        v.DrawPlots()
        s = v.SaveWindowAttributes()
        s.fileName = "BH_test_animation"
        s.format = s.PNG
        s.progressive = 1
        s.width = 772
        s.height = 702
        s.screenCapture = 1
        v.SetSaveWindowAttributes(s)
        # Save the window
        #v.SaveWindow()

# frame_name_pattern = "synthetic_BH_test_animation_%04d.png"
# movie_name = "streamline_crop_example.mp4"
# The following command worked to create a movie from the file of pngs
# ffmpeg -framerate 1000 -i synthetic_BH_test_animation%04d.png -vf "scale=770:-2" -c:v libx264 -r 1000 -pix_fmt yuv420p output.mp4