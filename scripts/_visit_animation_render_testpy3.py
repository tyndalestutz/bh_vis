import sys, csv, os
sys.path.append("/usr/local/3.3.1/linux-x86_64/lib/site-packages") 
import visit
visit.Launch()
import visit
v = visit

parent_directory = os.path.dirname(os.path.dirname(__file__))

# In order for the render to record, 'record_render' must be True
record_render = False


# object 0: Black hole 1
# object 1: Black hole 2
# bh1_x, bh2_x, etc are the values of the xyz coordinates of each black hole
# L1_X, L2_x, etc are the xyz magnitudes of angular momentum (L) vectors 1 and 2.


###File names

#Note: fix this db path
bh_database = "/home/guest/Documents/BH_Vis/data/mesh/spheres/iscos_sphere_sub8.obj"
bh_data = "/home/guest/Documents/BH_Vis/data/synthetic_coords/synthetic_data_ang_momentum.csv"
gw_database = "/home/guest/Documents/BH_Vis/data/mesh/synth_gw/test15/state*.vts database"
print(gw_database)
'''
#Background Image Citation: ESA/Hubble & NASA, https://www.nasa.gov/image-feature/goddard/2016/hubble-spots-an-irregular-island-in-a-sea-of-space
background_image = parent_directory + "/data/background_images/blue_grid.png"
frame_name = "synthetic_BH_test_animation"
movie_output_destination = parent_directory + "/movies/movie2"
'''

###Parameters
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

    #this doesn't work for some reason
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
    #set attributes
    PseudocolorAtts = v.PseudocolorAttributes()
    PseudocolorAtts.minFlag = 1
    PseudocolorAtts.min = -0.1
    PseudocolorAtts.useBelowMinColor = 1
    PseudocolorAtts.belowMinColor = (31, 31, 31, 255)
    PseudocolorAtts.maxFlag = 1
    PseudocolorAtts.max = 0.1
    PseudocolorAtts.useAboveMaxColor = 1
    PseudocolorAtts.aboveMaxColor = (0, 255, 0, 255)

    v.AddPlot("Pseudocolor", "mesh_quality/warpage", 1, 1)
    v.AddOperator("Transform")
    v.SetPlotOptions(PseudocolorAtts)
    v.AddPlot("Pseudocolor", "mesh_quality/warpage", 1, 1)
    v.AddOperator("Transform")
    v.SetPlotOptions(PseudocolorAtts)

def create_gw():
    visit.OpenDatabase(gw_database)

    #set attributes
    PseudocolorAtts = v.PseudocolorAttributes()
    PseudocolorAtts.minFlag = 1
    PseudocolorAtts.min = -0.1
    PseudocolorAtts.useBelowMinColor = 1
    PseudocolorAtts.belowMinColor = (31, 31, 31, 255)
    PseudocolorAtts.maxFlag = 1
    PseudocolorAtts.max = 0.1
    PseudocolorAtts.useAboveMaxColor = 1
    PseudocolorAtts.aboveMaxColor = (67, 91, 122, 255)

    v.AddPlot("Mesh", "mesh", 1, 1)
    v.AddPlot("Pseudocolor", "mesh_quality/degree", 1, 1)
    v.SetPlotOptions(PseudocolorAtts)

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
create_gw()

L1 = v.CreateAnnotationObject("Line3D")
L2 = v.CreateAnnotationObject("Line3D")

#to find annotation attributes
#print(L1) 

L1.useForegroundForLineColor = 0
L1.color = (255, 255, 255, 255)
L1.width = 3
L1.arrow2 = 1
L1.arrow2Height = 0.2
L1.arrow2Radius = 0.075

L2.useForegroundForLineColor = 0
L2.color = (255, 255, 255, 255)
L2.width = 3
L2.arrow2 = 1
L2.arrow2Height = 0.2
L2.arrow2Radius = 0.075


with open(bh_data, 'r') as file:
    reader = csv.reader(file)
    header = next(reader)  # skip the header row
    v.DrawPlots()
    for row in reader:
        # moves to the next gw vtk file
        v.TimeSliderNextState()
        
        t = float(row[0])
        bh1_x = float(row[1])
        bh1_y = float(row[2])
        bh1_z = 2 #float(row[3])
        L1_x = float(row[4])
        L1_y = float(row[5])
        L1_z = float(row[6])
        bh2_x = float(row[7])
        bh2_y = float(row[8])
        bh2_z = 2 #float(row[9])
        L2_x = float(row[10])
        L2_y = float(row[11])
        L2_z = float(row[12])
        
        L1.point1 = (bh1_x, bh1_y, bh1_z)
        L1.point2 = (bh1_x + L1_x, bh1_y + L1_y, bh1_z + L1_z)
        L2.point1 = (bh2_x, bh2_y, bh2_z)
        L2.point2 = (bh2_x + L2_x, bh2_y + L2_y, bh2_z + L2_z)
        set_coords(0, bh1_x, bh1_y, bh1_z)
        set_coords(1, bh2_x, bh2_y, bh2_z)
        v.DrawPlots()
        s = v.SaveWindowAttributes()
        s.fileName = "synthetic_BH_test_animation"
        s.format = s.PNG
        s.progressive = 1
        s.width = 772
        s.height = 702
        s.screenCapture = 1
        v.SetSaveWindowAttributes(s)
        # Save the window
        #v.SaveWindow()

#frame_name_pattern = "synthetic_BH_test_animation_%04d.png"
#movie_name = "streamline_crop_example.mp4"
# The following command worked to create a movie from the file of pngs
#ffmpeg -framerate 1000 -i synthetic_BH_test_animation%04d.png -vf "scale=770:-2" -c:v libx264 -r 1000 -pix_fmt yuv420p output.mp4

