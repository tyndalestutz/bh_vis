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
# bh1, bh2, are objects that hold the xyz coordinates of each black hole
# bhL1, bhL2, are the objects that hold the xyz magnitudes of angular momentum (L) vectors 1 and 2.


###File names
#Note: fix this db path
bh_database = "/home/guest/Documents/BH_Vis/data/mesh/spheres/iscos_sphere_sub8.obj"
bh_data = "/home/guest/Documents/BH_Vis/data/synthetic_coords/synthetic_data_ang_momentum.csv"
gw_database = "/home/guest/Documents/BH_Vis/data/mesh/gw_test/state*.vts database"

'''
#Background Image Citation: ESA/Hubble & NASA, https://www.nasa.gov/image-feature/goddard/2016/hubble-spots-an-irregular-island-in-a-sea-of-space
background_image = parent_directory + "/data/background_images/blue_grid.png"
frame_name = "synthetic_BH_test_animation"
movie_output_destination = parent_directory + "/movies/movie2"
'''

## object holding global parameters for the visulization settings, currently only background and camera settings
# to do: add more global parameters in this object (lighting, colors, etc.)
class VisSettings:
    green_screen = False
    show_axis_annotations = False

    # default camera view is top down
    class Camera:
        x_normal = 0 
        y_normal = 0
        z_normal = 1
        x_up = 0 
        z_up = 0
        y_up = 1
        zoom = 1
    cam = Camera()

# 2 different settings for testing purposes, the default is just settings
settings = VisSettings()
alt_settings = VisSettings()
alt_settings.green_screen = True
alt_settings.cam.x_normal = 1
alt_settings.cam.z_normal = 0.1
alt_settings.cam.x_up = -0.1
alt_settings.cam.y_up = 0.25
alt_settings.cam.z_up = 1
alt_settings.cam.zoom = 1.5

# change this to the setting configuration you want
active_setting = settings

# axis/boxes, background, camera, lighting, save window file settings
def default_atts():
    a = v.AnnotationAttributes()

    # determines if markers for axis, bounding box, etc. are visible
    if active_setting.show_axis_annotations == False:
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

    # creates background
    if active_setting.green_screen:
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
    
    # sets camera orientation
    View3DAtts = v.View3DAttributes()
    View3DAtts.viewNormal = (active_setting.cam.x_normal, active_setting.cam.y_normal, active_setting.cam.z_normal)
    View3DAtts.focus = (0, 0, 0)
    View3DAtts.viewUp = (active_setting.cam.x_up, active_setting.cam.y_up, active_setting.cam.z_up)
    View3DAtts.nearPlane = -28.4429
    View3DAtts.farPlane = 28.4429
    View3DAtts.imageZoom = active_setting.cam.zoom
    View3DAtts.windowValid = 1
    v.SetView3D(View3DAtts)

    # sets lighting settings
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

    # defines file settings for when we save a window
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

# updates the position of the black holes
def set_coords_mesh(objNum, x, y, z):
    v.SetActivePlots(objNum)
    trasnformAtts = v.TransformAttributes()
    trasnformAtts.doTranslate = 1
    trasnformAtts.translateX = x
    trasnformAtts.translateY = y
    trasnformAtts.translateZ = z
    v.SetOperatorOptions(trasnformAtts, 0, 0)

def create_BH_vectors():
    global bh1_vector, bh2_vector
    bh1_vector = v.CreateAnnotationObject("Line3D")
    bh2_vector = v.CreateAnnotationObject("Line3D")

    #to find annotation attributes
    #print(bh1_vector) 

    bh1_vector.useForegroundForLineColor = 0
    bh1_vector.color = (255, 255, 255, 255)
    bh1_vector.width = 3
    bh1_vector.arrow2 = 1
    bh1_vector.arrow2Height = 0.2
    bh1_vector.arrow2Radius = 0.075

    bh2_vector.useForegroundForLineColor = 0
    bh2_vector.color = (255, 255, 255, 255)
    bh2_vector.width = 3
    bh2_vector.arrow2 = 1
    bh2_vector.arrow2Height = 0.2
    bh2_vector.arrow2Radius = 0.075

default_atts()
create_spheres()
create_gw()
create_BH_vectors()

# general coordinate object - all method arguments are in the order of x, y, z coordinates
class CartesianCoords:
    def __init__(self, x_coord, y_coord, z_coord):
        self.x = x_coord
        self.y = y_coord
        self.z = z_coord
    def set_coords(self, set_x, set_y, set_z):
        self.x = float(set_x)
        self.y = float(set_y)
        self.z = float(set_z)

# creates coordinates for black holes and their vectors
bh1 = CartesianCoords(0, 0, 0)
bh2 = CartesianCoords(0, 0, 0)
bhL1 = CartesianCoords(0, 0, 0)
bhL2 = CartesianCoords(0, 0, 0)

with open(bh_data, 'r') as file:
    reader = csv.reader(file)
    header = next(reader)  # skip the header row
    v.DrawPlots()
    for row in reader:
        # moves to the next gw vtk file
        v.TimeSliderNextState()
        
        t = float(row[0])

        # updates positions of black holes and corresponding vectors
        bh1.set_coords(row[1], row[2], 2)
        bh2.set_coords(row[7], row[8], 2)
        bhL1.set_coords(row[4], row[5], row[6])
        bhL2.set_coords(row[10], row[11], row[12])

        # sets updated positions
        bh1_vector.point1 = (bh1.x, bh1.y, bh1.z)
        bh1_vector.point2 = (bh1.x + bhL1.x, bh1.y + bhL1.y, bh1.z + bhL1.z)
        bh2_vector.point1 = (bh2.x, bh2.y, bh2.z)
        bh2_vector.point2 = (bh2.x + bhL2.x, bh2.y + bhL2.y, bh2.z + bhL2.z)
        set_coords_mesh(0, bh1.x, bh1.y, bh1.z)
        set_coords_mesh(1, bh2.x, bh2.y, bh2.z)
        
        v.DrawPlots()

        # sets save window settings and saves the window to png
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

