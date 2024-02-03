"""

This script takes in a black hole mesh (i.e. a sphere mesh), a gravitational wave mesh,
and a csv file containing the coordinates of the black holes and the angular momenta
of the black holes. It starts by initializing VisIt's (very many) attributes, creates
the black holes, gravitational wave, and angular momentum vector objects, then
animates the collision. If record_render is set to True, it will save the frames of
the animation to the movie_output_destination directory, which will then be combined
into a movie.

"""

import sys
import csv
import os

sys.path.append("/usr/local/3.3.1/linux-x86_64/lib/site-packages")
import visit

visit.Launch()
import visit

v = visit


# Main function
def main():
    ##### File names #####
    parent_directory = os.path.dirname(os.path.dirname(__file__))
    bh_database = os.path.join(
        parent_directory, "..", "data", "dummy", "h.t277632.ah2_state0.obj"
    )
    gw_database = os.path.join(
        parent_directory, "..", "data", "dummy", "gw_dummy", "state*.vtu database"
    )
    bh_data = os.path.join(parent_directory, "..", "data", "dummy", "bh_xyz_dummy.csv")
    movie_output_destination = os.path.join(
        parent_directory, "..", "/movies/dummy"
    )

    ##### Parameters #####
    record_render = False

    ##### Initialize Objects #####
    default_atts(movie_output_destination)
    create_sphere(bh_database)  # Black hole 1
    create_sphere(bh_database)  # Black hole 2
    create_gw(gw_database)  # Gravitational wave
    L1 = create_vector()
    L2 = create_vector()

    ##### Animation #####
    with open(bh_data, "r") as file:
        reader = csv.reader(file)
        _ = next(reader)  # skip the header row
        v.DrawPlots()
        for row in reader:
            # Move to the next gw vtk file
            v.TimeSliderNextState()

            # DEBUG maybe make this look a little nicer using tuples,
            # don't want to mess with it until I can run it
            # Set black hole xyz coordinates and the xyz magnitues of the angular momenta
            _ = float(row[0]) # time
            bh1_x = float(row[1])
            bh1_y = float(row[2])
            bh1_z = float(row[3])
            L1_x = float(row[4])
            L1_y = float(row[5])
            L1_z = float(row[6])
            bh2_x = float(row[7])
            bh2_y = float(row[8])
            bh2_z = float(row[9])
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

            # Save the window
            if record_render:
                s = movie_atts()
                v.SetSaveWindowAttributes(s)
                v.SaveWindow()

    # This is bash code that will combine the frames into a movie,
    # lets incorporate it into the script
    # Either by calling it from the script or by writing it in python

    # frame_name_pattern = "synthetic_BH_test_animation_%04d.png"
    # movie_name = "streamline_crop_example.mp4"
    # ffmpeg -framerate 1000 -i synthetic_BH_test_animation%04d.png -vf "scale=770:-2" -c:v libx264 -r 1000 -pix_fmt yuv420p output.mp4


# Sets the default attributes for the visualization
def default_atts(movie_output_destination):
    ##### Annotation Parameters #####
    a = v.AnnotationAttributes()
    show_axis_annotations = False
    if not show_axis_annotations:
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

    # This doesn't work for some reason
    v.SetAnnotationAttributes(a)

    ##### View Parameters #####
    viewing_atts = v.View3DAttributes()
    viewing_atts.viewNormal = (0, 0, 1)
    viewing_atts.focus = (0, 0, 0)
    viewing_atts.viewUp = (0, 1, 0)
    viewing_atts.nearPlane = -28.4429
    viewing_atts.farPlane = 28.4429
    viewing_atts.windowValid = 1
    v.SetView3D(viewing_atts)

    light = v.LightAttributes()
    light.enabledFlag = 1
    light.type = light.Object
    light.direction = (0.666, -0.666, -0.666)
    v.SetLight(0, light)

    rendering_atts = v.RenderingAttributes()
    rendering_atts.specularFlag = 1
    rendering_atts.specularCoeff = 0.13
    rendering_atts.specularPower = 3.2
    rendering_atts.specularColor = (255, 255, 255, 255)
    v.SetRenderingAttributes(rendering_atts)

    ##### Recording Parameters #####
    green_screen = False
    # background_image = os.path.join(parent_directory, )

    if green_screen:
        a.backgroundColor = (0, 255, 0, 255)
        a.foregroundColor = (255, 255, 255, 255)
        a.backgroundMode = a.Solid  # Solid, Gradient, Image, ImageSphere
    else:
        a.gradientBackgroundStyle = a.Image
        a.gradientColor1 = (80, 80, 80, 255)
        a.gradientColor2 = (70, 70, 70, 255)
        # a.backgroundImage = background_image
        a.backgroundMode = a.Gradient  # Solid, Gradient, Image, ImageSphere
        a.imageRepeatX = 1
        a.imageRepeatY = 1

    s = v.SaveWindowAttributes()
    s.outputToCurrentDirectory = 0
    s.outputDirectory = movie_output_destination
    s.fileName = "BH_test_animation%04d.png"
    s.format = s.PNG
    s.progressive = 1
    s.width = 772
    s.height = 702
    s.screenCapture = 1
    v.SetSaveWindowAttributes(s)


# Creates a black hole sphere object
def create_sphere(bh_database):
    visit.OpenDatabase(bh_database)
    pseudocolor_atts = v.PseudocolorAttributes()
    pseudocolor_atts.minFlag = 1
    pseudocolor_atts.min = -0.1
    pseudocolor_atts.useBelowMinColor = 1
    pseudocolor_atts.belowMinColor = (31, 31, 31, 255)
    pseudocolor_atts.maxFlag = 1
    pseudocolor_atts.max = 0.1
    pseudocolor_atts.useAboveMaxColor = 1
    pseudocolor_atts.aboveMaxColor = (0, 255, 0, 255)

    v.AddPlot("Pseudocolor", "mesh_quality/warpage", 1, 1)
    v.AddOperator("Transform")
    v.SetPlotOptions(pseudocolor_atts)


# Creates a gravitational wave mesh object
def create_gw(gw_database):
    visit.OpenDatabase(gw_database)

    # set attributes
    pseudocolor_atts = v.PseudocolorAttributes()
    pseudocolor_atts.minFlag = 1
    pseudocolor_atts.min = -0.1
    pseudocolor_atts.useBelowMinColor = 1
    pseudocolor_atts.belowMinColor = (0, 0, 255, 255)
    pseudocolor_atts.maxFlag = 1
    pseudocolor_atts.max = 0.1
    pseudocolor_atts.useAboveMaxColor = 1
    pseudocolor_atts.aboveMaxColor = (0, 0, 255, 255)
    
    v.AddPlot("Mesh", "mesh", 1, 1)
    v.AddPlot("Pseudocolor", "mesh_quality/degree", 1, 1)
    v.SetPlotOptions(pseudocolor_atts)


# Creates a vector annotation
def create_vector(color=(255, 255, 255, 255)):
    vec = v.CreateAnnotationObject("Line3D")
    vec.useForegroundForLineColor = 0
    vec.color = color
    vec.width = 3
    vec.arrow2 = 1
    vec.arrow2Height = 0.2
    vec.arrow2Radius = 0.075

    return vec


# Sets the coordinates of object num obj_num
def set_coords(obj_num, x, y, z):
    v.SetActivePlots(obj_num)
    trasnform_atts = v.TransformAttributes()
    trasnform_atts.doTranslate = 1
    trasnform_atts.translateX = x
    trasnform_atts.translateY = y
    trasnform_atts.translateZ = z
    v.SetOperatorOptions(trasnform_atts, 0, 0)


# Sets the attributes for the movie
def movie_atts():
    s = v.SaveWindowAttributes()
    s.fileName = "frame"
    s.format = s.PNG
    s.progressive = 1
    s.width = 772
    s.height = 702
    s.screenCapture = 1

    return s


main()
