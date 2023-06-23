import sys, csv, os
sys.path.append("/usr/local/3.3.1/linux-x86_64/lib/site-packages") 
import visit
visit.Launch()
import visit
v = visit

parent_directory = os.path.dirname(os.path.dirname(__file__))
bh_database = "/home/guest/Documents/BH_Vis_local/data/mesh/spheres/iscos_sphere_sub8.obj"
bh_data = "/home/guest/Documents/BH_Vis_local/data/synthetic_coords/synthetic_data_ang_momentum.csv"
gw_database = "/home/guest/Documents/BH_Vis_local/data/mesh/gw_test8_polar_zeroR_r=200/state*.vtu database"

def EvalCubicSpline(t, allX, allY):
    n = len(allY)
    if((allX[0] > t) or (allX[n-1] < t)):
        raise Exception('t must be in the range between the first and last X')
    for i in range(1, n):
        if(allX[i] >= t):
            break
    i1 = max(i-2, 0)
    i2 = max(i-1, 0)
    i3 = i
    i4 = min(i+1, n-1)
    X = (allX[i1], allX[i2], allX[i3], allX[i4])
    Y = (allY[i1], allY[i2], allY[i3], allY[i4])
    dx = (X[2] - X[1])
    invdx = 1. / dx
    dy1   = (Y[2] + (Y[0] * -1.)) * (1. / (X[2] - X[0]))
    dy2   = (Y[2] + (Y[1] * -1.)) * invdx
    dy3   = (Y[3] + (Y[1] * -1.)) * (1. / (X[3] - X[1]))
    ddy2  = (dy2 + (dy1 * -1)) * invdx
    ddy3  = (dy3 + (dy2 * -1)) * invdx
    dddy3 = (ddy3 + (ddy2 * -1)) * invdx
    u = (t - X[1])
    return (Y[1] + dy1*u + ddy2*u*u + dddy3*u*u*(u-dx));

def fly():
    # Do a pseudocolor plot of u.
    v.DeleteAllPlots()
    v.OpenDatabase(gw_database)
    v.AddPlot("Pseudocolor", "Strain")
    pca = v.PseudocolorAttributes()
    pca.scaling = pca.Linear  # Linear, Log, Skew
    #pca.skewFactor = 1
    pca.limitsMode = pca.OriginalData  # OriginalData, ActualData
    #pca.minFlag = 1
    pca.min = -10
    pca.maxFlag = 1
    pca.max = 10
    pca.centering = pca.Nodal  # Natural, Nodal, Zonal
    pca.colorTableName = "viridis_light" # imola-seq and viridis_light are nice
    pca.invertColorTable = 0
    pca.opacityType = pca.FullyOpaque  
    pca.pointType = pca.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere
    pca.pointSizeVarEnabled = 0
    pca.pointSizeVar = "default"
    pca.pointSizePixels = 2
    pca.lineType = pca.Line  # Line, Tube, Ribbon
    pca.lineWidth = 0
    pca.renderSurfaces = 1
    pca.renderWireframe = 0
    pca.renderPoints = 0
    pca.smoothingLevel = 0
    pca.legendFlag = 1
    pca.lightingFlag = 1
    v.SetPlotOptions(pca)
    v.DrawPlots()
    v.ResetView()

    a = v.AnnotationAttributes()
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
    v.SetAnnotationAttributes(a)
 
    # Create the control points for the views.
    '''
    c0 = v.View3DAttributes()
    c0.viewNormal = (0, 0, 1)
    c0.focus = (0.5, 0.5, 0)
    c0.viewUp = (0, 1, 0)
    c0.viewAngle = 30
    c0.parallelScale = 140.719
    c0.nearPlane = -281.438
    c0.farPlane = 281.438
    c0.imagePan = (0, 0)
    c0.imageZoom = 1
    c0.perspective = 1
    c0.eyeAngle = 2
    c0.centerOfRotationSet = 0
    c0.centerOfRotation = (0.5, 0.5, 0)
    c0.axis3DScaleFlag = 0
    c0.axis3DScales = (1, 1, 1)
    c0.shear = (0, 0, 1)
    c0.windowValid = 1

    c1 = v.View3DAttributes()
    c1.viewNormal = (0.555427, -0.530184, 0.640629)
    c1.focus = (0.5, 0.5, 0)
    c1.viewUp = (-0.423414, 0.482733, 0.766609)
    c1.viewAngle = 30
    c1.parallelScale = 140.719
    c1.nearPlane = -281.438
    c1.farPlane = 281.438
    c1.imagePan = (0, 0)
    c1.imageZoom = 1
    c1.perspective = 1
    c1.eyeAngle = 2
    c1.centerOfRotationSet = 0
    c1.centerOfRotation = (0.5, 0.5, 0)
    c1.axis3DScaleFlag = 0
    c1.axis3DScales = (1, 1, 1)
    c1.shear = (0, 0, 1)
    c1.windowValid = 1

    c2 = v.View3DAttributes()
    c2.viewNormal = (0.841859, -0.0643092, 0.535852)
    c2.focus = (0.5, 0.5, 0)
    c2.viewUp = (-0.536986, -0.000412336, 0.843591)
    c2.viewAngle = 30
    c2.parallelScale = 140.719
    c2.nearPlane = -281.438
    c2.farPlane = 281.438
    c2.imagePan = (0, 0)
    c2.imageZoom = 1.4641
    c2.perspective = 1
    c2.eyeAngle = 2
    c2.centerOfRotationSet = 0
    c2.centerOfRotation = (0.5, 0.5, 0)
    c2.axis3DScaleFlag = 0
    c2.axis3DScales = (1, 1, 1)
    c2.shear = (0, 0, 1)
    c2.windowValid = 1

    c3 = v.View3DAttributes()
    c3.viewNormal = (0.583545, 0.609623, 0.536503)
    c3.focus = (0.5, 0.5, 0)
    c3.viewUp = (-0.343279, -0.413548, 0.84329)
    c3.viewAngle = 30
    c3.parallelScale = 140.719
    c3.nearPlane = -281.438
    c3.farPlane = 281.438
    c3.imagePan = (0, 0)
    c3.imageZoom = 2.14359
    c3.perspective = 1
    c3.eyeAngle = 2
    c3.centerOfRotationSet = 0
    c3.centerOfRotation = (0.5, 0.5, 0)
    c3.axis3DScaleFlag = 0
    c3.axis3DScales = (1, 1, 1)
    c3.shear = (0, 0, 1)
    c3.windowValid = 1


    c4 = v.View3DAttributes()
    c4.viewNormal = (-0.0704091, 0.867701, 0.492075)
    c4.focus = (0.5, 0.5, 0)
    c4.viewUp = (0.00193911, -0.493179, 0.869926)
    c4.viewAngle = 30
    c4.parallelScale = 140.719
    c4.nearPlane = -281.438
    c4.farPlane = 281.438
    c4.imagePan = (0, 0)
    c4.imageZoom = 2.59374
    c4.perspective = 1
    c4.eyeAngle = 2
    c4.centerOfRotationSet = 0
    c4.centerOfRotation = (0.5, 0.5, 0)
    c4.axis3DScaleFlag = 0
    c4.axis3DScales = (1, 1, 1)
    c4.shear = (0, 0, 1)
    c4.windowValid = 1


    c5 = v.View3DAttributes()
    c5.viewNormal = (-0.543602, 0.661879, 0.516151)
    c5.focus = (0.5, 0.5, 0)
    c5.viewUp = (0.273152, -0.441969, 0.854431)
    c5.viewAngle = 30
    c5.parallelScale = 140.719
    c5.nearPlane = -281.438
    c5.farPlane = 281.438
    c5.imagePan = (0, 0)
    c5.imageZoom = 3.7975
    c5.perspective = 1
    c5.eyeAngle = 2
    c5.centerOfRotationSet = 0
    c5.centerOfRotation = (0.5, 0.5, 0)
    c5.axis3DScaleFlag = 0
    c5.axis3DScales = (1, 1, 1)
    c5.shear = (0, 0, 1)
    c5.windowValid = 1

    c6 = v.View3DAttributes()
    c6.viewNormal = (-0.603814, 0.324425, 0.728119)
    c6.focus = (0.5, 0.5, 0)
    c6.viewUp = (0.72038, -0.168969, 0.672683)
    c6.viewAngle = 30
    c6.parallelScale = 140.719
    c6.nearPlane = -281.438
    c6.farPlane = 281.438
    c6.imagePan = (0, 0)
    c6.imageZoom = 2.59374
    c6.perspective = 1
    c6.eyeAngle = 2
    c6.centerOfRotationSet = 0
    c6.centerOfRotation = (0.5, 0.5, 0)
    c6.axis3DScaleFlag = 0
    c6.axis3DScales = (1, 1, 1)
    c6.shear = (0, 0, 1)
    c6.windowValid = 1
    '''

    c0 = v.View3DAttributes()
    c0.viewNormal = (0, 0, 1)
    c0.focus = (0, 0, -4.12473)
    c0.viewUp = (0, 1, 0)
    c0.viewAngle = 30
    c0.parallelScale = 282.904
    c0.nearPlane = -565.807
    c0.farPlane = 565.807
    c0.imagePan = (0, 0)
    c0.imageZoom = 2.59374
    c0.perspective = 1
    c0.eyeAngle = 2
    c0.centerOfRotationSet = 0
    c0.centerOfRotation = (0, 0, -4.12473)
    c0.axis3DScaleFlag = 0
    c0.axis3DScales = (1, 1, 1)
    c0.shear = (0, 0, 1)
    c0.windowValid = 1

    c1 = v.View3DAttributes()
    c1.viewNormal = (0.480004, -0.49164, 0.726558)
    c1.focus = (0, 0, -4.12473)
    c1.viewUp = (-0.532772, 0.494613, 0.686668)
    c1.viewAngle = 30
    c1.parallelScale = 282.904
    c1.nearPlane = -565.807
    c1.farPlane = 565.807
    c1.imagePan = (0, 0)
    c1.imageZoom = 3.13843
    c1.perspective = 1
    c1.eyeAngle = 2
    c1.centerOfRotationSet = 0
    c1.centerOfRotation = (0, 0, -4.12473)
    c1.axis3DScaleFlag = 0
    c1.axis3DScales = (1, 1, 1)
    c1.shear = (0, 0, 1)
    c1.windowValid = 1

    c2 = v.View3DAttributes()
    c2.viewNormal = (0.792583, -0.118684, 0.598102)
    c2.focus = (0, 0, -4.12473)
    c2.viewUp = (-0.595443, 0.0606849, 0.801102)
    c2.viewAngle = 30
    c2.parallelScale = 282.904
    c2.nearPlane = -565.807
    c2.farPlane = 565.807
    c2.imagePan = (0, 0)
    c2.imageZoom = 3.7975
    c2.perspective = 1
    c2.eyeAngle = 2
    c2.centerOfRotationSet = 0
    c2.centerOfRotation = (0, 0, -4.12473)
    c2.axis3DScaleFlag = 0
    c2.axis3DScales = (1, 1, 1)
    c2.shear = (0, 0, 1)
    c2.windowValid = 1

    c3 = v.View3DAttributes()
    c3.viewNormal = (0.82836, 0.17747, 0.531342)
    c3.focus = (0, 0, -4.12473)
    c3.viewUp = (-0.513939, -0.136656, 0.846872)
    c3.viewAngle = 30
    c3.parallelScale = 282.904
    c3.nearPlane = -565.807
    c3.farPlane = 565.807
    c3.imagePan = (0, 0)
    c3.imageZoom = 4.59497
    c3.perspective = 1
    c3.eyeAngle = 2
    c3.centerOfRotationSet = 0
    c3.centerOfRotation = (0, 0, -4.12473)
    c3.axis3DScaleFlag = 0
    c3.axis3DScales = (1, 1, 1)
    c3.shear = (0, 0, 1)
    c3.windowValid = 1

    c4 = v.View3DAttributes()
    c4.viewNormal = (0.474281, 0.56115, 0.678358)
    c4.focus = (0, 0, -4.12473)
    c4.viewUp = (-0.410203, -0.540923, 0.734259)
    c4.viewAngle = 30
    c4.parallelScale = 282.904
    c4.nearPlane = -565.807
    c4.farPlane = 565.807
    c4.imagePan = (0, 0)
    c4.imageZoom = 3.7975
    c4.perspective = 1
    c4.eyeAngle = 2
    c4.centerOfRotationSet = 0
    c4.centerOfRotation = (0, 0, -4.12473)
    c4.axis3DScaleFlag = 0
    c4.axis3DScales = (1, 1, 1)
    c4.shear = (0, 0, 1)
    c4.windowValid = 1

    c5 = v.View3DAttributes()
    c5.viewNormal = (0.474281, 0.56115, 0.678358)
    c5.focus = (0, 0, -4.12473)
    c5.viewUp = (-0.410203, -0.540923, 0.734259)
    c5.viewAngle = 30
    c5.parallelScale = 282.904
    c5.nearPlane = -565.807
    c5.farPlane = 565.807
    c5.imagePan = (0, 0)
    c5.imageZoom = 3.7975
    c5.perspective = 1
    c5.eyeAngle = 2
    c5.centerOfRotationSet = 0
    c5.centerOfRotation = (0, 0, -4.12473)
    c5.axis3DScaleFlag = 0
    c5.axis3DScales = (1, 1, 1)
    c5.shear = (0, 0, 1)
    c5.windowValid = 1

    # Create a tuple of camera values and x values. The x values are weights
    # that help to determine where in [0,1] the control points occur.
    cpts = (c0, c1, c2, c3, c4, c5)
    x=[]
    for i in range(7): # Personal Note: change this if you have more than 7 control points
        x = x + [float(i) / float(6.)]
 
    # Animate the camera. Note that we use the new built-in EvalCubicSpline
    # function which takes a t value from [0,1] a tuple of t values and a tuple
    # of control points. In this case, the control points are View3DAttributes
    # objects that we are using to animate the camera but they can be any object
    # that supports +, * operators.
    nsteps = 650
    for i in range(nsteps):
        t = float(i) / float(nsteps - 1)
        c = EvalCubicSpline(t, x, cpts)
        #c.nearPlane = -34.461
        #c.farPlane = 34.461
        v.SetView3D(c)
        AniAtts = v.AnimationAttributes()
        AniAtts.frameIncrement = 10
        v.SetAnimationAttributes(AniAtts)
        v.TimeSliderNextState()
        s = v.SaveWindowAttributes()
        s.fileName = "gw_mesh_camera_fly_test"
        s.format = s.PNG
        s.progressive = 1
        s.width = 772
        s.height = 702
        s.screenCapture = 1
        v.SetSaveWindowAttributes(s)
        # Save the window
        v.SaveWindow()
        # For moviemaking...
        # SaveWindow()
 
fly()