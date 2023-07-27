# This is a refactored version of the original main script. 
# Select configuration options are available in 'config-animation-bh_merger-2D.yml'

import sys, csv, os, yaml
sys.path.append("/usr/local/3.3.1/linux-x86_64/lib/site-packages") 
import visit
visit.Launch()
import visit
v = visit

# Load the YAML config file
with open('/home/guest/Documents/bh_vis/scripts/config-animation-bh_merger-2D.yml', 'r') as file:
    config = yaml.safe_load(file)

# Extract directories and parameters from the config
directories = config.get('directories', {})
parameters = config.get('parameters', {})

# Apply the directories and parameters to your variables
parent_directory = directories.get('parent_directory', '')
bh_database = directories.get('bh_database', '')
bh_data = directories.get('bh_data', '')
gw_database = directories.get('gw_database', '')

record_render = parameters.get('record_render', False)
green_screen = parameters.get('green_screen', False)
show_axis_annotations = parameters.get('show_axis_annotations', False)

def apply_attributes(cls, attributes):
    obj = cls()
    for key, value in attributes.items():
        print(f"Applying attribute: {key} with value: {value}")
        keys = key.split('.')
        if len(keys) > 1:  # this is a nested attribute
            try:
                nested_obj = obj
                for nested_key in keys[:-1]:  # traverse to the deepest object
                    nested_obj = getattr(nested_obj, nested_key)
                setattr(nested_obj, keys[-1], value)  # set the nested attribute
            except Exception as e:
                print(f"Error setting nested attribute {key} with value {value}: {e}")
        else:  # this is a top-level attribute
            if isinstance(value, list) and 'Color' in key:
                try:
                    setattr(obj, key, tuple(value))  # Convert list to tuple
                except Exception as e:
                    print(f"Error setting color attribute {key} with value {tuple(value)}: {e}")
            elif isinstance(value, dict) and key == "special":
                for special_key, special_value in value.items():
                    try:
                        special_value = getattr(obj, special_value) # Fetch property using getattr
                        setattr(obj, special_key, special_value)
                    except Exception as e:
                        print(f"Error setting special attribute {special_key} with value {special_value}: {e}")
            else:
                try:
                    setattr(obj, key, value)
                except Exception as e:
                    print(f"Error setting attribute {key} with value {value}: {e}")
    return obj

def set_visit_attributes():
    global AnnotationAttributes, RenderingAtts, View3DAtts, BHPseudocolorAtts, GWPseudocolorAtts
    AnnotationAttributes = apply_attributes(v.AnnotationAttributes, config['AnnotationAttributes'])
    RenderingAtts = apply_attributes(v.RenderingAttributes, config['RenderingAtts'])
    View3DAtts = apply_attributes(v.View3DAttributes, config['View3DAtts'])
    BHPseudocolorAtts = apply_attributes(v.PseudocolorAttributes, config['BHPseudocolorAtts'])
    GWPseudocolorAtts = apply_attributes(v.PseudocolorAttributes, config['GWPseudocolorAtts'])

# Use the attributes and parameters
def apply_vis_configuration():
    v.SetRenderingAttributes(RenderingAtts)
    v.SetView3D(View3DAtts)
    v.SetAnnotationAttributes(AnnotationAttributes)
    if record_render:
        v.SetSaveWindowAttributes(v.SaveWindowAttributes())

set_visit_attributes()
apply_vis_configuration()

# Create black holes
def create_spheres():
    v.OpenDatabase(bh_database)
    v.AddPlot("Pseudocolor", "mesh_quality/warpage", 1, 1)
    v.AddOperator("Transform")
    v.SetPlotOptions(BHPseudocolorAtts)
    v.AddPlot("Pseudocolor", "mesh_quality/warpage", 1, 1)
    v.AddOperator("Transform")
    v.SetPlotOptions(BHPseudocolorAtts)

# Create gravitation wave mesh
def create_gw():
    v.OpenDatabase(gw_database)
    v.AddPlot("Pseudocolor", "Strain", 1, 1)
    v.SetPlotOptions(GWPseudocolorAtts)

create_spheres()
create_gw()

# Provide trajectory function for bh
def set_coords(objNum, x, y, z):
    v.SetActivePlots(objNum)
    trasnformAtts = v.TransformAttributes()
    trasnformAtts.doTranslate = 1
    trasnformAtts.translateX = x
    trasnformAtts.translateY = y
    trasnformAtts.translateZ = z
    v.SetOperatorOptions(trasnformAtts, 0, 0)

with open(bh_data, 'r') as file:
    reader = csv.reader(file)
    header = next(reader)
    v.DrawPlots()
    for row in reader:
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