import vtk
import math
import numpy as np

output_directory = "/home/guest/Documents/BH_Vis_local/data/mesh/gw_test6_polar/"

# Define the dimensions of the grid
numTheta = 360
numRadius = 200
display_radius = 100

# Create vtkUnstructuredGrid
grid = vtk.vtkUnstructuredGrid()

# Define the points and their coordinates
points = vtk.vtkPoints()
index = 0
for j in range(numRadius):
    for i in range(numTheta):
        theta = 2.0 * vtk.vtkMath.Pi() * float(i) / float(numTheta)
        radius = display_radius * float(j) / float(numRadius - 1) # Scaling radius from 0 to 1
        x = radius * math.cos(theta)
        y = radius * math.sin(theta)
        z = np.sin(radius * 10) + np.cos(radius * 10) # Simple sine function
        points.InsertNextPoint(x, y, z)
        index += 1

# Set the points in the vtkUnstructuredGrid
grid.SetPoints(points)

# Define the connectivity of the cells
cellArray = vtk.vtkCellArray()
for j in range(numRadius - 1):
    for i in range(numTheta):
        cell = vtk.vtkQuad()
        cell.GetPointIds().SetId(0, i + j * numTheta)
        cell.GetPointIds().SetId(1, (i + 1) % numTheta + j * numTheta)
        cell.GetPointIds().SetId(2, (i + 1) % numTheta + (j + 1) * numTheta)
        cell.GetPointIds().SetId(3, i + (j + 1) * numTheta)
        cellArray.InsertNextCell(cell)


# Set the cells in the vtkUnstructuredGrid
grid.SetCells(vtk.VTK_QUAD, cellArray)

# Write to a .vts file
writer = vtk.vtkXMLUnstructuredGridWriter()
filename = output_directory + f"/state.vtu"
writer.SetFileName(filename)
writer.SetInputData(grid)
writer.Write()
