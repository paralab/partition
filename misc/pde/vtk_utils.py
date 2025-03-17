import vtk
import random


def part_to_color(ci):
    random.seed(ci**3)
    return (random.randint(0,255),random.randint(0,255),random.randint(0,255),125)

def export_points_with_parts_to_vtk(points, parts, file_name):
    vpoints = vtk.vtkPoints()
    vpoints.SetNumberOfPoints(len(points))
    for i in range(len(points)):
        vpoints.SetPoint(i, points[i])
    vpoly = vtk.vtkPolyData()
    vpoly.SetPoints(vpoints)

    part_colors = [part_to_color(p) for p in parts]

    vcolors = vtk.vtkUnsignedCharArray()
    vcolors.SetNumberOfComponents(3)
    vcolors.SetName("Colors")
    vcolors.SetNumberOfTuples(len(points))
    for i in range(len(points)):
        vcolors.SetTuple(i ,[part_colors[i][0],part_colors[i][1], part_colors[i][2]])
    vpoly.GetPointData().SetScalars(vcolors)


    vcells = vtk.vtkCellArray()
    
    for i in range(len(points)):
        vcells.InsertNextCell(1)
        vcells.InsertCellPoint(i)
        
    vpoly.SetVerts(vcells)

    writer = vtk.vtkPolyDataWriter()

    writer.SetFileName(file_name)
    writer.SetInputData(vpoly)
    writer.Write()


def export_points_scalar_to_vtk(points, scalar_values, file_name):
    vpoints = vtk.vtkPoints()
    vpoints.SetNumberOfPoints(len(points))
    for i in range(len(points)):
        vpoints.SetPoint(i, points[i])
    vpoly = vtk.vtkPolyData()
    vpoly.SetPoints(vpoints)

    scalars = vtk.vtkFloatArray()
    scalars.SetNumberOfComponents(1)
    scalars.SetName("diffusion_value")
    scalars.SetNumberOfTuples(len(points))
    for i in range(len(points)):
        scalars.SetTuple(i ,[scalar_values[i]])
    vpoly.GetPointData().SetScalars(scalars)


    vcells = vtk.vtkCellArray()
    
    for i in range(len(points)):
        vcells.InsertNextCell(1)
        vcells.InsertCellPoint(i)
        
    vpoly.SetVerts(vcells)

    writer = vtk.vtkPolyDataWriter()

    writer.SetFileName(file_name)
    writer.SetInputData(vpoly)
    writer.Write()