import vtk


def export_points_to_vtk(points,colors,file_name):
    vpoints = vtk.vtkPoints()
    vpoints.SetNumberOfPoints(len(points))
    for i in range(len(points)):
        vpoints.SetPoint(i, points[i])
    vpoly = vtk.vtkPolyData()
    vpoly.SetPoints(vpoints)

    if not colors is None:
        vcolors = vtk.vtkUnsignedCharArray()
        vcolors.SetNumberOfComponents(3)
        vcolors.SetName("Colors")
        vcolors.SetNumberOfTuples(len(points))
        for i in range(len(points)):
            vcolors.SetTuple(i ,[colors[i][0],colors[i][1], colors[i][2]])
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

    vcolors = vtk.vtkUnsignedIntArray()
    vcolors.SetNumberOfComponents(1)
    vcolors.SetName("pagerank")
    vcolors.SetNumberOfTuples(len(points))
    for i in range(len(points)):
        vcolors.SetTuple(i ,[scalar_values[i]])
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