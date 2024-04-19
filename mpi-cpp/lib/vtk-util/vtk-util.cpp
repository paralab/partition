#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataWriter.h>

#include <cassert>

#include <string>
#include <algorithm>

void PointsWithScalarsToVtk(std::vector<double>& point_coords, std::vector<double>& scalars, uint64_t count, std::string out_file_name){
    assert(point_coords.size() == (count*3));
    assert(scalars.size() == count);

    // Create vtkPoints object
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    // Create vtkDoubleArray for scalar values
    vtkSmartPointer<vtkDoubleArray> vtk_scalars = vtkSmartPointer<vtkDoubleArray>::New();
    vtk_scalars->SetName("scalar");
    vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();

    // Add points to vtkPoints
    for (size_t i = 0; i < point_coords.size(); i += 3) {
        double x = point_coords[i];
        double y = point_coords[i + 1];
        double z = point_coords[i + 2];
        points->InsertNextPoint(x, y, z);
        vtk_scalars->InsertNextValue(scalars[i/3]);
        vtkIdType pointId = i / 3; // Calculate the point ID
        vertices->InsertNextCell(1, &pointId);
    }

    // Create vtkPolyData object and set points
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);

    // Add scalar data to point data
    polydata->GetPointData()->SetScalars(vtk_scalars);
    
    polydata->SetVerts(vertices);
    // Write vtkPolyData to VTK file
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    char* out_file_name_ = new char[out_file_name.length() + 1];
    strcpy(out_file_name_, out_file_name.c_str()); 
    writer->SetFileName(out_file_name_);
    writer->SetInputData(polydata);
    writer->Write();
}

void PointsToVtk(std::vector<double>& point_coords, uint64_t count, std::string out_file_name){
    std::vector<double> dummy_scalars(count);
    std::fill(dummy_scalars.begin(), dummy_scalars.end(), 1);
    return PointsWithScalarsToVtk(point_coords, dummy_scalars, count,  out_file_name);
}