#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataWriter.h>
#include <vtkTuple.h>


#include <cassert>

#include <string>
#include <algorithm>

#include <vtk-util.hpp>
#include <util.hpp>

void PointsWithPartitionsToVtk(std::vector<double>& point_coords, std::vector<uint64_t>& partitions, uint64_t count, std::string out_file_name){
    assert(point_coords.size() == (count*3));
    assert(partitions.size() == count);

    std::vector<std::vector<unsigned char>> point_colors(count);
    for (size_t i = 0; i < count; i++)
    {
        std::srand(pow(partitions[i]+1,3));

        point_colors[i] = {(unsigned char)(rand()%256),(unsigned char)(rand()%256),(unsigned char)(rand()%256)};
    }
    


    // Create vtkPoints object
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    // Create vtkDoubleArray for scalar values
    // vtkSmartPointer<vtkDoubleArray> vtk_scalars = vtkSmartPointer<vtkDoubleArray>::New();
    // vtk_scalars->SetName("scalar");
    vtkSmartPointer<vtkUnsignedCharArray> vcolors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    vcolors->SetNumberOfComponents(3);
    vcolors->SetName("Colors");
    vcolors->SetNumberOfTuples(count);


    vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();

    // Add points to vtkPoints
    for (size_t i = 0; i < point_coords.size(); i += 3) {
        double x = point_coords[i];
        double y = point_coords[i + 1];
        double z = point_coords[i + 2];
        points->InsertNextPoint(x, y, z);
        auto partition = partitions[i/3];
        std::srand(pow(partition+1,3));
        vcolors->SetValue(i, (unsigned char)(rand()%256));
        vcolors->SetValue(i+1, (unsigned char)(rand()%256));
        vcolors->SetValue(i+2, (unsigned char)(rand()%256));


        // vcolors->SetTuple(i/3,i/3, point_colors.data());
        vtkIdType pointId = i / 3; // Calculate the point ID
        vertices->InsertNextCell(1, &pointId);
    }

    // Create vtkPolyData object and set points
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);

    // Add scalar data to point data
    polydata->GetPointData()->SetScalars(vcolors);
    
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
    std::vector<uint64_t> dummy_scalars(count);
    std::fill(dummy_scalars.begin(), dummy_scalars.end(), 1);
    return PointsWithPartitionsToVtk(point_coords, dummy_scalars, count,  out_file_name);
}


