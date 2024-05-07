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
#include <random>

#include "vtk-util.hpp"
#include "../util/util.hpp"


void PointsWithPartitionsToVtk(std::vector<double>& point_coords, std::vector<uint64_t>& partitions, uint64_t count, std::string out_file_name){
    assert(point_coords.size() == (count*3));
    assert(partitions.size() == count);



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
        // std::srand(pow(partition+2,3));
        std::mt19937 gen(pow(partition+2,3));
        std::uniform_int_distribution<unsigned int> dis(0, 255);
        vcolors->SetValue(i, static_cast<unsigned char>(dis(gen)));
        vcolors->SetValue(i+1, static_cast<unsigned char>(dis(gen)));
        vcolors->SetValue(i+2, static_cast<unsigned char>(dis(gen)));


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

void ElementsWithPartitionsToVtk(std::vector<ElementWithCoord>& elements, std::vector<uint16_t>& partitions, uint64_t count, std::string out_file_name){
    std::vector<double> coords(count*3);

    //TODO: can be parallelized
    for (size_t elem_i = 0; elem_i < count; elem_i++)
    {
        coords[3*elem_i] = elements[elem_i].x;
        coords[3*elem_i + 1] = elements[elem_i].y;
        coords[3*elem_i + 2] = elements[elem_i].z;
        
    }

    std::vector<uint64_t> partitions_(partitions.begin(), partitions.end());
    PointsWithPartitionsToVtk(coords,partitions_,count,out_file_name);
    
}



void PointsToVtk(std::vector<double>& point_coords, uint64_t count, std::string out_file_name){
    std::vector<uint64_t> dummy_scalars(count);
    std::fill(dummy_scalars.begin(), dummy_scalars.end(), 1);
    return PointsWithPartitionsToVtk(point_coords, dummy_scalars, count,  out_file_name);
}


