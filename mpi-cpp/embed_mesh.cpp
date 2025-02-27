#include <gmsh.h>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <limits>
#include <cmath>
#include <cassert>
#include <omp.h>
#include <cstdlib>

#define STRUCTURED_VOLUME_TAG   999

// Function to check if a point is inside any 3D element of the shape
bool is_point_inside_shape(double x, double y, double z, 
                          const std::vector<double>& shape_elem_bounding_boxes,
                          size_t shape_elem_count,
                          double shape_bounding_min_x, double shape_bounding_min_y, double shape_bounding_min_z,
                          double shape_bounding_max_x, double shape_bounding_max_y, double shape_bounding_max_z) {
    // checking whether the point is completely outside the global bbox
    if (x < shape_bounding_min_x || x > shape_bounding_max_x ||
        y < shape_bounding_min_y || y > shape_bounding_max_y ||
        z < shape_bounding_min_z || z > shape_bounding_max_z) {
        return false;
    }
    
    for (size_t i = 0; i < shape_elem_count; i++) {
        auto xmin = shape_elem_bounding_boxes[6*i + 0];
        auto xmax = shape_elem_bounding_boxes[6*i + 1];
        auto ymin = shape_elem_bounding_boxes[6*i + 2];
        auto ymax = shape_elem_bounding_boxes[6*i + 3];
        auto zmin = shape_elem_bounding_boxes[6*i + 4];
        auto zmax = shape_elem_bounding_boxes[6*i + 5];

        if (xmin <= x && x <= xmax &&
            ymin <= y && y <= ymax &&
            zmin <= z && z <= zmax) {
                // std::cout << "inside\n";
            return true;
        }
    }
    return false;
}

int main(int argc, char **argv) {
    if (argc != 4)
    {
        std::cout << "Usage:\n<exec> given_mesh_file output_mesh_file max_divisions_per_side\n";
        return 1;
    }
    

    std::string shape_mesh_file = argv[1];
    std::string output_file = argv[2];
    int max_divisions_per_dim = atoi(argv[3]);

    // Initialize Gmsh
    gmsh::initialize();
    gmsh::model::add("StructuredHexWithEmbeddedShape");

    // Import the unstructured mesh (shape)
    gmsh::merge(shape_mesh_file);
    gmsh::model::occ::synchronize();

    // Get bounding box of imported shape
    std::vector<std::pair<int, int>> volumes;
    gmsh::model::getEntities(volumes, 3);
    
    if (volumes.empty()) {
        std::cout << "No volume found in the imported mesh." << std::endl;
        gmsh::finalize();
        return 1;
    }

    double shape_bounding_min_x, shape_bounding_min_y, shape_bounding_min_z;
    double shape_bounding_max_x, shape_bounding_max_y, shape_bounding_max_z;
    
    gmsh::model::getBoundingBox(3, volumes[0].second, 
                               shape_bounding_min_x, shape_bounding_min_y, shape_bounding_min_z,
                               shape_bounding_max_x, shape_bounding_max_y, shape_bounding_max_z);

    // Parameters for structured mesh
    double box_size_x = shape_bounding_max_x - shape_bounding_min_x;
    double box_size_y = shape_bounding_max_y - shape_bounding_min_y;
    double box_size_z = shape_bounding_max_z - shape_bounding_min_z;

    // expanding the bounding box by 1% in each direction
    shape_bounding_min_x -= box_size_x * 0.01;
    shape_bounding_min_y -= box_size_y * 0.01;
    shape_bounding_min_z -= box_size_z * 0.01;

    shape_bounding_max_x += box_size_x * 0.01;
    shape_bounding_max_y += box_size_y * 0.01;
    shape_bounding_max_z += box_size_z * 0.01;

    box_size_x = shape_bounding_max_x - shape_bounding_min_x;
    box_size_y = shape_bounding_max_y - shape_bounding_min_y;
    box_size_z = shape_bounding_max_z - shape_bounding_min_z;

    double box_size_xyz_max = std::max(box_size_x, std::max(box_size_y, box_size_z));
    double elem_size = box_size_xyz_max / max_divisions_per_dim;

    // Now we create the structured hex mesh
    int p1 = gmsh::model::geo::addPoint(shape_bounding_min_x, shape_bounding_min_y, shape_bounding_min_z);
    int p2 = gmsh::model::geo::addPoint(shape_bounding_max_x, shape_bounding_min_y, shape_bounding_min_z);
    int p3 = gmsh::model::geo::addPoint(shape_bounding_max_x, shape_bounding_max_y, shape_bounding_min_z);
    int p4 = gmsh::model::geo::addPoint(shape_bounding_min_x, shape_bounding_max_y, shape_bounding_min_z);
    int p5 = gmsh::model::geo::addPoint(shape_bounding_min_x, shape_bounding_min_y, shape_bounding_max_z);
    int p6 = gmsh::model::geo::addPoint(shape_bounding_max_x, shape_bounding_min_y, shape_bounding_max_z);
    int p7 = gmsh::model::geo::addPoint(shape_bounding_max_x, shape_bounding_max_y, shape_bounding_max_z);
    int p8 = gmsh::model::geo::addPoint(shape_bounding_min_x, shape_bounding_max_y, shape_bounding_max_z);

    // Create lines
    // Bottom face
    int l1 = gmsh::model::geo::addLine(p1, p2);
    int l2 = gmsh::model::geo::addLine(p2, p3);
    int l3 = gmsh::model::geo::addLine(p3, p4);
    int l4 = gmsh::model::geo::addLine(p4, p1);

    // Top face
    int l5 = gmsh::model::geo::addLine(p5, p6);
    int l6 = gmsh::model::geo::addLine(p6, p7);
    int l7 = gmsh::model::geo::addLine(p7, p8);
    int l8 = gmsh::model::geo::addLine(p8, p5);

    // Vertical edges
    int l9 = gmsh::model::geo::addLine(p1, p5);
    int l10 = gmsh::model::geo::addLine(p2, p6);
    int l11 = gmsh::model::geo::addLine(p3, p7);
    int l12 = gmsh::model::geo::addLine(p4, p8);

    // Create surfaces
    // Bottom face
    std::vector<int> cl1 = {l1, l2, l3, l4};
    int s1 = gmsh::model::geo::addCurveLoop(cl1);
    int surf1 = gmsh::model::geo::addPlaneSurface({s1});

    // Top face
    std::vector<int> cl2 = {l5, l6, l7, l8};
    int s2 = gmsh::model::geo::addCurveLoop(cl2);
    int surf2 = gmsh::model::geo::addPlaneSurface({s2});

    // Side faces
    std::vector<int> cl3 = {l1, l10, -l5, -l9};
    int s3 = gmsh::model::geo::addCurveLoop(cl3);
    int surf3 = gmsh::model::geo::addPlaneSurface({s3});

    std::vector<int> cl4 = {l2, l11, -l6, -l10};
    int s4 = gmsh::model::geo::addCurveLoop(cl4);
    int surf4 = gmsh::model::geo::addPlaneSurface({s4});

    std::vector<int> cl5 = {l3, l12, -l7, -l11};
    int s5 = gmsh::model::geo::addCurveLoop(cl5);
    int surf5 = gmsh::model::geo::addPlaneSurface({s5});

    std::vector<int> cl6 = {l4, l9, -l8, -l12};
    int s6 = gmsh::model::geo::addCurveLoop(cl6);
    int surf6 = gmsh::model::geo::addPlaneSurface({s6});

    // Create volume
    std::vector<int> sl = {surf1, surf2, surf3, surf4, surf5, surf6};
    int surfLoop = gmsh::model::geo::addSurfaceLoop(sl);
    int vol = gmsh::model::geo::addVolume({surfLoop}, STRUCTURED_VOLUME_TAG);

    // Synchronize the model
    gmsh::model::geo::synchronize();

    // Set mesh size at each point
    std::vector<std::pair<int, int>> points;
    gmsh::model::getEntities(points, 0);
    // std::vector<std::pair<int, int>> pointTags;
    // for (auto& p : points) {
    //     pointTags.push_back(p);
    // }
    gmsh::model::mesh::setSize(points, elem_size);

    gmsh::model::geo::synchronize();

    // Set transfinite mesh to create structured grid
    // Calculate number of nodes in each direction
    int nx = std::max(2, static_cast<int>(box_size_x / elem_size) + 1);
    int ny = std::max(2, static_cast<int>(box_size_y / elem_size) + 1);
    int nz = std::max(2, static_cast<int>(box_size_z / elem_size) + 1);

    // Set transfinite curves
    gmsh::model::mesh::setTransfiniteCurve(l1, nx);
    gmsh::model::mesh::setTransfiniteCurve(l2, ny);
    gmsh::model::mesh::setTransfiniteCurve(l3, nx);
    gmsh::model::mesh::setTransfiniteCurve(l4, ny);

    gmsh::model::mesh::setTransfiniteCurve(l5, nx);
    gmsh::model::mesh::setTransfiniteCurve(l6, ny);
    gmsh::model::mesh::setTransfiniteCurve(l7, nx);
    gmsh::model::mesh::setTransfiniteCurve(l8, ny);

    gmsh::model::mesh::setTransfiniteCurve(l9, nz);
    gmsh::model::mesh::setTransfiniteCurve(l10, nz);
    gmsh::model::mesh::setTransfiniteCurve(l11, nz);
    gmsh::model::mesh::setTransfiniteCurve(l12, nz);

    // Set transfinite surfaces with recombination to get quadrangular elements
    std::vector<int> surfaces = {surf1, surf2, surf3, surf4, surf5, surf6};
    for (int s : surfaces) {
        gmsh::model::mesh::setTransfiniteSurface(s);
        gmsh::model::mesh::setRecombine(2, s);  // This ensures quad faces
    }

    // Set transfinite volume
    gmsh::model::mesh::setTransfiniteVolume(vol);

    // Force hex mesh generation
    gmsh::option::setNumber("Mesh.Algorithm3D", 1);  // Use Delaunay for 3D mesh
    gmsh::option::setNumber("Mesh.RecombinationAlgorithm", 0);  // Standard recombination algorithm
    gmsh::option::setNumber("Mesh.Recombine3DAll", 1);  // Recombine all tetrahedra into hexahedra if possible
    gmsh::option::setNumber("Mesh.SubdivisionAlgorithm", 1);  // All quads

    // Generate 3D mesh
    gmsh::model::mesh::generate(3);

    gmsh::model::geo::synchronize();

    // Get structured mesh elements
    std::vector<int> elementTypes;
    std::vector<std::vector<size_t>> elementTags;
    std::vector<std::vector<size_t>> elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags, 3, STRUCTURED_VOLUME_TAG);
    
    std::vector<size_t> structured_elemTags = elementTags[0];
    std::vector<size_t> structured_elemNodeTags = elementNodeTags[0];
    
    std::vector<double> structured_elems_barycenters;
    gmsh::model::mesh::getBarycenters(5, STRUCTURED_VOLUME_TAG, false, true, structured_elems_barycenters);

    // Get unstructured shape elements
    std::vector<std::pair<int, int>> shape_volumes;
    for (auto& vol : volumes) {
        if (vol.second != STRUCTURED_VOLUME_TAG) {
            shape_volumes.push_back(vol);
        }
    }
    
    std::vector<std::vector<int>> shape_elemTypes_grouped;
    std::vector<std::vector<std::vector<size_t>>> shape_elemTags_grouped;
    std::vector<std::vector<std::vector<size_t>>> shape_elemNodeTags_grouped;
    
    for (auto& vol : shape_volumes) {
        std::vector<int> types;
        std::vector<std::vector<size_t>> tags;
        std::vector<std::vector<size_t>> nodeTags;
        gmsh::model::mesh::getElements(types, tags, nodeTags, vol.first, vol.second);
        shape_elemTypes_grouped.push_back(types);
        shape_elemTags_grouped.push_back(tags);
        shape_elemNodeTags_grouped.push_back(nodeTags);
    }

    std::vector<size_t> shape_elemTags;
    std::vector<size_t> shape_elemNodeTags;
    
    for (size_t i = 0; i < shape_elemTypes_grouped.size(); i++) {
        for (size_t j = 0; j < shape_elemTypes_grouped[i].size(); j++) {
            if (shape_elemTypes_grouped[i][j] == 5) { // for now, getting only hex elements
                shape_elemTags.insert(shape_elemTags.end(), shape_elemTags_grouped[i][j].begin(), shape_elemTags_grouped[i][j].end());
                shape_elemNodeTags.insert(shape_elemNodeTags.end(), shape_elemNodeTags_grouped[i][j].begin(), shape_elemNodeTags_grouped[i][j].end());
            }
        }
    }

    std::vector<double> shape_elem_bounding_boxes(shape_elemTags.size() * 6);

    for (size_t i = 0; i < shape_elemTags.size(); i++) {
        double* bbox = &shape_elem_bounding_boxes[6*i];
        bbox[0] = std::numeric_limits<double>::infinity();
        bbox[1] = -std::numeric_limits<double>::infinity();
        bbox[2] = std::numeric_limits<double>::infinity();
        bbox[3] = -std::numeric_limits<double>::infinity();
        bbox[4] = std::numeric_limits<double>::infinity();
        bbox[5] = -std::numeric_limits<double>::infinity();
        
        for (size_t j = 0; j < 8; j++) {  // 8 nodes per hex element
            size_t node = shape_elemNodeTags[8*i + j];
            std::vector<double> coord;
            std::vector<double> paramCoord;
            int entity;
            int dim;
            gmsh::model::mesh::getNode(node, coord, paramCoord, entity, dim);
            
            bbox[0] = std::min(coord[0], bbox[0]);
            bbox[1] = std::max(coord[0], bbox[1]);
            bbox[2] = std::min(coord[1], bbox[2]);
            bbox[3] = std::max(coord[1], bbox[3]);
            bbox[4] = std::min(coord[2], bbox[4]);
            bbox[5] = std::max(coord[2], bbox[5]);
        }
    }

    std::cout << "checking structured elements to remove/add" << std::endl;

    std::vector<bool> elemtags_to_remove_flag(structured_elemTags.size(), false);
    // Remove structured elements outside the shape

    #pragma omp parallel for
    for (size_t i = 0; i < structured_elemTags.size(); i++) {
        // std::cout << 100.0 * (i+1) / structured_elemTags.size() << " %" << std::endl;
        double x = structured_elems_barycenters[3*i];
        double y = structured_elems_barycenters[3*i+1];
        double z = structured_elems_barycenters[3*i+2];
        
        if (!is_point_inside_shape(x, y, z, shape_elem_bounding_boxes, shape_elemTags.size(), 
                                   shape_bounding_min_x, shape_bounding_min_y, shape_bounding_min_z,
                                   shape_bounding_max_x, shape_bounding_max_y, shape_bounding_max_z)) 
        {
            // elemtags_to_remove.push_back(structured_elemTags[i]);
            elemtags_to_remove_flag[i] = true;
            
        }
    }

    std::vector<size_t> elemtags_to_remove;

    for (size_t i = 0; i < structured_elemTags.size(); i++)
    {
        if (elemtags_to_remove_flag[i])
        {
            elemtags_to_remove.push_back(structured_elemTags[i]);
        }
        
    }
    
    gmsh::model::mesh::setVisibility(elemtags_to_remove, 0);

    // Convert vector of size_t to vector of int for setVisibility
    // std::vector<int> shape_elemTags_int(shape_elemTags.size());
    // for (size_t i = 0; i < shape_elemTags.size(); i++) {
    //     shape_elemTags_int[i] = static_cast<int>(shape_elemTags[i]);
    // }
    gmsh::model::mesh::setVisibility(shape_elemTags, 0);

    gmsh::plugin::setNumber("Invisible", "DeleteElements", 1);
    gmsh::plugin::run("Invisible");

    gmsh::model::occ::synchronize();

    gmsh::write(output_file);

    // gmsh::fltk::run();
    gmsh::finalize();

    return 0;
}