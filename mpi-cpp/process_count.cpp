/**
 * Given a mesh file path and grain size, returns the process count (= mesh element count / grain size)
 * 
 */

#include <iostream>
#include <stdexcept>
#include "gmsh.h"

int main(int argc, char const *argv[])
{
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <mesh file path> <gain size>" << std::endl;
        return 1; // indicating an error
    }
    const std::string mesh_file_path = argv[1];
    int grain_size = static_cast<size_t>(std::stoi(argv[2]));

    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 0);
    gmsh::open(mesh_file_path);

    std::vector<std::pair<int, int>> dimTags;
    gmsh::model::getEntities(dimTags);

    std::vector<int> elementTypes;
    gmsh::model::mesh::getElementTypes(elementTypes, 3);

    

    if (elementTypes.size() == 0)
    {
        throw std::invalid_argument("no 3D elements\t exiting...");
    }

    if (elementTypes.size() > 1)
    {
        throw std::invalid_argument("more than 1 element type\t exiting...");
    }
    int gmsh_element_type = elementTypes[0];
    double element_count;

    switch (gmsh_element_type)
    {
    case 4: // linear tet
    {
        gmsh::option::getNumber("Mesh.NbTetrahedra", element_count);
        break;
    }
    case 5: // linear hexahedra
    {
        gmsh::option::getNumber("Mesh.NbHexahedra", element_count);

        break;
    }

    default:
    {
        throw std::invalid_argument("unknown element type\t exiting...");
        break;
    }
    }
    gmsh::finalize();

    std::cout << static_cast<size_t>(element_count)/grain_size;

    return 0;
}
