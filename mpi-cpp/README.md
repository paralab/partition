Cmake dependencies

- Gmsh SDK https://gmsh.info/#Download
- VTK https://vtk.org/download/
- METIS https://github.com/KarypisLab/METIS/tags
- ParMETIS (depends on METIS) https://github.com/KarypisLab/ParMETIS
- GKLIB (required for METIS) https://github.com/KarypisLab/GKlib


Note:
Build GKLIB before building METIS. On METIS cmake file, You might have to set `set(GKLIB_PATH "/path-to/GKLIB/build-dir")` if GKLIB is not installed globally.