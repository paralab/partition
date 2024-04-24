add_library(metis-util metis-util.cpp)
target_include_directories(metis-util PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${METIS_INCLUDE_DIR})


target_link_libraries(metis-util PUBLIC util ${METIS_LIBRARY})
