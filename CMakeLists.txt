################################################################################
cmake_minimum_required(VERSION 3.1)
project(WildMeshing)
################################################################################


set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# option(FLOAT_TETWILD_WITH_EXACT_ENVELOPE    "Use exact envelope"        ON)
option(LIBIGL_WITH_PREDICATES        "Use exact predicates"         ON)


# Color output
include(UseColors)

# Prepend function
include(PrependCurrentPath)


# Extra warnings
include(Warnings)

# Use C++11/14
include(CXXFeatures)

# Sort projects inside the solution
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Generate position independent code by default
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

################################################################################


# Setup dependencies
include(WildMeshingDependencies)

# triwild
wildmeshing_download_triwild()
add_subdirectory(${THIRD_PARTY_DIR}/triwild)


# tetwild
wildmeshing_download_tetwild()
add_subdirectory(${THIRD_PARTY_DIR}/tetwild)


# pybind11
wildmeshing_download_pybind11()
add_subdirectory(${THIRD_PARTY_DIR}/pybind11)


#for testing purpose
wildmeshing_download_data()


################################################################################
# Subdirectories
################################################################################
add_library(wildmeshing MODULE
	src/binding.cpp
	src/Utils.cpp src/Utils.hpp
	src/triangulate_data.hpp src/triangulate_data.cpp
	src/triangulate.cpp src/triangulate.hpp
	src/tetrahedralize.cpp src/tetrahedralize.hpp)
target_include_directories(wildmeshing PRIVATE ${THIRD_PARTY_DIR}/triwild/src)
target_compile_features(wildmeshing PUBLIC cxx_std_14)


target_link_libraries(wildmeshing PRIVATE triwild_lib)
target_link_libraries(wildmeshing PRIVATE FloatTetwild)
target_link_libraries(wildmeshing PRIVATE pybind11::module)
set_target_properties(wildmeshing PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}" SUFFIX "${PYTHON_MODULE_EXTENSION}")




