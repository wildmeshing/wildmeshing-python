################################################################################
cmake_minimum_required(VERSION 3.1)
project(WildMeshing)
################################################################################


set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

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

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	foreach(config ${CMAKE_CONFIGURATION_TYPES})
		string(TOUPPER ${config} config)
		string(REPLACE /MD /MT CMAKE_C_FLAGS_${config} "${CMAKE_C_FLAGS_${config}}")
		string(REPLACE /MD /MT CMAKE_CXX_FLAGS_${config} "${CMAKE_CXX_FLAGS_${config}}")
	endforeach()
endif()

################################################################################


# Setup dependencies
include(WildMeshingDependencies)

wildmeshing_download_svg_parser()

# triwild
wildmeshing_download_triwild()
add_subdirectory(${THIRD_PARTY_DIR}/triwild)


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
add_library(wildmeshing MODULE src/binding.cpp)
target_include_directories(wildmeshing PRIVATE ${THIRD_PARTY_DIR}/triwild/src)
target_link_libraries(wildmeshing PRIVATE triwild_lib)
target_link_libraries(wildmeshing PRIVATE FloatTetwild)
target_link_libraries(wildmeshing PRIVATE pybind11::module)
set_target_properties(wildmeshing PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}" SUFFIX "${PYTHON_MODULE_EXTENSION}")




