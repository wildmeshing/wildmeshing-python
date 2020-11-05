################################################################################
include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
    set(WILDMESHING_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
    set(WILDMESHING_EXTRA_OPTIONS "")
endif()

# Shortcut function
function(wildmeshing_download_project name)
    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${WILDMESHING_EXTERNAL}/${name}
        DOWNLOAD_DIR ${WILDMESHING_EXTERNAL}/.cache/${name}
        QUIET
        ${WILDMESHING_EXTRA_OPTIONS}
        ${ARGN}
    )
endfunction()

################################################################################


function(wildmeshing_download_triwild)
    wildmeshing_download_project(triwild
        GIT_REPOSITORY  https://github.com/wildmeshing/TriWild
        GIT_TAG         a69ef125e5c6e336f1bd982afadefbc87a92398e
    )
endfunction()

function(wildmeshing_download_tetwild)
    wildmeshing_download_project(tetwild
        GIT_REPOSITORY  https://github.com/wildmeshing/fTetWild
        GIT_TAG         9e1bb27188a5780fe42398403a60953aef25ebb4
    )
endfunction()

function(wildmeshing_download_pybind11)
    wildmeshing_download_project(pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11
        GIT_TAG        085a29436a8c472caaaf7157aa644b571079bcaa
    )
endfunction()

# data
function(wildmeshing_download_data)
    wildmeshing_download_project(data
        GIT_REPOSITORY https://github.com/wildmeshing/data
        GIT_TAG        1484054abbac36e9c8340c3b32d87ad6eee45016
    )
endfunction()