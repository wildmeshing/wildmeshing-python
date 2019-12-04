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
        GIT_TAG         90e45669acdf8c453f287da62757a8acf07fbb57
    )
endfunction()

function(wildmeshing_download_tetwild)
    wildmeshing_download_project(tetwild
        GIT_REPOSITORY  https://github.com/wildmeshing/fTetWild
        GIT_TAG         cdca6ec6059be0b585377b19aae5fa04643262e9
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
        GIT_TAG        2a255e554a8e8410144514d127753c8263a3bbee
    )
endfunction()