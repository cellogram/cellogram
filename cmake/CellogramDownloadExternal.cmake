################################################################################
include(DownloadProject)

# Shortcut function
function(cellogram_download_project name)
    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${THIRD_PARTY_DIR}/${name}
        DOWNLOAD_DIR ${THIRD_PARTY_DIR}/.cache/${name}
        ${ARGN}
    )
endfunction()

################################################################################

## Polyfem
function(cellogram_download_polyfem)
    if(CELLOGRAM_USE_HTTPS)
        set(POLYFEM_URL https://github.com/geometryprocessing/polyfem.git)
    else()
        set(POLYFEM_URL git@github.com:geometryprocessing/polyfem.git)
    endif()

    cellogram_download_project(polyfem
        GIT_REPOSITORY ${POLYFEM_URL}
        GIT_TAG        f3b9f6bd5693887a801d58c72e739f5e54287304
    )
endfunction()
