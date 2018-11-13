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
        GIT_TAG        f2435e0d1909e3bd121f493432c169379fd5698f
    )
endfunction()


## pngwriter
function(cellogram_download_pngwriter)
    if(CELLOGRAM_USE_HTTPS)
        set(PNGWRITER_URL https://github.com/pngwriter/pngwriter.git)
    else()
        set(PNGWRITER_URL git@github.com:pngwriter/pngwriter.git)
    endif()

    cellogram_download_project(pngwriter
        GIT_REPOSITORY ${PNGWRITER_URL}
        GIT_TAG        79bb3ef5ebd91996bc4cd02cdc031c5cd560b897
    )
endfunction()
