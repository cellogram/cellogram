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
    set(POLYFEM_URL https://github.com/polyfem/polyfem.git)

    cellogram_download_project(polyfem
        GIT_REPOSITORY ${POLYFEM_URL}
        GIT_TAG        85b6d275b4f7443344dcde4113ab410310ee1786
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
