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
    cellogram_download_project(polyfem
        GIT_REPOSITORY https://github.com/polyfem/polyfem.git
        GIT_TAG        e4c857dbae4d032ffb980255326c3d0177b2223a
    )
endfunction()


## pngwriter
function(cellogram_download_pngwriter)
    cellogram_download_project(pngwriter
        GIT_REPOSITORY https://github.com/pngwriter/pngwriter.git
        GIT_TAG        79bb3ef5ebd91996bc4cd02cdc031c5cd560b897
    )
endfunction()
