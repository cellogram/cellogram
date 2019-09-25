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

## tinytiff
function(cellogram_download_tinytiff)
    cellogram_download_project(TinyTIFF
        GIT_REPOSITORY https://github.com/jkriege2/TinyTIFF.git
        GIT_TAG        e392a3e05bed63686cc536b6365f572a67c2e63f
    )
endfunction()

## gsl
function(cellogram_download_gsl)
    cellogram_download_project(gsl
        GIT_REPOSITORY https://github.com/jdumas/gsl.git
        GIT_TAG        b16754c5d0afe6010c78d0aecee357e0d5c90cb9
    )
endfunction()

## mmg
function(cellogram_download_mmg)
    cellogram_download_project(mmg
        GIT_REPOSITORY https://github.com/MmgTools/mmg.git
        GIT_TAG        88e2dd6cc773c43141b137fd0972c0eb2f4bbd2a
    )
endfunction()
