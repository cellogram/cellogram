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
        GIT_TAG        0dc1e7e536432ab2192b4ce55201a8df0662dadd
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


## tinyfiledialogs
function(cellogram_download_tinyfiledialogs)
    cellogram_download_project(tinyfiledialogs
        GIT_REPOSITORY https://git.code.sf.net/p/tinyfiledialogs/code
        GIT_TAG        511e6500fa9184923d4859e06ee9a6a4e70820c4
    )
endfunction()

## CLI11 3-Clause BSD license optional
function(cellogram_download_cli11)
    cellogram_download_project(cli11
        URL     https://github.com/CLIUtils/CLI11/archive/v1.8.0.tar.gz
        URL_MD5 5e5470abcb76422360409297bfc446ac
    )
endfunction()
