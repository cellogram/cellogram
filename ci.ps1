param(
    # Force using MSVC
    [switch]
    $ForceMSVC,
    # Remove the build directory
    [switch]
    $Clean
)

$ErrorActionPreference = "Stop"

$cmake = (Get-Command cmake).Source

$source_dir = $PSScriptRoot
$bin_dir = Join-Path $source_dir build

if ($Clean -and (Test-Path $bin_dir)) {
    Write-Host "Removing old directory $bin_dir"
    Remove-Item $bin_dir -Recurse -Force
}

function Check-ExitCode {
    if ($LASTEXITCODE) {
        throw "Command failed [$LASTEXITCODE]"
    }
}

if ($ForceMSVC) {
    $env:CC = "cl"
    $env:CXX = "cl"
}

& $cmake -GNinja "-H$source_dir" "-B$bin_dir"
Check-ExitCode

& $cmake --build $bin_dir
Check-ExitCode
