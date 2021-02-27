////////////////////////////////////////////////////////////////////////////////
// Keep this file empty, and implement unit tests in separate compilation units!
////////////////////////////////////////////////////////////////////////////////

#define CATCH_CONFIG_MAIN
#include <catch.hpp>
#include <zebrafish/Cylinder.h>
#include <zebrafish/autodiff.h>

DECLARE_DIFFSCALAR_BASE();
double zebrafish::cylinder::alpha{0.5};
double zebrafish::cylinder::K{sqrt(2)};
double zebrafish::cylinder::H{2.5};
