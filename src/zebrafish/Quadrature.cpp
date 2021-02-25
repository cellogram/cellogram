#include <zebrafish/Quadrature.h>

namespace zebrafish {

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// struct quadrature

void quadrature::LoadDiskParas(int diskQuadMethod) {

    switch (diskQuadMethod) {
        #include <zebrafish/Quad.ipp>

        default:
            assert(false);
    }

    numDiskPts = xyArray.rows();
}


quadrature::quadrature(int diskQuadMethod) {

    // Quadrature method (by default use #6)
    LoadDiskParas(diskQuadMethod);
    // height weight array (1D gaussian quadrature)
    heightLocArray.resize(4, 1);
    heightWeightArray.resize(4, 1);
    heightLocArray << -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526;
    heightWeightArray << 0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538;
    // heightLayers
    heightLayers = 4;  // NOTE: this is not the cylinder height. This is how many different depths we sample in height.
}


quadrature::~quadrature() {

}

}  // namespace zebrafish
