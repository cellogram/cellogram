// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2018 Marco Tarini
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.

// #include <igl/orient2D.h>

template<typename Scalar>
IGL_INLINE short igl::orient2D(
    const Scalar A[2],
    const Scalar B[2],
    const Scalar C[2])
{
    auto sign = [](Scalar s){return (s<0)?-1:(s>0)?+1:0;};

    Scalar ACx = A[0]-C[0], ACy = A[1]-C[1];
    Scalar CBx = C[0]-B[0], CBy = C[1]-B[1];
    return sign( CBx*ACy - CBy*ACx );
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template short igl::orient2D( const double* pa, const double* pb, const double* pc );
template short igl::orient2D( const float* pa, const float* pb, const float* pc );
#endif
