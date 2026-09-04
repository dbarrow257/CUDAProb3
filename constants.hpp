/*
This file is part of CUDAProb3++.

CUDAProb3++ is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

CUDAProb3++ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with CUDAProb3++.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CUDAPROB3_CONSTANTS_HPP
#define CUDAPROB3_CONSTANTS_HPP

#include "hpc_helpers.cuh"

#include <cmath>

namespace cudaprob3{
    static double EarthRadius_h = 6371.0;

    #ifdef __CUDACC__
    __device__ __constant__ double EarthRadius_d = 6371.0;
    #endif

    template<typename FLOAT_T>
    struct Constants{
        HOSTDEVICEQUALIFIER
        static constexpr FLOAT_T tworttwoGf(){ return 1.52588e-4; }

        HOSTDEVICEQUALIFIER
        static constexpr FLOAT_T km2cm(){ return 1.0e5; }

        static void SetEarthRadius(FLOAT_T EarthRadius_)
        {
            EarthRadius_h = double(EarthRadius_);
            #ifdef __CUDACC__
            cudaMemcpyToSymbol(EarthRadius_d, &EarthRadius_h, sizeof(EarthRadius_h)); CUERR;
            #endif
        }

        HOSTDEVICEQUALIFIER
        static FLOAT_T REarth(){
            #ifdef __CUDA_ARCH__
            return FLOAT_T(EarthRadius_d);
            #else
            return FLOAT_T(EarthRadius_h);
            #endif
        }

        HOSTDEVICEQUALIFIER
        static FLOAT_T REarthcm(){ return REarth() * km2cm(); }

        HOSTDEVICEQUALIFIER
        static FLOAT_T REarthcm2(){ return REarthcm() * REarthcm(); }

        HOSTDEVICEQUALIFIER
        static constexpr FLOAT_T TwoPiOverThree() { return 2.0 * M_PI / 3.0; }

        HOSTDEVICEQUALIFIER
        static constexpr FLOAT_T density_convert(){ return 0.5; }

        HOSTDEVICEQUALIFIER
        static constexpr int MaxProdHeightBins(){ return 28; }

        HOSTDEVICEQUALIFIER
        static constexpr int MaxNLayers(){ return 11; }
      
        HOSTDEVICEQUALIFIER
        static constexpr FLOAT_T Epsilon(){ return 1e-6; }
    };
}


#endif
