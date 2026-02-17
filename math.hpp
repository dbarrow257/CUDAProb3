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

#ifndef CUDAPROB3_MATH_HPP
#define CUDAPROB3_MATH_HPP

#include "hpc_helpers.cuh"
#include "constants.hpp"
#include <cmath>

namespace cudaprob3{

    namespace math{

        template<typename T>
        struct ComplexNumber{
            T re;
            T im;
        };

        template<typename T>
        HOSTDEVICEQUALIFIER
        constexpr T ct_sqr(T x){
            return x * x;
        }

        template<typename T>
        HOSTDEVICEQUALIFIER
        constexpr T ct_cube(T x){
            return x * x * x;
        }

        template<typename FLOAT_T>
        HOSTDEVICEQUALIFIER
        FLOAT_T defined_sinc(FLOAT_T A) {
            if (abs(A) >= Constants<FLOAT_T>::Epsilon()) {
                return sin(A)/A;
            } else {
                const FLOAT_T A2 = A*A;
                return FLOAT_T(1) - A2/6. + A2*A2/120.;
            }
        }
      
        template<typename FLOAT_T>
        HOSTDEVICEQUALIFIER
        void multiply_phase_matrix(const FLOAT_T Phase, const ComplexNumber<FLOAT_T> A[3][3], ComplexNumber<FLOAT_T> B[3][3]) {

        #ifdef __CUDACC__
            FLOAT_T c,s;
            sincos(Phase, &s, &c);
        #else
            const FLOAT_T s = sin(Phase);
            const FLOAT_T c = cos(Phase);
        #endif
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    const FLOAT_T ar = A[i][j].re;
                    const FLOAT_T ai = A[i][j].im;

                    B[i][j].re += c * ar - s * ai;
                    B[i][j].im += c * ai + s * ar;
                }
            }
        }

        /*
        *   multiply complex 3x3 matrix
        *        C = A X B
        */
        template<typename FLOAT_T>
        HOSTDEVICEQUALIFIER
        void multiply_complex_matrix(const ComplexNumber<FLOAT_T> A[3][3], const ComplexNumber<FLOAT_T> B[3][3], ComplexNumber<FLOAT_T> C[3][3]){

            for (int i=0; i<3; i++) {

                for (int j=0; j<3; j++) {
                    FLOAT_T cr = 0.;
                    FLOAT_T ci = 0.;
                    for (int k=0; k<3; k++) {
                        const FLOAT_T ar = A[i][k].re;
                        const FLOAT_T ai = A[i][k].im;
                        const FLOAT_T br = B[k][j].re;
                        const FLOAT_T bi = B[k][j].im;

                        cr += ar*br - ai*bi;
                        ci += ai*br + ar*bi;
                    }
                    C[i][j].re += cr;
                    C[i][j].im += ci;
                }
            }
        }

        /*
         *   multiply complex 3x3 matrix
         *        C = A X B
         */
        template<typename FLOAT_T>
        HOSTDEVICEQUALIFIER
        void multiply_complex_matrix_cleared(const ComplexNumber<FLOAT_T> A[3][3], const ComplexNumber<FLOAT_T> B[3][3], ComplexNumber<FLOAT_T> C[3][3]){
            for (int i=0; i<3; ++i) {
                for (int j=0; j<3; ++j) {
                    FLOAT_T cr = 0.;
                    FLOAT_T ci = 0.;
                    for (int k=0; k<3; ++k) {
                        const FLOAT_T ar = A[i][k].re;
                        const FLOAT_T ai = A[i][k].im;
                        const FLOAT_T br = B[k][j].re;
                        const FLOAT_T bi = B[k][j].im;

                        cr += ar*br - ai*bi;
                        ci += ai*br + ar*bi;
                    }
                    C[i][j].re = cr;
                    C[i][j].im = ci;
                }
            }
        }
        /*
        *   multiply complex 3x3 matrix and 3 vector
        *        W = A X V
        */
        template<typename FLOAT_T>
        HOSTDEVICEQUALIFIER
        void multiply_complex_matvec(ComplexNumber<FLOAT_T> A[][3], ComplexNumber<FLOAT_T> V[3], ComplexNumber<FLOAT_T> W[3]){

            for(int i=0;i<3;i++) {
                W[i].re = A[i][0].re*V[0].re-A[i][0].im*V[0].im+
                    A[i][1].re*V[1].re-A[i][1].im*V[1].im+
                    A[i][2].re*V[2].re-A[i][2].im*V[2].im ;
                W[i].im = A[i][0].re*V[0].im+A[i][0].im*V[0].re+
                    A[i][1].re*V[1].im+A[i][1].im*V[1].re+
                    A[i][2].re*V[2].im+A[i][2].im*V[2].re ;
            }
        }

        /*
        *   copy complex 3x3 matrix
        *        A --> B
        */
        template<typename FLOAT_T>
        HOSTDEVICEQUALIFIER
        void copy_complex_matrix(const ComplexNumber<FLOAT_T> A[3][3], ComplexNumber<FLOAT_T> B[3][3]) {
            B[0][0].re = A[0][0].re;
            B[0][0].im = A[0][0].im;

            B[0][1].re = A[0][1].re;
            B[0][1].im = A[0][1].im;

            B[0][2].re = A[0][2].re;
            B[0][2].im = A[0][2].im;

            B[1][0].re = A[1][0].re;
            B[1][0].im = A[1][0].im;

            B[1][1].re = A[1][1].re;
            B[1][1].im = A[1][1].im;

            B[1][2].re = A[1][2].re;
            B[1][2].im = A[1][2].im;

            B[2][0].re = A[2][0].re;
            B[2][0].im = A[2][0].im;

            B[2][1].re = A[2][1].re;
            B[2][1].im = A[2][1].im;

            B[2][2].re = A[2][2].re;
            B[2][2].im = A[2][2].im;
        }

        /*
        *   clear complex 3x3 matrix
        *
        */
        template<typename FLOAT_T>
        HOSTDEVICEQUALIFIER
        void clear_complex_matrix(ComplexNumber<FLOAT_T> A[3][3]){
            A[0][0].re = 0;
            A[0][0].im = 0;
            A[0][1].re = 0;
            A[0][1].im = 0;
            A[0][2].re = 0;
            A[0][2].im = 0;
            A[1][0].re = 0;
            A[1][0].im = 0;
            A[1][1].re = 0;
            A[1][1].im = 0;
            A[1][2].re = 0;
            A[1][2].im = 0;
            A[2][0].re = 0;
            A[2][0].im = 0;
            A[2][1].re = 0;
            A[2][1].im = 0;
            A[2][2].re = 0;
            A[2][2].im = 0;
            //memset(A,0,sizeof(ComplexNumber<FLOAT_T>)*9);
        }

    }

}


#endif
