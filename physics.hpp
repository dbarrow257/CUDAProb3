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

#ifndef CUDAPROB3_PHYSICS_HPP
#define CUDAPROB3_PHYSICS_HPP

#include "constants.hpp"
#include "math.hpp"

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <omp.h>


/*
 * This file contains the Barger et al physics which are used by Prob3++ to compute oscillation probabilities.
 *
 * Core function to loop over energys and cosine is function:
 *
 * template<typename FLOAT_T>
 * __host__ __device__
 * void calculate(NeutrinoType type, const FLOAT_T* const cosinelist, int n_cosines, const FLOAT_T* const energylist, int n_energies,
 *                       const FLOAT_T* const radii, const FLOAT_T* const rhos, const int* const maxlayers, FLOAT_T ProductionHeightinCentimeter, FLOAT_T* const result)
 *
 * It can either be called directly on the CPU, or on the GPU via kernel
 *
 * template<typename FLOAT_T>
 * __global__
 * void calculateKernel(NeutrinoType type, const FLOAT_T* const cosinelist, int n_cosines, const FLOAT_T* const energylist, int n_energies,
 *                       const FLOAT_T* const radii, const FLOAT_T* const rhos, const int* const maxlayers, FLOAT_T ProductionHeightinCentimeter, FLOAT_T* const result)
 *
 *
 * Both host and device code is combined in function void calculate(..), such that only one function has to be maintained for host and device.
 *
 *
 * Before using function void calculate(..) (or the kernel), neutrino mixing matrix and neutrino mass differences have to be set.
 * Use
 *
 * template<typename FLOAT_T>
 * void setMixMatrix(math::ComplexNumber<FLOAT_T>* U);
 *
 * and
 *
 * template<typename FLOAT_T>
 * void setMassDifferences(FLOAT_T* dm);
 *
 * before GPU calculation.
 *
 * Use
 *
 * template<typename FLOAT_T>
 * void setMixMatrix_host(math::ComplexNumber<FLOAT_T>* U);
 *
 * and
 *
 * template<typename FLOAT_T>
 * void setMassDifferences_host(FLOAT_T* dm);
 *
 * before CPU calculation.
 *
 *
 *
 *
 * NVCC macro __CUDA_ARCH__ is used for gpu exclusive code inside __host__ __device__ functions
 *
 */





// in device code, we need to access the device global constants instead of host global constants
#ifdef __CUDA_ARCH__

#define U(i,j) ((math::ComplexNumber<FLOAT_T>*)cudaprob3::physics::mix_data_device)[( i * 3 + j)]
#define DM(i,j) ((FLOAT_T*)cudaprob3::physics::mass_data_device)[( i * 3 + j)]
#define AXFAC(a,b,c,d,e) ((FLOAT_T*)cudaprob3::physics::A_X_factor_device)[a * 3 * 3 * 3 * 4 + b * 3 * 3 * 4 + c * 3 * 4 + d * 4 + e]
#define ORDER(i) cudaprob3::physics::mass_order_device[i]

#else

#define U(i,j) ((math::ComplexNumber<FLOAT_T>*)cudaprob3::physics::mix_data)[( i * 3 + j)]
#define DM(i,j) ((FLOAT_T*)cudaprob3::physics::mass_data)[( i * 3 + j)]
#define AXFAC(a,b,c,d,e) ((FLOAT_T*)cudaprob3::physics::A_X_factor)[a * 3 * 3 * 3 * 4 + b * 3 * 3 * 4 + c * 3 * 4 + d * 4 + e]
#define ORDER(i) cudaprob3::physics::mass_order[i]
#endif


namespace cudaprob3{

        namespace physics{

            /*
            * Constant global data
            */

            #ifdef __NVCC__
                __constant__ double mix_data_device [9 * sizeof(math::ComplexNumber<double>)] ;
                __constant__ double mass_data_device[9];
                __constant__ double A_X_factor_device[81 * 4]; //precomputed factors which only depend on the mixing matrix for faster calculation
                __constant__ int mass_order_device[3];
            #endif

            static double mix_data [9 * sizeof(math::ComplexNumber<double>)] ;
            static double mass_data[9];
            static double A_X_factor[81 * 4]; //precomputed factors for faster calculation
            static int mass_order[3];

            /*
             * Set global 3x3 pmns mixing matrix
             */
            template<typename FLOAT_T>
            void setMixMatrix(math::ComplexNumber<FLOAT_T>* U){
                memcpy((FLOAT_T*)mix_data, U, sizeof(math::ComplexNumber<FLOAT_T>) * 9);

                //precomputed factors for faster calculation
                for (int n=0; n<3; n++) {
                    for (int m=0; m<3; m++) {
                        for (int i=0; i<3; i++) {
                            for (int j=0; j<3; j++) {
                                AXFAC(n,m,i,j,0) = U[n * 3 + i].re * U[m * 3 + j].re + U[n * 3 + i].im * U[m * 3 + j].im;
                                AXFAC(n,m,i,j,1) = U[n * 3 + i].re * U[m * 3 + j].im - U[n * 3 + i].im * U[m * 3 + j].re;
                                AXFAC(n,m,i,j,2) = U[n * 3 + i].im * U[m * 3 + j].im + U[n * 3 + i].re * U[m * 3 + j].re;
                                AXFAC(n,m,i,j,3) = U[n * 3 + i].im * U[m * 3 + j].re - U[n * 3 + i].re * U[m * 3 + j].im;
                            }
                        }
                    }
                }
                #ifdef __NVCC__
                    //copy to constant memory on GPU
                    cudaMemcpyToSymbol(mix_data_device, U, sizeof(math::ComplexNumber<FLOAT_T>) * 9, 0, H2D); CUERR;
                    cudaMemcpyToSymbol(A_X_factor_device, A_X_factor, sizeof(FLOAT_T) * 81 * 4, 0, H2D); CUERR;
                #endif
            }

            /*
             * Set global 3x3 pmns mixing matrix on host only
             */
            template<typename FLOAT_T>
            void setMixMatrix_host(math::ComplexNumber<FLOAT_T>* U){
                memcpy((FLOAT_T*)mix_data, U, sizeof(math::ComplexNumber<FLOAT_T>) * 9);

                //precomputed factors for faster calculation
                for (int n=0; n<3; n++) {
                    for (int m=0; m<3; m++) {
                        for (int i=0; i<3; i++) {
                            for (int j=0; j<3; j++) {
                                AXFAC(n,m,i,j,0) = U[n * 3 + i].re * U[m * 3 + j].re + U[n * 3 + i].im * U[m * 3 + j].im;
                                AXFAC(n,m,i,j,1) = U[n * 3 + i].re * U[m * 3 + j].im - U[n * 3 + i].im * U[m * 3 + j].re;
                                AXFAC(n,m,i,j,2) = U[n * 3 + i].im * U[m * 3 + j].im + U[n * 3 + i].re * U[m * 3 + j].re;
                                AXFAC(n,m,i,j,3) = U[n * 3 + i].im * U[m * 3 + j].re - U[n * 3 + i].re * U[m * 3 + j].im;
                            }
                        }
                    }
                }
            }

            /*
             * Set global 3x3 neutrino mass difference matrix
             */
            /// \brief set mass differences to constant memory
            template<typename FLOAT_T>
            void setMassDifferences(FLOAT_T* dm){
                memcpy((FLOAT_T*)mass_data, dm, sizeof(FLOAT_T) * 9);
                #ifdef __NVCC__
                cudaMemcpyToSymbol(mass_data_device, dm, sizeof(FLOAT_T) * 9 , 0, cudaMemcpyHostToDevice); CUERR;
                #endif
            }

            /*
             * Set global 3x3 neutrino mass difference matrix on host only
             */
            template<typename FLOAT_T>
            void setMassDifferences_host(FLOAT_T* dm){
                memcpy((FLOAT_T*)mass_data, dm, sizeof(FLOAT_T) * 9);
            }


            //
            template<typename FLOAT_T>
            void prepare_getMfast(NeutrinoType type) {
                FLOAT_T alphaV, betaV, gammaV, argV, tmpV;
                FLOAT_T theta0V, theta1V, theta2V;
                FLOAT_T mMatV[3];

                /* The strategy to sort out the three roots is to compute the vacuum
                * mass the same way as the "matter" masses are computed then to sort
                * the results according to the input vacuum masses
                */
                alphaV = DM(0,1) + DM(0,2);

                betaV = DM(0,1) * DM(0,2);

                gammaV = 0.0;

                /* Compute the argument of the arc-cosine */
                tmpV = alphaV*alphaV-3.0*betaV;

                /* Equation (21) */
                argV = (2.0*alphaV*alphaV*alphaV-9.0*alphaV*betaV+27.0*gammaV)/
                    (2.0*sqrt(tmpV*tmpV*tmpV));
                if (fabs(argV)>1.0) argV = argV/fabs(argV);

                /* These are the three roots the paper refers to */
                theta0V = acos(argV)/3.0;
                theta1V = theta0V-(2.0*M_PI/3.0);
                theta2V = theta0V+(2.0*M_PI/3.0);

                mMatV[0] = mMatV[1] = mMatV[2] = -(2.0/3.0)*sqrt(tmpV);
                mMatV[0] *= cos(theta0V); mMatV[1] *= cos(theta1V); mMatV[2] *= cos(theta2V);
                tmpV = DM(0,0) - alphaV/3.0;
                mMatV[0] += tmpV; mMatV[1] += tmpV; mMatV[2] += tmpV;

                /* Sort according to which reproduce the vaccum eigenstates */
                int order[3];

                for (int i=0; i<3; i++) {
                    tmpV = fabs(DM(i,0)-mMatV[0]);
                    int k = 0;

                    for (int j=1; j<3; j++) {
                        FLOAT_T tmp = fabs(DM(i,0)-mMatV[j]);
                        if (tmp<tmpV) {
                            k = j;
                            tmpV = tmp;
                        }
                    }
                    order[i] = k;
                }
                memcpy(mass_order, order, sizeof(int) * 3);

                #ifdef __NVCC__
                cudaMemcpyToSymbol(mass_order_device, order, sizeof(int) * 3, 0, cudaMemcpyHostToDevice); CUERR;
                #endif
            }

           /*
            * Return induced neutrino mass difference matrix d_dmMatMat,
            * and d_dmMatVac, which is the mass difference matrix between induced masses and vacuum masses
            *
            * The strategy to sort out the three roots is to compute the vacuum
            * mass the same way as the "matter" masses are computed then to sort
            * the results according to the input vacuum masses. Subsequently, the "matter" masses
            * are calculated, using the found sorting for vacuum masses
            *
            * In the original implementation the order of vacuum masses is computed for each bin.
            * However, the ordering of vacuum masses does only depend on the constant neutrino mixing matrix.
            * Thus, the ordering can be precomputed, which is done in prepare_getMfast
            */
            template<typename FLOAT_T>
            HOSTDEVICEQUALIFIER
            void getMfast(const FLOAT_T Enu, const FLOAT_T rho,
                const NeutrinoType type,
                FLOAT_T d_dmMatMat[][3], FLOAT_T d_dmMatVac[][3]) {

                FLOAT_T mMatU[3], mMat[3];

                /* Equations (22) fro Barger et.al.*/
                const FLOAT_T fac = [&](){
                    if(type == Antineutrino)
		      return Constants<FLOAT_T>::tworttwoGf()*Enu*rho;
                    else
		      return -Constants<FLOAT_T>::tworttwoGf()*Enu*rho;
                }();

                const FLOAT_T alpha  = fac + DM(0,1) + DM(0,2);

                const FLOAT_T beta = DM(0,1)*DM(0,2) +
                    fac*(DM(0,1)*(1.0 -
                            U(0,1).re*U(0,1).re -
                            U(0,1).im*U(0,1).im ) +
                    DM(0,2)*(1.0-
                            U(0,2).re*U(0,2).re -
                            U(0,2).im*U(0,2).im));


                const FLOAT_T gamma = fac*DM(0,1)*DM(0,2)*(U(0,0).re * U(0,0).re + U(0,0).im * U(0,0).im);

                /* Compute the argument of the arc-cosine */
                const FLOAT_T tmp = alpha*alpha-3.0*beta < 0 ? 0 : alpha*alpha-3.0*beta;

                /* Equation (21) */
                const FLOAT_T argtmp = (2.0*alpha*alpha*alpha-9.0*alpha*beta+27.0*gamma)/
                    (2.0*sqrt(tmp*tmp*tmp));
                const FLOAT_T arg = [&](){
                    if (fabs(argtmp)>1.0)
                        return argtmp/fabs(argtmp);
                    else
                        return argtmp;
                }();

                /* These are the three roots the paper refers to */
                const FLOAT_T theta0 = acos(arg)/3.0;
                const FLOAT_T theta1 = theta0-(2.0*M_PI/3.0);
                const FLOAT_T theta2 = theta0+(2.0*M_PI/3.0);

                mMatU[0] = -(2.0/3.0)*sqrt(tmp);
                mMatU[1] = -(2.0/3.0)*sqrt(tmp);
                mMatU[2] = -(2.0/3.0)*sqrt(tmp);
                mMatU[0] *= cos(theta0);
                mMatU[1] *= cos(theta1);
                mMatU[2] *= cos(theta2);
                const FLOAT_T tmp2 = DM(0,0) - alpha/3.0;
                mMatU[0] += tmp2;
                mMatU[1] += tmp2;
                mMatU[2] += tmp2;

                /* Sort according to which reproduce the vaccum eigenstates */

                UNROLLQUALIFIER
                for (int i=0; i<3; i++) {
                    mMat[i] = mMatU[ORDER(i)];
                }

                UNROLLQUALIFIER
                for (int i=0; i<3; i++) {
                    UNROLLQUALIFIER
                    for (int j=0; j<3; j++) {
                        d_dmMatMat[i][j] = mMat[i] - mMat[j];
                        d_dmMatVac[i][j] = mMat[i] - DM(j,0);
                    }
                }
            }

            /*
                Calculate the product of Eq. (11)
            */
            template<typename FLOAT_T>
            HOSTDEVICEQUALIFIER
            void get_product(const FLOAT_T L, const FLOAT_T E, const FLOAT_T rho, const FLOAT_T d_dmMatVac[][3], const FLOAT_T d_dmMatMat[][3],
                const NeutrinoType type, math::ComplexNumber<FLOAT_T> product[][3][3]){

                math::ComplexNumber<FLOAT_T> twoEHmM[3][3][3];

                const FLOAT_T fac = [&](){
                    if(type == Antineutrino)
                        return Constants<FLOAT_T>::tworttwoGf()*E*rho;
                    else
                        return -Constants<FLOAT_T>::tworttwoGf()*E*rho;
                }();

                /* Calculate the matrix 2EH-M_j */
                UNROLLQUALIFIER
                for (int n=0; n<3; n++) {
                    UNROLLQUALIFIER
                    for (int m=0; m<3; m++) {
                        twoEHmM[n][m][0].re = -fac*(U(0,n).re*U(0,m).re+U(0,n).im*U(0,m).im);
                        twoEHmM[n][m][0].im = -fac*(U(0,n).re*U(0,m).im-U(0,n).im*U(0,m).re);
                        twoEHmM[n][m][1].re = -fac*(U(0,n).re*U(0,m).re+U(0,n).im*U(0,m).im);
                        twoEHmM[n][m][1].im = -fac*(U(0,n).re*U(0,m).im-U(0,n).im*U(0,m).re);
                        twoEHmM[n][m][2].re = -fac*(U(0,n).re*U(0,m).re+U(0,n).im*U(0,m).im);
                        twoEHmM[n][m][2].im = -fac*(U(0,n).re*U(0,m).im-U(0,n).im*U(0,m).re);
                    }
                }

                UNROLLQUALIFIER
                for (int j=0; j<3; j++){
                    twoEHmM[0][0][j].re-= d_dmMatVac[j][0];
                    twoEHmM[1][1][j].re-= d_dmMatVac[j][1];
                    twoEHmM[2][2][j].re-= d_dmMatVac[j][2];
                }

                /* Calculate the product in eq.(11) of twoEHmM for j!=k */
                //memset(product, 0, 3*3*3*sizeof(math::ComplexNumber<FLOAT_T>));
                UNROLLQUALIFIER
                for (int i=0; i<3; i++) {
                    UNROLLQUALIFIER
                    for (int j=0; j<3; j++) {
                        UNROLLQUALIFIER
                        for (int k=0; k<3; k++) {
                            product[i][j][k].re = 0;
                            product[i][j][k].im = 0;
                        }
                    }
                }

                UNROLLQUALIFIER
                for (int i=0; i<3; i++) {
                    UNROLLQUALIFIER
                    for (int j=0; j<3; j++) {
                        UNROLLQUALIFIER
                        for (int k=0; k<3; k++) {
                            product[i][j][0].re +=
                                twoEHmM[i][k][1].re*twoEHmM[k][j][2].re -
                                twoEHmM[i][k][1].im*twoEHmM[k][j][2].im;
                            product[i][j][0].im +=
                                twoEHmM[i][k][1].re*twoEHmM[k][j][2].im +
                                twoEHmM[i][k][1].im*twoEHmM[k][j][2].re;

                            product[i][j][1].re +=
                                twoEHmM[i][k][2].re*twoEHmM[k][j][0].re -
                                twoEHmM[i][k][2].im*twoEHmM[k][j][0].im;
                            product[i][j][1].im +=
                                twoEHmM[i][k][2].re*twoEHmM[k][j][0].im +
                                twoEHmM[i][k][2].im*twoEHmM[k][j][0].re;

                            product[i][j][2].re +=
                                twoEHmM[i][k][0].re*twoEHmM[k][j][1].re -
                                twoEHmM[i][k][0].im*twoEHmM[k][j][1].im;
                            product[i][j][2].im +=
                                twoEHmM[i][k][0].re*twoEHmM[k][j][1].im +
                                twoEHmM[i][k][0].im*twoEHmM[k][j][1].re;
                        }

                        product[i][j][0].re /= (d_dmMatMat[0][1]*d_dmMatMat[0][2]);
                        product[i][j][0].im /= (d_dmMatMat[0][1]*d_dmMatMat[0][2]);
                        product[i][j][1].re /= (d_dmMatMat[1][2]*d_dmMatMat[1][0]);
                        product[i][j][1].im /= (d_dmMatMat[1][2]*d_dmMatMat[1][0]);
                        product[i][j][2].re /= (d_dmMatMat[2][0]*d_dmMatMat[2][1]);
                        product[i][j][2].im /= (d_dmMatMat[2][0]*d_dmMatMat[2][1]);
                    }
                }
            }
	  
	  /***********************************************************************
  getArg
  
  Transition matrix expanded as A = sum_k C_k exp(i arg_k).
  This function returns arg, whereas getC returns C.
  arg has index [k], where k is this expansion index.
	  ***********************************************************************/
          template<typename FLOAT_T>
          HOSTDEVICEQUALIFIER
	  void getArg(const FLOAT_T L, const FLOAT_T E, const FLOAT_T dmMatVac[][3], const FLOAT_T phase_offset, FLOAT_T arg[3]) {

	    /* (1/2)*(1/(h_bar*c)) in units of GeV/(eV^2-km) */
	    const FLOAT_T LoEfac = 2.534;

	    for (int k=0; k<3; k++) {
	      arg[k] = -LoEfac*dmMatVac[k][0]*L/E;
	      if ( k==2 ) arg[k] += phase_offset ;
	    }
	  }
	  
	  /***********************************************************************
  getC
  
  Transition matrix expanded as A = sum_k C_k exp(i arg_k).
  This function returns C, whereas getArg returns arg.
  C_{re,im} have indices [k][row][col] where k is this expansion index.
	  ***********************************************************************/
          template<typename FLOAT_T>
          HOSTDEVICEQUALIFIER
          void getC(FLOAT_T E, FLOAT_T rho, FLOAT_T dmMatVac[][3], FLOAT_T dmMatMat[][3],
		    cudaprob3::NeutrinoType type, FLOAT_T phase_offset, math::ComplexNumber<FLOAT_T> C[3][3][3]) {
		    
	    const int nExp = 3;
	    const int nNuFlav = 3;
	  
	    math::ComplexNumber<FLOAT_T> product[nNuFlav][nNuFlav][nExp];

            if (phase_offset == 0.0) {
	      FLOAT_T L = NAN;
              get_product(L, E, rho, dmMatVac, dmMatMat, type, product);
            }

            /* Compute the product with the mixing matrices */
            for (int iExp=0;iExp<nNuFlav;iExp++) {
              for (int iNuFlav=0;iNuFlav<nNuFlav;iNuFlav++) {

                FLOAT_T RR_nk[nNuFlav] = { 0., 0., 0. };
                FLOAT_T RI_nk[nNuFlav] = { 0., 0., 0. };
                FLOAT_T IR_nk[nNuFlav] = { 0., 0., 0. };
                FLOAT_T II_nk[nNuFlav] = { 0., 0., 0. };

                for (int i=0; i<nNuFlav; i++) {
                  for (int j=0; j<nNuFlav; j++) {
                    RR_nk[j] += U(iNuFlav,i).re * product[i][j][iExp].re;
                    RI_nk[j] += U(iNuFlav,i).re * product[i][j][iExp].im;
                    IR_nk[j] += U(iNuFlav,i).im * product[i][j][iExp].re;
                    II_nk[j] += U(iNuFlav,i).im * product[i][j][iExp].im;
                  }
                }

                for (int jNuFlav=0;jNuFlav<nNuFlav;jNuFlav++) {
                  FLOAT_T ReSum=0., ImSum=0.;

                  for (int j=0; j<nNuFlav; j++) {
                    ReSum += RR_nk[j] * U(jNuFlav,j).re;
                    ReSum += RI_nk[j] * U(jNuFlav,j).im;
                    ReSum += IR_nk[j] * U(jNuFlav,j).im;
                    ReSum -= II_nk[j] * U(jNuFlav,j).re;

                    ImSum += II_nk[j] * U(jNuFlav,j).im;
                    ImSum += IR_nk[j] * U(jNuFlav,j).re;
                    ImSum += RI_nk[j] * U(jNuFlav,j).re;
                    ImSum -= RR_nk[j] * U(jNuFlav,j).im;
                  }
		
                  C[iExp][iNuFlav][jNuFlav].re = ReSum;
                  C[iExp][iNuFlav][jNuFlav].im = ImSum;

                }
              }
            }
	        
          }

	  template<typename FLOAT_T>
            HOSTDEVICEQUALIFIER
	    void getA(const FLOAT_T L, const FLOAT_T E, const FLOAT_T rho, const FLOAT_T d_dmMatVac[][3], const FLOAT_T d_dmMatMat[][3],
		      const NeutrinoType type,  const FLOAT_T phase_offset, math::ComplexNumber<FLOAT_T> A[3][3]){
	      
	      math::ComplexNumber<FLOAT_T> X[3][3];
	      math::ComplexNumber<FLOAT_T> product[3][3][3];
	      /* (1/2)*(1/(h_bar*c)) in units of GeV/(eV^2-km) */
	      const FLOAT_T LoEfac = 2.534;
	      
	      if (phase_offset == 0.0) {
		get_product(L, E, rho, d_dmMatVac, d_dmMatMat, type, product);
	      }
	      
	      
	      /* Make the sum with the exponential factor in Eq. (11) */
	      //memset(X, 0, 3*3*sizeof(math::ComplexNumber<FLOAT_T>));
	      UNROLLQUALIFIER
		for (int i=0; i<3; i++) {
		  UNROLLQUALIFIER
		    for (int j=0; j<3; j++) {
		      X[i][j].re = 0;
		      X[i][j].im = 0;
		    }
		}
	      
	      UNROLLQUALIFIER
		for (int k=0; k<3; k++) {
		  const FLOAT_T arg = [&](){
		    if( k == 2)
		      return -LoEfac * d_dmMatVac[k][0] * L/E + phase_offset;
		    else
		      return -LoEfac * d_dmMatVac[k][0] * L/E;
		  }();
		  
#ifdef __CUDACC__
		  FLOAT_T c,s;
		  sincos(arg, &s, &c);
#else
		  const FLOAT_T s = sin(arg);
		  const FLOAT_T c = cos(arg);
#endif
		  UNROLLQUALIFIER
		    for (int i=0; i<3; i++) {
		      UNROLLQUALIFIER
			for (int j=0; j<3; j++) {
			  X[i][j].re += c*product[i][j][k].re - s*product[i][j][k].im;
			  X[i][j].im += c*product[i][j][k].im + s*product[i][j][k].re;
			}
		    }
		}	      
	      
	      /* Eq. (10)*/
	      //memset(A, 0, 3*3*2*sizeof(FLOAT_T));
	      
	      UNROLLQUALIFIER
		for (int n=0; n<3; n++) {
		  UNROLLQUALIFIER
		    for (int m=0; m<3; m++) {
		      A[n][m].re = 0;
		      A[n][m].im = 0;
		    }
		}	      
	      UNROLLQUALIFIER
		for (int n=0; n<3; n++) {
		  UNROLLQUALIFIER
		    for (int m=0; m<3; m++) {
		      UNROLLQUALIFIER
			for (int i=0; i<3; i++) {
			  UNROLLQUALIFIER
			    for (int j=0; j<3; j++) {
			      // use precomputed factors
			      A[n][m].re +=
				AXFAC(n,m,i,j,0) * X[i][j].re +
				AXFAC(n,m,i,j,1) * X[i][j].im;
			      A[n][m].im +=
				AXFAC(n,m,i,j,2) * X[i][j].im +
				AXFAC(n,m,i,j,3) * X[i][j].re;
			    }
			}
		    }
		}
	    }
	 
	  /*
	   * Get 3x3 transition amplitude Aout for neutrino with energy E travelling Len kilometers through matter of constant density rho
	   */
	  template<typename FLOAT_T>
	  HOSTDEVICEQUALIFIER
	  void get_transition_matrix(NeutrinoType nutype, FLOAT_T Enu, FLOAT_T rho, FLOAT_T Len, math::ComplexNumber<FLOAT_T> Aout[3][3], FLOAT_T phase_offset){
	    
	    FLOAT_T d_dmMatVac[3][3], d_dmMatMat[3][3];
	    getMfast(Enu, rho, nutype, d_dmMatMat, d_dmMatVac);
	    
	    getA(Len, Enu, rho, d_dmMatVac, d_dmMatMat, nutype, phase_offset, Aout);
	  }

	  //##########################################################
	  /*
	   *   Obtain transition matrix expanded as A = sum_k C_k exp(i arg_k).
	   *   Cout_{re,im} have indices [row][col][k] where k is this expansion index.
	   *   Cout only depends on nutypei, Enuf and rhof, but not on Lenf.
	   *   Similarly argout has index [k].
	   */
	  template<typename FLOAT_T>
	  HOSTDEVICEQUALIFIER
	  void get_transition_matrix_expansion(NeutrinoType nutype, FLOAT_T Enu, FLOAT_T rho, FLOAT_T Len, math::ComplexNumber<FLOAT_T> Cout[3][3][3], FLOAT_T Arg[3], FLOAT_T phase_offset)
	  {
	    FLOAT_T d_dmMatVac[3][3], d_dmMatMat[3][3];
	    getMfast(Enu, rho, nutype, d_dmMatMat, d_dmMatVac);
	    
	    getArg(Len, Enu, d_dmMatVac, phase_offset, Arg);
	    getC(Enu, rho, d_dmMatVac, d_dmMatMat, nutype, phase_offset, Cout);
	  }

            /*
                Find density in layer
            */
            template<typename FLOAT_T>
            HOSTDEVICEQUALIFIER
            FLOAT_T getDensityOfLayer(const FLOAT_T* const rhos, int layer, int max_layer){
                if(layer == 0) return 0.0;
                int i;
                if(layer <= max_layer){
                    i = layer-1;
                }else{
                    i = 2 * max_layer - layer - 1;
                }

                return rhos[i];
            }

            /*
                Find distance in layer
            */
            template<typename FLOAT_T>
            HOSTDEVICEQUALIFIER
            FLOAT_T getTraversedDistanceOfLayer(const FLOAT_T* const radii,
                                                int layer,
                                                int max_layer,
                                                FLOAT_T PathLength,
                                                FLOAT_T TotalEarthLength,
                                                FLOAT_T cosine_zenith){

                if(cosine_zenith >= 0) return PathLength;
                if(layer == 0) return PathLength - TotalEarthLength;

                int i;
                if(layer >= max_layer)
                    i = -layer - 1 + 2 * max_layer;
                else{
                    i = layer-1;
                }

                const FLOAT_T CrossThis = 2.0*sqrt( radii[i] * radii[i]  - (Constants<FLOAT_T>::REarth())*(Constants<FLOAT_T>::REarth())*( 1 - cosine_zenith*cosine_zenith ) );
                const FLOAT_T CrossNext = 2.0*sqrt( radii[i+1] * radii[i+1] - (Constants<FLOAT_T>::REarth())*((FLOAT_T)Constants<FLOAT_T>::REarth())*( 1 -cosine_zenith*cosine_zenith ) );

                if(i < max_layer - 1){
                    return 0.5*( CrossThis-CrossNext )*(Constants<FLOAT_T>::km2cm());
                }else{
                    return CrossThis*(Constants<FLOAT_T>::km2cm());
                }
            }

	  template<typename FLOAT_T>
	  HOSTDEVICEQUALIFIER
	  void calculate(NeutrinoType type,
			 const FLOAT_T* const cosinelist,
			 int n_cosines,
			 const FLOAT_T* const energylist,
			 int n_energies,
			 const FLOAT_T* const radii,
			 const FLOAT_T* const rhos,
			 const int* const maxlayers,
			 FLOAT_T ProductionHeightinCentimeter,
			 const FLOAT_T* const productionHeight_prob_list, // 20 (nBins) * 2 (nu,nubar) * 3 (e,mu,tau) * n_energies * n_cosines
			 const FLOAT_T* const productionHeight_binedges_list, // 21 (BinEdges) in cm
			 FLOAT_T* const result){
	    
            //prepare more constant data. For the kernel, this is done by the wrapper function callCalculateKernelAsync
#ifndef __CUDA_ARCH__
	    prepare_getMfast<FLOAT_T>(type);
#endif
	    
#ifdef __CUDA_ARCH__
	    // on the device, we use the global thread Id to index the data
	    const int max_energies_per_path = SDIV(n_energies, blockDim.x) * blockDim.x;
	    for(unsigned index = blockIdx.x * blockDim.x + threadIdx.x; index < n_cosines * max_energies_per_path; index += blockDim.x * gridDim.x){
	      const unsigned index_energy = index % max_energies_per_path;
	      const unsigned index_cosine = index / max_energies_per_path;
#else
	      // on the host, we use OpenMP to parallelize looping over cosines
#pragma omp parallel for schedule(dynamic)
	      for(int index_cosine = 0; index_cosine < n_cosines; index_cosine += 1){
#endif
		
		const FLOAT_T cosine_zenith = cosinelist[index_cosine];
		
		const int iLayerAtm = 0;
		const FLOAT_T TotalEarthLength =  -2.0*cosine_zenith*Constants<FLOAT_T>::REarthcm(); // in [cm]
		const int MaxLayer = maxlayers[index_cosine];

		const int nExp = 3;
		const int nNuFlav = 3;

		FLOAT_T phaseOffset = 0.;

		const int nMaxLayers = 8;
		//DB Should check this when maxlayers is set
		if (MaxLayer > nMaxLayers) {
		  printf("Error: MaxLayer=%d given by nMaxLayers=%d. Please increase nMaxLayers in physics.cc:%d\n", MaxLayer, nMaxLayers, __LINE__);
		  std::exit(-1);
		}

		math::ComplexNumber<FLOAT_T> TransitionMatrix[nNuFlav][nNuFlav];
		math::ComplexNumber<FLOAT_T> TransitionMatrixCoreToMantle[nNuFlav][nNuFlav];
		math::ComplexNumber<FLOAT_T> finalTransitionMatrix[nNuFlav][nNuFlav];
		math::ComplexNumber<FLOAT_T> TransitionTemp[nNuFlav][nNuFlav];

		math::ComplexNumber<FLOAT_T> ExpansionMatrix[nMaxLayers][nExp][nNuFlav][nNuFlav];
		FLOAT_T arg[nMaxLayers][nNuFlav];

		// DB Uncomment for debugging get_transition_matrix against get_transition_matrix_expansion
		math::ComplexNumber<FLOAT_T> TransitionMatrix_getA[nNuFlav][nNuFlav];

		FLOAT_T Prob[nNuFlav][nNuFlav];

		math::ComplexNumber<FLOAT_T> totalLenShiftFactor[nNuFlav][nNuFlav][nExp];
		FLOAT_T PathLengths[NPRODHEIGHTBINS];
		FLOAT_T darg0_ddistance[nNuFlav];

		math::ComplexNumber<FLOAT_T> Product[nExp][nNuFlav][nNuFlav];

#ifndef __CUDA_ARCH__
		for(int index_energy = 0; index_energy < n_energies; index_energy += 1){
#else
		  if(index_energy < n_energies){
#endif
		    const FLOAT_T energy = energylist[index_energy];
		    		    
		    //============================================================================================================
		    //DB Reset all the values

		    for (int iLayer=0;iLayer<nMaxLayers;iLayer++) {
		      for (int iNuFlav=0;iNuFlav<nNuFlav;iNuFlav++) {
			clear_complex_matrix(ExpansionMatrix[iLayer][iNuFlav]);
		      }
		    }

		    UNROLLQUALIFIER
		      for (int iLayer=0;iLayer<nMaxLayers;iLayer++) {
			UNROLLQUALIFIER
			  for (int iNuFlav=0;iNuFlav<nNuFlav;iNuFlav++) {
			    arg[iLayer][iNuFlav] = 0.;
			  }
		      }
		    
		    UNROLLQUALIFIER
		      for (int iNuFlav=0;iNuFlav<nNuFlav;iNuFlav++) {
			UNROLLQUALIFIER
			  for (int jNuFlav=0;jNuFlav<nNuFlav;jNuFlav++) {
			    Prob[iNuFlav][jNuFlav] = 0.;
			  }
		      }
		    
		    UNROLLQUALIFIER
		    for (int iNuFlav=0;iNuFlav<nNuFlav;iNuFlav++) {
		      darg0_ddistance[iNuFlav] = 0.;
		    }

		    for (int iExp=0;iExp<nExp;iExp++) {
		      clear_complex_matrix(Product[iExp]);
		    }

		    // set TransitionMatrixCoreToMantle to unit matrix
		    UNROLLQUALIFIER
		      for (int iNuFlav=0;iNuFlav<nNuFlav;iNuFlav++) {
			UNROLLQUALIFIER
			for (int jNuFlav=0;jNuFlav<nNuFlav;jNuFlav++) {
			    TransitionMatrixCoreToMantle[iNuFlav][jNuFlav].re = (iNuFlav == jNuFlav ? 1.0 : 0.0);
			    TransitionMatrixCoreToMantle[iNuFlav][jNuFlav].im = 0.0;
			  }
		      }

		    //============================================================================================================
		    //DB Calculate Path Lengths for given production heights

		    //DB PathLength is used to calculate the distance traversed
		    const FLOAT_T PathLength = sqrt((Constants<FLOAT_T>::REarthcm() + ProductionHeightinCentimeter )*(Constants<FLOAT_T>::REarthcm() + ProductionHeightinCentimeter)
						    - (Constants<FLOAT_T>::REarthcm()*Constants<FLOAT_T>::REarthcm())*( 1 - cosine_zenith*cosine_zenith)) - Constants<FLOAT_T>::REarthcm()*cosine_zenith;

		    UNROLLQUALIFIER
		    for (int iProductionHeight=0;iProductionHeight<NPRODHEIGHTBINS;iProductionHeight++) {
		      FLOAT_T ProdHeightInCentimeter = Constants<FLOAT_T>::km2cm() * (productionHeight_binedges_list[iProductionHeight]+productionHeight_binedges_list[iProductionHeight+1])/2.0;
		      PathLengths[iProductionHeight] = sqrt((Constants<FLOAT_T>::REarthcm() + ProdHeightInCentimeter )*(Constants<FLOAT_T>::REarthcm() + ProdHeightInCentimeter)
					     - (Constants<FLOAT_T>::REarthcm()*Constants<FLOAT_T>::REarthcm())*( 1 - cosine_zenith*cosine_zenith)) - Constants<FLOAT_T>::REarthcm()*cosine_zenith;
		    }

		    //DB Set unit matrix
		    UNROLLQUALIFIER
		      for (int iNuFlav=0;iNuFlav<nNuFlav;iNuFlav++) {
			UNROLLQUALIFIER
			  for (int jNuFlav=0;jNuFlav<nNuFlav;jNuFlav++) {
			    UNROLLQUALIFIER
			      for (int iExp=0;iExp<nExp;iExp++) {
				totalLenShiftFactor[iNuFlav][jNuFlav][iExp].re = (iNuFlav == jNuFlav ? 1.0 : 0.0);
				totalLenShiftFactor[iNuFlav][jNuFlav][iExp].im = 0.;
			      }
			  }
		      }

		    //============================================================================================================
		    //DB Loop over layers		    

		    // loop from vacuum layer to innermost crossed layer
		    for (int iLayer=0;iLayer<=MaxLayer;iLayer++){
		      const FLOAT_T distance = getTraversedDistanceOfLayer(radii, iLayer, MaxLayer, PathLength, TotalEarthLength, cosine_zenith);
		      const FLOAT_T density = getDensityOfLayer(rhos, iLayer, MaxLayer);
		     
		      //DB Uncomment for debugging get_transition_matrix against get_transition_matrix_expansion
		      get_transition_matrix(type,
					    energy,
					    density,
					    distance / Constants<FLOAT_T>::km2cm(),
					    TransitionMatrix_getA,
					    phaseOffset
					    );

		      get_transition_matrix_expansion(type,
						      energy,
						      density,
						      distance / Constants<FLOAT_T>::km2cm(),
						      ExpansionMatrix[iLayer],
						      arg[iLayer],
						      phaseOffset
						      ); 
		
		      //DB For each layer, A = sum_k C[k] exp(i a[k])
		      clear_complex_matrix(TransitionMatrix);

		      for (int iNuFlav=0;iNuFlav<nNuFlav;iNuFlav++) {
			multiply_phase_matrix(arg[iLayer][iNuFlav],ExpansionMatrix[iLayer][iNuFlav],TransitionMatrix);
		      }
		      
		      //DB Uncomment for debugging get_transition_matrix against get_transition_matrix_expansion
		      for (int iNuFlav=0;iNuFlav<nNuFlav;iNuFlav++) {
			for (int jNuFlav=0;jNuFlav<nNuFlav;jNuFlav++) {
			  
			  if ( fabs(TransitionMatrix[iNuFlav][jNuFlav].re - TransitionMatrix_getA[iNuFlav][jNuFlav].re)>1e-9 || fabs(TransitionMatrix[iNuFlav][jNuFlav].im - TransitionMatrix_getA[iNuFlav][jNuFlav].im)>1e-9 ) {
			    printf("TransitionMatrix[iNuFlav][jNuFlav].re: %4.2f \n",TransitionMatrix[iNuFlav][jNuFlav].re);
			    printf("TransitionMatrix[iNuFlav][jNuFlav].im: %4.2f \n",TransitionMatrix[iNuFlav][jNuFlav].im);
			    printf("TransitionMatrix_getA[iNuFlav][jNuFlav].re: %4.2f \n",TransitionMatrix_getA[iNuFlav][jNuFlav].re);
			    printf("TransitionMatrix_getA[iNuFlav][jNuFlav].im: %4.2f \n",TransitionMatrix_getA[iNuFlav][jNuFlav].im);

			    std::cout << "------------ Arg[i] -------------" << std::endl;
			    for (int kNuFlav=0;kNuFlav<nNuFlav;kNuFlav++) {
			      std::cout << "arg[" << kNuFlav << "]:" << arg[iLayer][kNuFlav] << std::endl;
			    }

			    std::cout << "------------ ExpansionMatrix[iLayer,iExp,kNuFlav,mNuFlav] -------------" << std::endl;
			    for (int iExp=0;iExp<nExp;iExp++) {
			      for (int kNuFlav=0;kNuFlav<nNuFlav;kNuFlav++) {
				for (int mNuFlav=0;mNuFlav<nNuFlav;mNuFlav++) {
				  std::cout << "ExpansionMatrix[" << iLayer << "," << iExp << "," << kNuFlav << "," << mNuFlav << "].re:" << ExpansionMatrix[iLayer][iExp][kNuFlav][mNuFlav].re << std::endl;
				  std::cout << "ExpansionMatrix[" << iLayer << "," << iExp << "," << kNuFlav << "," << mNuFlav << "].im:" << ExpansionMatrix[iLayer][iExp][kNuFlav][mNuFlav].im << std::endl;
				}
			      }
			    }

			    std::cout << "------------ TransitionMatrix[kNuFlav,mNuFlav] -------------" << std::endl;
			    for (int kNuFlav=0;kNuFlav<nNuFlav;kNuFlav++) {
			      for (int mNuFlav=0;mNuFlav<nNuFlav;mNuFlav++) {

				if (fabs(TransitionMatrix[kNuFlav][mNuFlav].re) < 1e-9) {
				  std::cout << "TransitionMatrix[" << kNuFlav << "," << mNuFlav << "].re:" << 0 << std::endl;
				} else {
				  std::cout << "TransitionMatrix[" << kNuFlav << "," << mNuFlav << "].re:" << TransitionMatrix[kNuFlav][mNuFlav].re << std::endl;
				}

				if (fabs(TransitionMatrix[kNuFlav][mNuFlav].im) < 1e-9) {
				  std::cout << "TransitionMatrix[" << kNuFlav << "," << mNuFlav << "].im:" << 0 << std::endl;
				} else {
				  std::cout << "TransitionMatrix[" << kNuFlav << "," << mNuFlav << "].im:" << TransitionMatrix[kNuFlav][mNuFlav].im << std::endl;
				}
			      }
			    }

			    std::cout << "------------ TransitionMatrix_getA[kNuFlav,mNuFlav] -------------" << std::endl;
			    for (int kNuFlav=0;kNuFlav<nNuFlav;kNuFlav++) {
			      for (int mNuFlav=0;mNuFlav<nNuFlav;mNuFlav++) {

				if (fabs(TransitionMatrix_getA[kNuFlav][mNuFlav].re) < 1e-9) {
				  std::cout << "TransitionMatrix_getA[" << kNuFlav << "," << mNuFlav << "].re:" << 0 << std::endl;
				} else {
				  std::cout << "TransitionMatrix_getA[" << kNuFlav << "," << mNuFlav << "].re:" << TransitionMatrix_getA[kNuFlav][mNuFlav].re << std::endl;
				}

				if (fabs(TransitionMatrix_getA[kNuFlav][mNuFlav].im) < 1e-9) {
				  std::cout << "TransitionMatrix_getA[" << kNuFlav << "," << mNuFlav << "].im:" << 0 << std::endl;
				} else {
				  std::cout << "TransitionMatrix_getA[" << kNuFlav << "," << mNuFlav << "].im:" << TransitionMatrix_getA[kNuFlav][mNuFlav].im << std::endl;
				}
			      }
			    }

			    std::exit(-1);
			  }
			}
		      }

		      if (iLayer == iLayerAtm) { // atmosphere

			copy_complex_matrix(TransitionMatrix , finalTransitionMatrix);
			
			if (distance==0.) {
			  for (int iNuFlav=0;iNuFlav<nNuFlav;iNuFlav++) {
			    darg0_ddistance[iNuFlav] = 0.;
			  }
			}
			else {
			  for (int iNuFlav=0;iNuFlav<nNuFlav;iNuFlav++) {
			    darg0_ddistance[iNuFlav] = arg[iLayerAtm][iNuFlav]/distance;
			  }
			}

		      } else if (iLayer < MaxLayer) { // not the innermost layer, can reuse current TransitionMatrix
			clear_complex_matrix( TransitionTemp );
			multiply_complex_matrix( TransitionMatrix, finalTransitionMatrix, TransitionTemp );
			copy_complex_matrix( TransitionTemp, finalTransitionMatrix );
			
			clear_complex_matrix( TransitionTemp );
			multiply_complex_matrix( TransitionMatrixCoreToMantle, TransitionMatrix, TransitionTemp );
			copy_complex_matrix( TransitionTemp, TransitionMatrixCoreToMantle );
		      } else { // innermost layer
			clear_complex_matrix( TransitionTemp );
			multiply_complex_matrix( TransitionMatrix, finalTransitionMatrix, TransitionTemp );
			copy_complex_matrix( TransitionTemp, finalTransitionMatrix );
		      }

		    }
		    
		    // calculate final transition matrix
		    clear_complex_matrix( TransitionTemp );
		    multiply_complex_matrix( TransitionMatrixCoreToMantle, finalTransitionMatrix, TransitionTemp );
		    copy_complex_matrix( TransitionTemp, finalTransitionMatrix );

		    //============================================================================================================
		    //DB Calculate totalLenShiftFactors using atmospheric layer
		    
		    //for (int iPathLength=0;iPathLength<NPRODHEIGHTBINS;iPathLength++) {
		    UNROLLQUALIFIER
		      for (int iPathLength=0;iPathLength<(20-1);iPathLength++) {
			FLOAT_T h0 = PathLengths[iPathLength];
			FLOAT_T h1 = PathLengths[iPathLength+1];
			FLOAT_T hm = (h1+h0)/2.;
			FLOAT_T hw = (h1-h0);
		      
			UNROLLQUALIFIER
			  for (int ieig0=0;ieig0<nNuFlav;ieig0++) { 
			    UNROLLQUALIFIER
			      for (int jeig0=0;jeig0<ieig0;jeig0++) { 
				FLOAT_T darg_distance = darg0_ddistance[ieig0]-darg0_ddistance[jeig0];
				
				//factor.re = 0
				//factor.im = darg_distance*hm
				//
				//exp(factor).re = exp(f.re)*cos(f.im) = cos(darg_distance*hm)
				//exp(factor).im = exp(f.re)*sin(f.im) = sin(darg_distance*hm)
				//
				//sinc(0.5 * darg_distance * hw) * exp(factor).re = sinc(0.5 * darg_distance * hw) * cos(darg_distance*hm)
				//sinc(0.5 * darg_distance * hw) * exp(factor).re = sinc(0.5 * darg_distance * hw) * sin(darg_distance*hm)
				
				math::ComplexNumber<FLOAT_T> sinc_exp_factor; 
				FLOAT_T Sinc_Arg = 0.5 * darg_distance * hw;
				
				sinc_exp_factor.re = cudaprob3::math::defined_sinc(Sinc_Arg) * cos(darg_distance * hm);
				sinc_exp_factor.im = cudaprob3::math::defined_sinc(Sinc_Arg) * sin(darg_distance * hm);
				
				UNROLLQUALIFIER
				  for (int iNuFlav=0;iNuFlav<nNuFlav;iNuFlav++) { //In flav
				    
				    int ProbIndex = type*nNuFlav*n_energies*n_cosines*NPRODHEIGHTBINS + iNuFlav*n_energies*n_cosines*NPRODHEIGHTBINS
				      + index_energy*n_cosines*NPRODHEIGHTBINS + index_cosine*NPRODHEIGHTBINS + iPathLength;
				    
				    /*
				    //Temporary debug : 0.05 = 1.0/20 bins
				    if (productionHeight_prob_list[ProbIndex] != 0.05) {
				      std::cout << "Probability:" << productionHeight_prob_list[ProbIndex] << std::endl;
				      std::cout << "Index:" << ProbIndex << std::endl;
				      
				      std::cout << "type:" << type << std::endl;
				      std::cout << "iNuFlav:" << iNuFlav << std::endl;
				      std::cout << "index_energy:" << index_energy << std::endl;
				      std::cout << "index_cosine:" << index_cosine << std::endl;
				      std::cout << "iPathLength:" << iPathLength << std::endl;
				      std::cout << "NPRODHEIGHTBINS:" << NPRODHEIGHTBINS << std::endl;
				      
				      
				      std::cout << "nNuFlav:" << nNuFlav << std::endl;
				      
				      throw;
				    }
				    */

				    totalLenShiftFactor[ieig0][jeig0][iNuFlav].re += productionHeight_prob_list[ProbIndex] * sinc_exp_factor.re;
				    totalLenShiftFactor[ieig0][jeig0][iNuFlav].im += productionHeight_prob_list[ProbIndex] * sinc_exp_factor.im;
				    
				    totalLenShiftFactor[jeig0][ieig0][iNuFlav].re += productionHeight_prob_list[ProbIndex] * sinc_exp_factor.re;
				    totalLenShiftFactor[jeig0][ieig0][iNuFlav].im += productionHeight_prob_list[ProbIndex] * sinc_exp_factor.im;
				  }
			      }
			  }
		      }
		    
		    //============================================================================================================
		    //DB Calculate Probability from finalTransitionMatrix

		    /*
		    UNROLLQUALIFIER
		      for (int iNuFlav=0;iNuFlav<nNuFlav;iNuFlav++) { //Flavour after osc
			UNROLLQUALIFIER
			  for (int jNuFlav=0;jNuFlav<nNuFlav;jNuFlav++) { //Flavour before osc
			    Prob[jNuFlav][iNuFlav] += finalTransitionMatrix[jNuFlav][iNuFlav].re * finalTransitionMatrix[jNuFlav][iNuFlav].re + finalTransitionMatrix[jNuFlav][iNuFlav].im * finalTransitionMatrix[jNuFlav][iNuFlav].im;
			  }
		      }
		    */

		    //DB B[iExp]       = A(layer = nLayers-1)  * C(layer=0)[iExp]
		    //   Product[iExp] = finalTransitionMatrix * ExpansionMatrix[iLayerAtm][iExp]
		    //
		    //   math::ComplexNumber<FLOAT_T> Product[nExp][nNuFlav][nNuFlav];
		    //   math::ComplexNumber<FLOAT_T> finalTransitionMatrix[nNuFlav][nNuFlav];
		    //   math::ComplexNumber<FLOAT_T> ExpansionMatrix[nMaxLayers][nExp][nNuFlav][nNuFlav];
		    for (int iExp=0;iExp<nExp;iExp++) {
		      multiply_complex_matrix(finalTransitionMatrix,ExpansionMatrix[iLayerAtm][iExp],Product[iExp]);
		    }

		    //If I comment out the below, I have matching oscillograms to prior implementation
		    //From here
		    UNROLLQUALIFIER
		      for (int iExp=0;iExp<nExp;iExp++) {
			
			UNROLLQUALIFIER
			for (int jNuFlav=0;jNuFlav<nNuFlav;jNuFlav++) { //Flavour before osc 
			  UNROLLQUALIFIER
			  for (int iNuFlav=0;iNuFlav<nNuFlav;iNuFlav++) { //Flavour after osc 
			    Prob[iNuFlav][jNuFlav] += Product[iExp][iNuFlav][jNuFlav].re * Product[iExp][iNuFlav][jNuFlav].re + Product[iExp][iNuFlav][jNuFlav].im * Product[iExp][iNuFlav][jNuFlav].im;
			  }
			}

			UNROLLQUALIFIER
			  for (int jExp=0;jExp<iExp;jExp++) { //Expansion * Expansion terms
			    UNROLLQUALIFIER
			      for (int jNuFlav=0;jNuFlav<nNuFlav;jNuFlav++) { //Flavour before osc
				UNROLLQUALIFIER
				  for (int iNuFlav=0;iNuFlav<nNuFlav;iNuFlav++) { //Flavour after osc
				    Prob[iNuFlav][jNuFlav] +=  2. * Product[jExp][iNuFlav][jNuFlav].re * Product[iExp][iNuFlav][jNuFlav].re * totalLenShiftFactor[iExp][jExp][jNuFlav].re
				                            +  2. * Product[jExp][iNuFlav][jNuFlav].im * Product[iExp][iNuFlav][jNuFlav].im * totalLenShiftFactor[iExp][jExp][jNuFlav].re
				                            +  2. * Product[jExp][iNuFlav][jNuFlav].im * Product[iExp][iNuFlav][jNuFlav].re * totalLenShiftFactor[iExp][jExp][jNuFlav].im
				                            -  2. * Product[jExp][iNuFlav][jNuFlav].re * Product[iExp][iNuFlav][jNuFlav].im * totalLenShiftFactor[iExp][jExp][jNuFlav].im;
				  }
			      }
			  }
		      }
		    //To here

		    //============================================================================================================
		    //DB Fill Arrays

		    UNROLLQUALIFIER
		      for (int iNuFlav=0;iNuFlav<nNuFlav;iNuFlav++){ //Flavour before osc
                        UNROLLQUALIFIER
                          for (int jNuFlav=0;jNuFlav<nNuFlav;jNuFlav++){  //Flavour after osc
#ifdef __CUDA_ARCH__
			    const unsigned long long resultIndex = (unsigned long long)(n_energies) * (unsigned long long)(index_cosine) + (unsigned long long)(index_energy);
			    result[resultIndex + (unsigned long long)(n_energies) * (unsigned long long)(n_cosines) * (unsigned long long)((iNuFlav * nNuFlav + jNuFlav))] = Prob[jNuFlav][iNuFlav];
#else
			    const unsigned long long resultIndex = (unsigned long long)(index_cosine) * (unsigned long long)(n_energies) * (unsigned long long)(9)
			      + (unsigned long long)(index_energy) * (unsigned long long)(9);
			    result[resultIndex + (unsigned long long)((iNuFlav * nNuFlav + jNuFlav))] = Prob[jNuFlav][iNuFlav];
#endif
			  }
		      }
		  }
		}
	      }
	  
            #ifdef __NVCC__
            template<typename FLOAT_T>
            KERNEL
            __launch_bounds__( 64, 8 )
            void calculateKernel(NeutrinoType type,
                                const FLOAT_T* const cosinelist,
                                int n_cosines,
                                const FLOAT_T* const energylist,
                                int n_energies,
                                const FLOAT_T* const radii,
                                const FLOAT_T* const rhos,
                                const int* const maxlayers,
				FLOAT_T ProductionHeightinCentimeter,
				const FLOAT_T* const productionHeight_prob_list,
				const FLOAT_T* const productionHeight_binedges_list,
                                FLOAT_T* const result){

	      calculate(type, cosinelist, n_cosines, energylist, n_energies, radii, rhos, maxlayers, ProductionHeightinCentimeter, productionHeight_prob_list, productionHeight_binedges_list, result);
            }

            template<typename FLOAT_T>
            void callCalculateKernelAsync(dim3 grid,
                                        dim3 block,
                                        cudaStream_t stream,
                                        NeutrinoType type,
                                        const FLOAT_T* const cosinelist,
                                        int n_cosines,
                                        const FLOAT_T* const energylist,
                                        int n_energies,
                                        const FLOAT_T* const radii,
                                        const FLOAT_T* const rhos,
					const int* const maxlayers,
				        FLOAT_T ProductionHeightinCentimeter,
		                        const FLOAT_T* const productionHeight_prob_list,
		                        const FLOAT_T* const productionHeight_binedges_list,
                                        FLOAT_T* const result){

                prepare_getMfast<FLOAT_T>(type);

                calculateKernel<FLOAT_T><<<grid, block, 0, stream>>>(type, cosinelist, n_cosines, energylist, n_energies, radii, rhos, maxlayers, ProductionHeightinCentimeter, productionHeight_prob_list, productionHeight_binedges_list,  result);
                CUERR;
	    }
            #endif

	  }// namespace physics

} // namespace cudaprob3





#endif
