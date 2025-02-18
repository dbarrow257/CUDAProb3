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

#ifndef CUDAPROB3_PROPAGATOR_HPP
#define CUDAPROB3_PROPAGATOR_HPP

#include "constants.hpp"
#include "types.hpp"
#include "math.hpp"


#include <algorithm>
#include <array>
#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <iostream>



namespace cudaprob3{


    /// \class Propagator
    /// \brief Abstract base class of the library which sets up input parameter on the host.
    /// Concrete implementation of calcuations is provided in derived classes
    /// @param FLOAT_T The floating point type to use for calculations, i.e float, double
    template<class FLOAT_T>
    class Propagator{
    public:
        /// \brief Constructor
        ///
        /// @param n_cosines Number cosine bins
        /// @param n_energies Number of energy bins
        Propagator(int n_cosines, int n_energies) : n_cosines(n_cosines), n_energies(n_energies){
            energyList.resize(n_energies);
            cosineList.resize(n_cosines);
            maxlayers.resize(n_cosines);
        }

        /// \brief Copy constructor
        /// @param other
        Propagator(const Propagator& other){
            *this = other;
        }

        /// \brief Move constructor
        /// @param other
        Propagator(Propagator&& other){
            *this = std::move(other);
        }

        virtual ~Propagator(){}

        /// \brief Copy assignment operator
        /// @param other
        Propagator& operator=(const Propagator& other){
            energyList = other.energyList;
            cosineList = other.cosineList;
            maxlayers = other.maxlayers;
            radii = other.radii;
            rhos = other.rhos;
            as = other.as;
            bs = other.bs;
            cs = other.cs;
	    yps = other.yps;
            coslimit = other.coslimit;
            Mix_U = other.Mix_U;
            dm = other.dm;

            ProductionHeightinCentimeter = other.ProductionHeightinCentimeter;
            isSetCosine = other.isSetCosine;
            isSetProductionHeight = other.isSetProductionHeight;
            isInit = other.isInit;

	    nProductionHeightBins = other.nProductionHeightBins;
	    useProductionHeightAveraging = other.useProductionHeightAveraging;

	    productionHeightList_prob = other.productionHeightList_prob;
	    productionHeightList_bins = other.productionHeightList_bins;
	    isSetProductionHeightArray = other.isSetProductionHeightArray;

            return *this;
        }

        /// \brief Move assignment operator
        /// @param other
        Propagator& operator=(Propagator&& other){
            energyList = std::move(other.energyList);
            cosineList = std::move(other.cosineList);
            maxlayers = std::move(other.maxlayers);
            radii = std::move(other.radii);
            rhos = std::move(other.rhos);
            as = std::move(other.as);
            bs = std::move(other.bs);
            cs = std::move(other.cs);
	    yps = std::move(other.yps);
            coslimit = std::move(other.coslimit);
            Mix_U = std::move(other.Mix_U);
            dm = std::move(other.dm);

            ProductionHeightinCentimeter = other.ProductionHeightinCentimeter;
            isSetCosine = other.isSetCosine;
            isSetProductionHeight = other.isSetProductionHeight;
            isInit = other.isInit;

	    nProductionHeightBins = other.nProductionHeightBins;
	    useProductionHeightAveraging = other.useProductionHeightAveraging;

            productionHeightList_prob = other.productionHeightList_prob;
            productionHeightList_bins = other.productionHeightList_bins;
            isSetProductionHeightArray = other.isSetProductionHeightArray;

            other.isInit = false;

            return *this;
        }

    public:

      // Returns the number of LAYER BOUNDARIES
      // i.e. number of layers is this number *MINUS ONE*
      int getNlayerBoundaries() { return n_layers; };

      void SetNumberOfProductionHeightBinsForAveraging(int nProductionHeightBins_) {
	if (nProductionHeightBins_ > Constants<FLOAT_T>::MaxProdHeightBins()) {
	  std::cerr << "Invalid number of production height averages:" << nProductionHeightBins_ << std::endl;
	  std::cerr << "Need to increase value of Constants<FLOAT_T>::MaxProdHeightBins() in $CUDAPROB3/constants.hpp" << std::endl;
	  throw std::runtime_error("SetNumberOfProductionHeightBinsForAveraging : invalid number of production height bins");
	}

	this->nProductionHeightBins = nProductionHeightBins_;

	if (this->nProductionHeightBins >= 1) {
	  this->useProductionHeightAveraging = true;
	}

	if (this->useProductionHeightAveraging == true) {
	  std::cout << "Set " << this->nProductionHeightBins << " Production height bins" << std::endl;
	} else {
	  std::cout << "Using fixed production height" << std::endl;
	}
      }

        /// \brief Set density information from arrays.
        /// \details radii_ and rhos_ must be same size. both radii_ and rhos_ must be sorted, in the same order.
        /// The density (g/cm^3) at a distance (km) from the center of the sphere between radii_[i], exclusive,
        /// and radii_[j], inclusive, i < j  is assumed to be rhos_[j]
        /// @param radii_ List of radii
        /// @param rhos_ List of densities
        /// @param yps_ List of chemical compositions
      virtual void setDensity(
          const std::vector<FLOAT_T>& radii_, 
          const std::vector<FLOAT_T>& rhos_, 
          const std::vector<FLOAT_T>& yps_){

          UsePolyDensity = false;

            if(rhos_.size() != radii_.size()){
                throw std::runtime_error("setDensity : rhos.size() != radii.size()");
            }

            if(rhos_.size() != yps_.size()){
	      throw std::runtime_error("setDensity : rhos.size() != yps.size()");
            }

            if(rhos_.size() == 0 || radii_.size() == 0 || yps_.size() == 0){
	      throw std::runtime_error("setDensity : vectors must not be empty");
            }

            bool needFlip = false;

            if(radii_.size() >= 2){
                int sign = (radii_[1] - radii_[0] > 0 ? 1 : -1);

                for(size_t i = 1; i < radii_.size(); i++){
                    if((radii_[i] - radii_[i-1]) * sign < 0)
                        throw std::runtime_error("radii order messed up");
                }

                if(sign == 1)
                    needFlip = true;
            }

            radii = radii_;
            rhos = rhos_;
	    yps = yps_;

            if(needFlip){
                std::reverse(radii.begin(), radii.end());
                std::reverse(rhos.begin(), rhos.end());
		std::reverse(yps.begin(), yps.end());
            }

            coslimit.clear();

            // first element of _Radii is largest radius
            for(size_t i=0; i < radii.size() ; i++ ) {
                // Using a cosine threshold
                FLOAT_T x = -1* sqrt( 1 - (radii[i] * radii[i] / ( Constants<FLOAT_T>::REarth()*Constants<FLOAT_T>::REarth())) );
                if ( i  == 0 ) x = 0;
                coslimit.push_back(x);
            }

            setMaxlayers();
        }

        /// \brief Set density information from arrays including polynomials for a non-constant density in each layer
        /// \details radii_ and rhos_ must be same size. both radii_ and rhos_ must be sorted, in the same order.
        /// The density (g/cm^3) at a distance (km) from the center of the sphere between radii_[i], exclusive,
        /// and radii_[j], inclusive, i < j  is assumed to be rhos_[j]
        /// @param radii_ List of radii
        /// @param a_ List of densities coefficient a
        /// @param b_ List of densities coefficient b
        /// @param c_ List of densities coefficient c
        /// @param yps_ List of chemical compositions
      virtual void setDensity(
          const std::vector<FLOAT_T>& radii_, 
          const std::vector<FLOAT_T>& a_, 
          const std::vector<FLOAT_T>& b_, 
          const std::vector<FLOAT_T>& c_, 
          const std::vector<FLOAT_T>& yps_) {

        UsePolyDensity = true;

        if(a_.size() != radii_.size()){
          throw std::runtime_error("setDensity : a.size() != radii.size()");
        }

        if(a_.size() != yps_.size()){
          throw std::runtime_error("setDensity : a.size() != yps.size()");
        }

        if(a_.size() == 0 || b_.size() == 0 || c_.size() == 0 || radii_.size() == 0 || yps_.size() == 0){
          throw std::runtime_error("setDensity : vectors must not be empty");
        }

        bool needFlip = false;

        if(radii_.size() >= 2){
          int sign = (radii_[1] - radii_[0] > 0 ? 1 : -1);

          for(size_t i = 1; i < radii_.size(); i++){
            if((radii_[i] - radii_[i-1]) * sign < 0)
              throw std::runtime_error("radii order messed up");
          }

          if(sign == 1) needFlip = true;
        }


        // Copy over the content (probably unnecessary...)
        radii = radii_;
        as = a_;
        bs = b_;
        cs = c_;
        yps = yps_;

        if(needFlip){
          std::reverse(radii.begin(), radii.end());
          std::reverse(yps.begin(), yps.end());
          std::reverse(as.begin(), as.end());
          std::reverse(bs.begin(), bs.end());
          std::reverse(cs.begin(), cs.end());
        }

        coslimit.clear();

        // first element of _Radii is largest radius
        for(size_t i=0; i < radii.size() ; i++ ) {
          // Using a cosine threshold
          FLOAT_T x = -1* sqrt( 1 - (radii[i] * radii[i] / ( Constants<FLOAT_T>::REarth()*Constants<FLOAT_T>::REarth())) );
          if ( i  == 0 ) x = 0;
          coslimit.push_back(x);
        }

        setMaxlayers();
      }

      /// \brief Set density information from file
      /// \details File must contain two columns where the first column contains the radius (km)
      /// and the second column contains the density (g/cm³).
      /// The first row must have the radius 0. The last row must have to radius of the sphere
      ///
      /// @param filename File with density information
      virtual void setDensityFromFile(const std::string& filename){
        std::ifstream file(filename);
        if(!file)
          throw std::runtime_error("could not open density file " + filename);

        std::vector<FLOAT_T> radii_temp;
        std::vector<FLOAT_T> rhos_temp;
        std::vector<FLOAT_T> yps_temp;

        // First check if the file contains rho or polynomial coefficient
        // reading the first line should suffice

        std::string line;
        int nentries_old = 0;
        if (file.is_open()) {
          while (std::getline(file, line)) {
            // Allow for comments or empty lines
            if (line[0] == '#' || line.empty()) continue;
            // Check how many entries we have per line
            int nentries = 0;
            while (line.find_first_of(" ") != std::string::npos) {
              int newpos = line.find_first_of(" ");
              std::string substring = line.substr(0, newpos);
              // Check repeated spaces
              while (line[newpos] == line[newpos+1]) newpos++;
              line = line.substr(newpos+1, line.size());
              nentries++;
            }
            nentries++;
            if (nentries_old == 0) nentries_old = nentries;
            if (nentries != nentries_old) std::cout << "Inconsitent number of entries" << std::endl;
          }
        }

        std::cout << "Found " << nentries_old << " entries in file " << filename << std::endl;
        // Reset the file reader
        file.clear();
        file.seekg(0);

        // If the file is formatted as radius, density, electron fraction
        if (nentries_old == 3) {
          FLOAT_T r;
          FLOAT_T d;
          FLOAT_T yp;
          while (file >> r >> d >> yp){
            radii_temp.push_back(r);
            rhos_temp.push_back(d);
            yps_temp.push_back(yp);
          }

          setDensity(radii_temp, rhos_temp, yps_temp);
        } 
        else if (nentries_old == 5) {

          // Coefficients of density
          std::vector<FLOAT_T> a_temp;
          std::vector<FLOAT_T> b_temp;
          std::vector<FLOAT_T> c_temp;

          if (file.is_open()) {
            while (std::getline(file, line)) {
              std::vector<std::string> entries;
              // Allow for comments or empty lines
              if (line[0] == '#' || line.empty()) continue;
              // Check how many entries we have per line
              while (line.find_first_of(" ") != std::string::npos) {
                int newpos = line.find_first_of(" ");
                std::string substring = line.substr(0, newpos);
                entries.push_back(substring);
                while (line[newpos] == line[newpos+1]) newpos++;
                line = line.substr(newpos+1, line.size());
              }
              entries.push_back(line);

              // Now push back into our main vectors
              radii_temp.push_back(std::atof(entries[0].c_str()));
              a_temp.push_back(std::atof(entries[1].c_str()));
              b_temp.push_back(std::atof(entries[2].c_str()));
              c_temp.push_back(std::atof(entries[3].c_str()));
              yps_temp.push_back(std::atof(entries[4].c_str()));
            }
          }
          //for (int i = 0; i < nentries_old; ++i) {
            //std::cout << radii_temp[i] << " " << a_temp[i] << " " << b_temp[i] << " " << c_temp[i] << " " << yps_temp[i] << std::endl;
          //}

          setDensity(radii_temp, a_temp, b_temp, c_temp, yps_temp);

        } else {
          std::cout << "Unsupported earty model in " << filename << std::endl;
          std::cout << "  Number of entries per line: " << nentries_old << std::endl;
          throw;
        }
        n_layers = radii.size();
      }

      virtual void ModifyEarthModelPoly(std::vector<FLOAT_T> list_radii, std::vector<FLOAT_T> list_weights){
        
        int nBoundaries(list_radii.size());
        int nWeights(list_weights.size());

        if(nBoundaries!=radii.size()-1 || nWeights!=radii.size()-1){
          std::cout<<"Error in ModifyEarthModel(list_radii, list_weights)"<<std::endl;
          std::cout<<"Expected "<<radii.size()-1<<" radii, got "<< nBoundaries <<std::endl;
          std::cout<<"Expected "<<radii.size()-1<<" weights, got "<< nWeights <<std::endl;
          throw;
        }

        for(int i=0;i<nBoundaries;i++){
          radii[i] = list_radii[nBoundaries-i-1];
        }

        for(int i=0;i<nWeights;i++){
          as[i]*= list_weights[nWeights-i-1];
          bs[i]*= list_weights[nWeights-i-1];
          cs[i]*= list_weights[nWeights-i-1];
        }
        as[nWeights]*= list_weights[0];
        bs[nWeights]*= list_weights[0];
        cs[nWeights]*= list_weights[0];

        setDensity(radii, as, bs, cs, yps);
      }

      virtual void ModifyEarthModel(std::vector<FLOAT_T> list_radii, std::vector<FLOAT_T> list_weights){

        int nBoundaries(list_radii.size());
        int nWeights(list_weights.size());

        if(nBoundaries!=radii.size()-1 || nWeights!=radii.size()-1){
          std::cout<<"Error in ModifyEarthModel(list_radii, list_weights)"<<std::endl;
          std::cout<<"Expected "<<radii.size()-1<<" radii, got "<< nBoundaries <<std::endl;
          std::cout<<"Expected "<<radii.size()-1<<" weights, got "<< nWeights <<std::endl;
          throw;
        }

        for(int i=0;i<nBoundaries;i++){
          radii[i] = list_radii[nBoundaries-i-1];
        }
        
        for(int i=0;i<nWeights;i++){
          rhos[i]*= list_weights[nWeights-i-1];
        }
        rhos[nWeights]*= list_weights[0];

        setDensity(radii, rhos, yps);
      }

      virtual bool PolynomialDensity(){
        return UsePolyDensity;
      }
      
      /// \brief Set mixing angles and cp phase in radians
      /// @param theta12
      /// @param theta13
      /// @param theta23
      /// @param dCP
      virtual void setMNSMatrix(FLOAT_T theta12, FLOAT_T theta13, FLOAT_T theta23, FLOAT_T dCP){

        const FLOAT_T s12 = sin(theta12);
        const FLOAT_T s13 = sin(theta13);
        const FLOAT_T s23 = sin(theta23);
        const FLOAT_T c12 = cos(theta12);
        const FLOAT_T c13 = cos(theta13);
        const FLOAT_T c23 = cos(theta23);

        const FLOAT_T sd  = sin(dCP);
        const FLOAT_T cd  = cos(dCP);

        U(0,0).re =  c12*c13;
        U(0,0).im =  0.0;
        U(0,1).re =  s12*c13;
        U(0,1).im =  0.0;
        U(0,2).re =  s13*cd;
        U(0,2).im = -s13*sd;
        U(1,0).re = -s12*c23-c12*s23*s13*cd;
        U(1,0).im =         -c12*s23*s13*sd;
        U(1,1).re =  c12*c23-s12*s23*s13*cd;
        U(1,1).im =         -s12*s23*s13*sd;
        U(1,2).re =  s23*c13;
        U(1,2).im =  0.0;
        U(2,0).re =  s12*s23-c12*c23*s13*cd;
        U(2,0).im =         -c12*c23*s13*sd;
        U(2,1).re = -c12*s23-s12*c23*s13*cd;
        U(2,1).im  =         -s12*c23*s13*sd;
        U(2,2).re =  c23*c13;
        U(2,2).im  =  0.0;
      }

      /// \brief Set neutrino mass differences (m_i_j)^2 in (eV)^2. no assumptions about mass hierarchy are made
      /// @param dm12sq
      /// @param dm23sq
      virtual void setNeutrinoMasses(FLOAT_T dm12sq, FLOAT_T dm23sq){
        FLOAT_T mVac[3];

        mVac[0] = 0.0;
        mVac[1] = dm12sq;
        mVac[2] = dm12sq + dm23sq;

        const FLOAT_T delta = 5.0e-9;
        /* Break any degeneracies */
        if (dm12sq == 0.0) mVac[0] -= delta;
        if (dm23sq == 0.0) mVac[2] += delta;

        DM(0,0) = 0.0;
        DM(1,1) = 0.0;
        DM(2,2) = 0.0;
        DM(0,1) = mVac[0]-mVac[1];
        DM(1,0) = -DM(0,1);
        DM(0,2) = mVac[0]-mVac[2];
        DM(2,0) = -DM(0,2);
        DM(1,2) = mVac[1]-mVac[2];
        DM(2,1) = -DM(1,2);
      }

      /// \brief Set the energy bins. Energies are given in GeV
      /// @param list Energy list
      virtual void setEnergyList(const std::vector<FLOAT_T>& list){
        if(list.size() != size_t(n_energies))
          throw std::runtime_error("Propagator::setEnergyList. Propagator was not created for this number of energy nodes");

        energyList = list;
      }

      /// \brief Set cosine bins. Cosines are given in radians
      /// @param list Cosine list
      virtual void setCosineList(const std::vector<FLOAT_T>& list){
        if(list.size() != size_t(n_cosines))
          throw std::runtime_error("Propagator::setCosineList. Propagator was not created for this number of cosine nodes");

        cosineList = list;

        if(isSetProductionHeight){
          setProductionHeight(ProductionHeightinCentimeter / 100000.0);
        }

        setMaxlayers();

        isSetCosine = true;
      }

      /// \brief Set production height in km of neutrinos
      /// \details Adds a layer of length heightKM with zero density to the density model
      /// @param heightKM Set neutrino production height
      virtual void setProductionHeight(FLOAT_T heightKM){
        if(!isSetCosine)
          throw std::runtime_error("must set cosine list before production height");

        ProductionHeightinCentimeter = heightKM * 100000.0;

        isSetProductionHeight = true;
      }

      virtual void setProductionHeightList(const std::vector<FLOAT_T>& list_prob, const std::vector<FLOAT_T>& list_bins) {
        if (!this->useProductionHeightAveraging) {
          throw std::runtime_error("Propagator::setProductionHeightList. Trying to set Production Height information but propagator is not expecting to use it");
        }

        if (list_prob.size() != this->nProductionHeightBins*2*3*n_energies*n_cosines) {
          throw std::runtime_error("Propagator::setProductionHeightList. Prob array is not the expected size");
        }

        if (list_bins.size()-1 != this->nProductionHeightBins) {
          throw std::runtime_error("Propagator::setProductionHeightList. ProductionHeightBins array is not expected size");
        }

        int MaxSize = Constants<FLOAT_T>::MaxProdHeightBins()*2*3*n_energies*n_cosines;
        productionHeightList_prob = std::vector<FLOAT_T>(MaxSize);
        for (unsigned int i=0;i<MaxSize;i++) {
          productionHeightList_prob[i] = 0.;
        }

        for (unsigned int i=0;i<list_prob.size();i++) {
          productionHeightList_prob[i] = list_prob[i];
        }

        productionHeightList_bins = std::vector<FLOAT_T>(Constants<FLOAT_T>::MaxProdHeightBins()+1);
        for (unsigned int i=0;i<Constants<FLOAT_T>::MaxProdHeightBins();i++) {
          productionHeightList_bins[i] = 0.;
        }

        for (unsigned int i=0;i<list_bins.size();i++) {
          productionHeightList_bins[i] = list_bins[i];
        }

        isSetProductionHeightArray = true;
      }

      /// \brief Set chemical composition of each layer in the Earth model
      /// \details Set chemical composition of each layer in the Earth model
      /// @param
      virtual void setChemicalComposition(const std::vector<FLOAT_T>& list) = 0;

      /// \brief Calculate the probability of each cell
      /// @param type Neutrino or Antineutrino
      virtual void calculateProbabilities(NeutrinoType type) = 0;

      /// \brief get oscillation weight for specific cosine and energy
      /// @param index_cosine Cosine bin index (zero based)
      /// @param index_energy Energy bin index (zero based)
      /// @param t Specify which probability P(i->j)
      virtual FLOAT_T getProbability(int index_cosine, int index_energy, ProbType t) = 0;

      /// \brief get oscillation weight
      /// @param probArr Cosine bin index (zero based)
      /// @param t Specify which probability P(i->j)
      virtual void getProbabilityArr(FLOAT_T* probArr, ProbType t) = 0;

    protected:
      // for each cosine bin, determine the number of layers which will be crossed by the neutrino path
      // the atmospheric layers is excluded
      virtual void setMaxlayers(){
        for(int index_cosine = 0; index_cosine < n_cosines; index_cosine++){
          FLOAT_T c = cosineList[index_cosine];
          const int maxLayer = std::count_if(coslimit.begin(), coslimit.end(), [c](FLOAT_T limit){ return c < limit;});

          if (maxLayer > Constants<FLOAT_T>::MaxNLayers()) {
            std::cerr << "Invalid number of maxLayer:" << maxLayer << std::endl;
            std::cerr << "Need to increase value of Constants<FLOAT_T>::MaxNLayers() in $CUDAPROB3/constants.hpp" << std::endl;
            throw std::runtime_error("setMaxlayers : invalid number of maxLayer");
          }

          maxlayers[index_cosine] = maxLayer;
        }
      }

      cudaprob3::math::ComplexNumber<FLOAT_T>& U(int i, int j){
        return Mix_U[( i * 3 + j)];
      }

      FLOAT_T& DM(int i, int j){
        return dm[( i * 3 + j)];
      }

      std::vector<FLOAT_T> energyList;
      std::vector<FLOAT_T> cosineList;
      std::vector<int> maxlayers;
      //std::vector<FLOAT_T> pathLengths;

      std::vector<FLOAT_T> productionHeightList_prob;
      std::vector<FLOAT_T> productionHeightList_bins;

      std::vector<FLOAT_T> radii;
      std::vector<FLOAT_T> rhos;
      std::vector<FLOAT_T> as;
      std::vector<FLOAT_T> bs;
      std::vector<FLOAT_T> cs;
      std::vector<FLOAT_T> yps;
      std::vector<FLOAT_T> coslimit;

      std::array<cudaprob3::math::ComplexNumber<FLOAT_T>, 9> Mix_U; // MNS mixing matrix
      std::array<FLOAT_T, 9> dm; // mass differences;

      FLOAT_T ProductionHeightinCentimeter;

      bool useProductionHeightAveraging = false;
      int nProductionHeightBins = 0;

      bool isSetProductionHeightArray = false;
      bool isSetProductionHeight = false;
      bool isSetCosine = false;
      bool isInit = true;

      int n_cosines;
      int n_energies;
      int n_layers;

      // Use polynomial density for density averaging each track?
      bool UsePolyDensity;
    };





} // namespace cudaprob3


#endif
