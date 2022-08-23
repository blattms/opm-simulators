/*
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 Statoil AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef OPM_MILU_HEADER_INCLUDED
#define OPM_MILU_HEADER_INCLUDED

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>
#include <functional>  // so we can use std::function


namespace Opm
{

enum class MILU_VARIANT{
    /// \brief Do not perform modified ILU
    ILU = 0,
    /// \brief \f$U_{ii} = U_{ii} +\f$  sum(dropped entries)
    MILU_1 = 1,
    /// \brief \f$U_{ii} = U_{ii} + sign(U_{ii}) * \f$ sum(dropped entries)
    MILU_2 = 2,
    /// \brief \f$U_{ii} = U_{ii} sign(U_{ii}) * \f$ sum(|dropped entries|)
    MILU_3 = 3,
    /// \brief \f$U_{ii} = U_{ii} + (U_{ii}>0?1:0) * \f$ sum(dropped entries)
    MILU_4 = 4
};

MILU_VARIANT convertString2Milu(const std::string& milu);


  namespace detail 
{
  /*  class AuxFunctors
{
public:
  static std::function<double(const double&)> identityFunctor(); 
  static std::function<double(const double&)> oneFunctor(); 
  static std::function<double(const double&)> signFunctor(); 
  static std::function<double(const double&)> isPositiveFunctor(); 
  static std::function<double(const double&)> absFunctor(); 
  };*/
  
  template <typename T>
  T identityFunctor(const T&);
  //std::function<T(const T&)> identityFunctor();

  template <typename T>
  T oneFunctor(const T&);
  // std::function<T(const T&)> oneFunctor();

  template <typename T>
  T signFunctor(const T&);
  // std::function<T(const T&)> signFunctor();

  template <typename T>
  T isPositiveFunctor(const T&);
  // std::function<T(const T&)> isPositiveFunctor();

  template <typename T>
  T absFunctor(const T&);
  // std::function<T(const T&)> absFunctor();

struct Reorderer
{
    virtual std::size_t operator[](std::size_t i) const = 0;
    virtual ~Reorderer() {}
};

struct NoReorderer : public Reorderer
{
    virtual std::size_t operator[](std::size_t i) const
    {
        return i;
    }
};

struct RealReorderer : public Reorderer
{
    RealReorderer(const std::vector<std::size_t>& ordering)
        : ordering_(&ordering)
    {}
    virtual std::size_t operator[](std::size_t i) const
    {
        return (*ordering_)[i];
    }
    const std::vector<std::size_t>* ordering_;
};

  
template<class M>
void milu0_decomposition(M& A, double absFunctor(const double&) = signFunctor<Opm::detail<double>>(const double&),
			 double signFunctor(const double&) = oneFunctor<Opm::detail<double>>(const double&),
                         std::vector<typename M::block_type>* diagonal = nullptr);

template<class M>
void  milu0_decomposition(M& A, std::vector<typename M::block_type>* diagonal)    
 {
   milu0_decomposition(A, identityFunctor<Opm::detail<double>>(const double&), oneFunctor<Opm::detail<double>>(const double&), diagonal);
 };


template<class M>
void milun_decomposition(const M& A, int n, MILU_VARIANT milu, M& ILU,
			 Reorderer& ordering, Reorderer& inverseOrdering);    

} // end namespace details

} // end namespace Opm

#endif
