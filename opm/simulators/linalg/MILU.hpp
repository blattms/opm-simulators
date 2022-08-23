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



class auxiliaryFunctors
{
public:
  static std::function<double(const double&)> IdentityFunctor(); 
  static std::function<double(const double&)> OneFunctor(); 
  static std::function<double(const double&)> SignFunctor(); 
  static std::function<double(const double&)> IsPositiveFunctor(); 
  static std::function<double(const double&)> AbsFunctor(); 
};

//std::function<double(const double&)>
static std::function<double(const double&)> auxiliaryFunctors::IdentityFunctor(const double& t){return t;};
//std::function<double(const double&)> 
static std::function<double(const double&)> auxiliaryFunctors::OneFunctor(const double&){return 1.0;};
//std::function<double(const double&)>
static std::function<double(const double&)> auxiliaryFunctors::SignFunctor(const double& t){if (t< 0){return -1;} else{return 1;}};
//std::function<double(const double&)> 
static std::function<double(const double&)> auxiliaryFunctors::IsPositiveFunctor(const double& t){if (t<0){return 0;} else{return 1;}};
//std::function<double(const double&)>
static std::function<double(const double&)> auxiliaryFunctors::AbsFunctor(const double& t){return std::abs(t);};

 

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


  namespace detail  // : public auxiliaryFunctors
{

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

  ///////////// Anto-modifications start here ////////////////////////
  /*
struct IdentityFunctor
{
    template<class T>
    auto operator() = [](const T& t){return t;};  // when call it, operator()(<variable type const T& t>)?
  T operator()(const T& t)
    {
        return t;
	}
};

struct OneFunctor
{
    template<class T>
    auto operator()= [](const T&){return 1.0;};
  T operator()(const T&)
    {
        return 1.0;
	}
};
struct SignFunctor
{
    template<class T>
    auto operator() = [](const T& t){if (t< 0){return -1;} else{return 1;}};
     double operator()(const T& t)
    {
        if (t < 0.0)
        {
            return -1;
        }
        else
        {
            return 1.0;
        }
	}
    
};

struct IsPositiveFunctor
{
    template<class T>
    auto operator() = [](const T& t){if (t<0){return 0;} else{return 1;}};
    double operator()(const T& t)
     {
        if (t < 0.0)
        {
            return 0;
        }
        else
        {
            return 1;
        }
	}
};
struct AbsFunctor
{
    template<class T>
    auto operator() = [](const T& t){return std::abs(t);};
    T operator()(const T& t)
    {
        return std::abs(t);
	}
};
  
*/
  
  //template<class M, class F1=IdentityFunctor, class F2=OneFunctor >
 template<class M>
//void milu0_decomposition(M& A, F1 absFunctor = F1(), F2 signFunctor = F2(),
 void milu0_decomposition(M& A, std::function<double(double)> absFunctor = auxiliaryFunctors::SignFunctor(),
			  std::function<double(double)> signFunctor = auxiliaryFunctors::OneFunctor(),
                         std::vector<typename M::block_type>* diagonal = nullptr);

template<class M>
void  milu0_decomposition(M& A, std::vector<typename M::block_type>* diagonal)    
 {
   milu0_decomposition(A, auxiliaryFunctors::IdentityFunctor(), auxiliaryFunctors::OneFunctor(), diagonal);
 };
/*
{ milu0_decomposition(A, detail::IdentityFunctor(), detail::OneFunctor(),
                        diagonal); }*/

template<class M>
void milun_decomposition(const M& A, int n, MILU_VARIANT milu, M& ILU,
			 Reorderer& ordering, Reorderer& inverseOrdering);     // no changes, it's just a declaration

} // end namespace details

} // end namespace Opm

#endif
