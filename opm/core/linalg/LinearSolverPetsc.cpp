/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2014 STATOIL ASA.

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

#include "config.h"
#include <cstring>
#include <opm/core/linalg/LinearSolverPetsc.hpp>
#include <unordered_map>
#define PETSC_CLANGUAGE_CXX 1 //enable CHKERRXX macro.
#include <petsc.h>
#include <opm/core/utility/ErrorMacros.hpp>

namespace Opm
{
    
namespace{

    class KSPTypeMap {
    public:
        explicit
        KSPTypeMap(const std::string& default_type = "gmres")
            : default_type_(default_type)
        {
            type_map_.insert(std::make_pair("richardson", KSPRICHARDSON));
            type_map_.insert(std::make_pair("chebyshev", KSPCHEBYSHEV));
            type_map_.insert(std::make_pair("cg", KSPCG));
            type_map_.insert(std::make_pair("bicgs", KSPBICG));
            type_map_.insert(std::make_pair("gmres", KSPGMRES));
            type_map_.insert(std::make_pair("fgmres", KSPFGMRES));
            type_map_.insert(std::make_pair("dgmres", KSPDGMRES));
            type_map_.insert(std::make_pair("gcr", KSPGCR));
            type_map_.insert(std::make_pair("bcgs", KSPBCGS));
            type_map_.insert(std::make_pair("cgs", KSPCGS));
            type_map_.insert(std::make_pair("tfqmr", KSPTFQMR));
            type_map_.insert(std::make_pair("tcqmr", KSPTCQMR));
            type_map_.insert(std::make_pair("cr", KSPCR));
            type_map_.insert(std::make_pair("preonly", KSPPREONLY));
        }
          
        KSPType 
        find(const std::string& type) const
        {
            Map::const_iterator it = type_map_.find(type);

            if (it == type_map_.end()) {
                it = type_map_.find(default_type_);
            }
        
            if (it == type_map_.end()) {
                OPM_THROW(std::runtime_error, "Unknown KSPType: '" << type << "'");
            }
            return it->second;
        }
    private:
        typedef std::unordered_map<std::string, KSPType> Map;

        std::string default_type_;
        Map type_map_;
    };


    class PCTypeMap {
    public:
        explicit
        PCTypeMap(const std::string& default_type = "jacobi")
            : default_type_(default_type)
        {
            type_map_.insert(std::make_pair("jacobi", PCJACOBI));
            type_map_.insert(std::make_pair("bjacobi", PCBJACOBI));
            type_map_.insert(std::make_pair("sor", PCSOR));
            type_map_.insert(std::make_pair("eisenstat", PCEISENSTAT));
            type_map_.insert(std::make_pair("icc", PCICC));
            type_map_.insert(std::make_pair("ilu", PCILU));
            type_map_.insert(std::make_pair("asm", PCASM));
            type_map_.insert(std::make_pair("gamg", PCGAMG));
            type_map_.insert(std::make_pair("ksp", PCKSP));
            type_map_.insert(std::make_pair("composite", PCCOMPOSITE));
            type_map_.insert(std::make_pair("lu", PCLU));
            type_map_.insert(std::make_pair("cholesky", PCCHOLESKY));
            type_map_.insert(std::make_pair("none", PCNONE));
        }
            
        PCType 
        find(const std::string& type) const
        {
            Map::const_iterator it = type_map_.find(type);

            if (it == type_map_.end()) {
                it = type_map_.find(default_type_);
            }
       
            if (it == type_map_.end()) {
                OPM_THROW(std::runtime_error, "Unknown PCType: '" << type << "'");
            }
            return it->second;
        }
    private:
        typedef std::unordered_map<std::string, PCType> Map;

        std::string default_type_;
        Map type_map_;
    };

    struct OEM_DATA {
        /* Convenience struct to handle automatic (de)allocation of some useful
         * variables, as well as group them up for easier parameter passing
         */
        Vec rhs;
        Vec solution;
        Mat A;
        KSP ksp;
        PC preconditioner;

        OEM_DATA( const int size ) {
            CHKERRXX( VecCreate( PETSC_COMM_WORLD, &solution ) );
            CHKERRXX( VecSetSizes( solution, PETSC_DECIDE, size ) );
            CHKERRXX( VecSetFromOptions( solution ) );
            CHKERRXX( VecDuplicate( solution, &rhs ) );

            CHKERRXX( MatCreate( PETSC_COMM_WORLD, &A ) );
            CHKERRXX( MatSetSizes( A, PETSC_DECIDE, PETSC_DECIDE, size, size ) );
            CHKERRXX( MatSetFromOptions( A ) );
            CHKERRXX( MatSetUp( A ) );
        }

        ~OEM_DATA() {
            CHKERRXX( VecDestroy( &rhs ) );
            CHKERRXX( VecDestroy( &solution ) );
            CHKERRXX( MatDestroy( &A ) );
            CHKERRXX( KSPDestroy( &ksp ) );
        }
    };

    void to_petsc_vec( const double* x, Vec v ) {
        if( !v ) OPM_THROW( std::runtime_error,
                    "PETSc CopySolution: Invalid PETSc vector." );

        PetscScalar* vec;
        PetscInt size;

        CHKERRXX( VecGetLocalSize( v, &size ) );
        CHKERRXX( VecGetArray( v, &vec ) );

        std::memcpy( vec, x,  size * sizeof( double ) );
        CHKERRXX( VecRestoreArray( v, &vec ) );
    }

    void from_petsc_vec( double* x, Vec v ) {
        if( !v ) OPM_THROW( std::runtime_error,
                    "PETSc CopySolution: Invalid PETSc vector." );

        PetscScalar* vec;
        PetscInt size;

        CHKERRXX( VecGetLocalSize( v, &size ) );
        CHKERRXX( VecGetArray( v, &vec ) );

        std::memcpy( x, vec, size * sizeof( double ) );
        CHKERRXX( VecRestoreArray( v, &vec ) );
    }

    void to_petsc_mat( const int size, const int nonzeros,
            const int* ia, const int* ja, const double* sa, Mat A ) {

        for( int i = 0; i < size; ++i ) {
            int nrows = ia[ i + 1 ] - ia[ i ];

            if( nrows == 0 ) continue;

            for( int j = ia[ i ]; j < ia[ i + 1 ]; ++j )
                CHKERRXX( MatSetValues( A, 1, &i, 1, &ja[ j ], &sa[ j ], INSERT_VALUES ) );

            CHKERRXX( MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY ) );
            CHKERRXX( MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY ) );
        }
    }


    void solve_system( OEM_DATA& t, KSPType method, PCType pcname,
            double rtol, double atol, double dtol, int maxits, int ksp_view ) {
        PetscInt its;
        PetscReal residual;
        KSPConvergedReason reason;

        CHKERRXX( KSPCreate( PETSC_COMM_WORLD, &t.ksp ) );
        CHKERRXX( KSPSetOperators( t.ksp, t.A, t.A, DIFFERENT_NONZERO_PATTERN ) );
        CHKERRXX( KSPGetPC( t.ksp, &t.preconditioner ) );
        CHKERRXX( KSPSetType( t.ksp, method ) );
        CHKERRXX( PCSetType( t.preconditioner, pcname ) );
        CHKERRXX( KSPSetTolerances( t.ksp, rtol, atol, dtol, maxits ) );
        CHKERRXX( KSPSetFromOptions( t.ksp ) );
        CHKERRXX( KSPSetInitialGuessNonzero( t.ksp, PETSC_TRUE ) );
        CHKERRXX( KSPSolve( t.ksp, t.rhs, t.solution ) );
        CHKERRXX( KSPGetConvergedReason( t.ksp, &reason ) );
        CHKERRXX( KSPGetIterationNumber( t.ksp, &its ) );
        CHKERRXX( KSPGetResidualNorm( t.ksp, &residual ) );

        if( ksp_view )
            CHKERRXX( KSPView( t.ksp, PETSC_VIEWER_STDOUT_WORLD ) );

        CHKERRXX( PetscPrintf( PETSC_COMM_WORLD, "KSP Iterations %D, Final Residual %G\n", its, residual ) );
    }

} // anonymous namespace.

    LinearSolverPetsc::LinearSolverPetsc(const parameter::ParameterGroup& param)
        : ksp_type_( param.getDefault( std::string( "ksp_type" ), std::string( "gmres" ) ) )
        , pc_type_( param.getDefault( std::string( "pc_type" ), std::string( "sor" ) ) )
        , ksp_view_( param.getDefault( std::string( "ksp_view" ), int( false ) ) )
        , rtol_( param.getDefault( std::string( "ksp_rtol" ), 1e-5 ) )
        , atol_( param.getDefault( std::string( "ksp_atol" ), 1e-50 ) )
        , dtol_( param.getDefault( std::string( "ksp_dtol" ), 1e5 ) )
        , maxits_( param.getDefault( std::string( "ksp_max_it" ), 1e5 ) )
    {
        int argc = 0;
        char** argv = NULL;
        PetscInitialize(&argc, &argv, (char*)0, "Petsc interface for OPM!\n");
    }


    LinearSolverPetsc::~LinearSolverPetsc()
    {
       PetscFinalize();
    }


    LinearSolverInterface::LinearSolverReport
    LinearSolverPetsc::solve(const int size,
                               const int nonzeros,
                               const int* ia,
                               const int* ja,
                               const double* sa,
                               const double* rhs,
                               double* solution,
                               const boost::any&) const
    {
        KSPTypeMap ksp(ksp_type_);
        KSPType ksp_type = ksp.find(ksp_type_);
        PCTypeMap pc(pc_type_);
        PCType pc_type = pc.find(pc_type_);

        OEM_DATA t( size );
        to_petsc_mat( size, nonzeros, ia, ja, sa, t.A );
        to_petsc_vec( rhs, t.rhs );

        solve_system( t, ksp_type, pc_type, rtol_, atol_, dtol_, maxits_, ksp_view_ );
        from_petsc_vec( solution, t.solution );

        LinearSolverReport rep = {};
        rep.converged = true;
        return rep;
    }

    void LinearSolverPetsc::setTolerance(const double /*tol*/)
    {
    }

    double LinearSolverPetsc::getTolerance() const
    {
        return -1.;
    }

} // namespace Opm

