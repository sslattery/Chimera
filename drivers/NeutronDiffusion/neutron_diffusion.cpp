//---------------------------------------------------------------------------//
/*!
 * \file neutron_diffusion.cpp
 * \author Stuart R. Slattery
 * \brief Neutron diffusion driver.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <string>
#include <cassert>

#include "DiffusionProblem.hpp"
#include "Partitioner.hpp"
#include "VTKOutput.hpp"

#include <Chimera_LinearSolver.hpp>
#include <Chimera_LinearSolverFactory.hpp>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>

//---------------------------------------------------------------------------//
int main( int argc, char * argv[] )
{
    // Initialize parallel communication.
    Teuchos::GlobalMPISession mpi_session( &argc, &argv );
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();
    Teuchos::FancyOStream out( Teuchos::rcpFromRef( std::cout ) );
    out.setOutputToRootOnly( 0 );
    out.setShowProcRank( true );

    // Build the parameter list.
    Teuchos::RCP<Teuchos::ParameterList> plist =
	Teuchos::rcp( new Teuchos::ParameterList() );

    // Mesh.
    plist->set<int>(         "I_BLOCKS",            1         );
    plist->set<int>(         "J_BLOCKS",            1         );
    plist->set<double>(      "X_MIN",               0.0       );
    plist->set<double>(      "X_MAX",               1.0       );
    plist->set<double>(      "Y_MIN",               0.0       );
    plist->set<double>(      "Y_MAX",               1.0       );
    plist->set<int>(         "X_NUM_CELLS",         10        );
    plist->set<int>(         "Y_NUM_CELLS",         10        );
    plist->set<std::string>( "GRID_TYPE",           "UNIFORM" );

    // Physics.
    plist->set<double>(      "ABSORPTION XS"        0.5       );
    plist->set<double>(      "SCATTERING XS"        0.5       );

    // Solver.
    plist->set<std::string>( "SOLVER TYPE",         "MCSA"    );
    plist->set<std::string>( "RNG TYPE",            "MT19937" );
    plist->set<double>(      "TOLERANCE",           1.0e-8    );
    plist->set<int>(         "MAX ITERS",           100       );
    plist->set<std::string>( "SPLIT TYPE",          "JACOBI"  );
    plist->set<std::string>( "MC TYPE",             "ADJOINT" );
    plist->set<double>(      "WEIGHT CUTOFF",       1.0e-4    );
    plist->set<int>(         "HISTORIES PER STAGE", 100       );
    plist->set<int>(         "NUM OVERLAP",         0         );

    // Output.
    plist->set<std::string>( "OUTPUT_FILENAME", "diffusion_out.vtk" );

    // Build and partition the mesh.
    Teuchos::RCP<Chimera::Partitioner> partitioner = Teuchos::rcp(
	new Chimera::Partitioner( comm, plist ) );

    // Build the diffusion problem.
    Teuchos::RCP<Chimera::DiffusionProblem> diffusion_problem = Teuchos::rcp(
	new Chimera::DiffusionProblem( comm, partitioner, plist ) );

    // Build the solver.
    Teuchos::RCP<Chimera::LinearSolver<double,int,int> > solver =
	Chimera::LinearSolverFactory::create( 
	    plist, diffusion_problem->getProblem() );

    // Solve.
    solver->iterate();

    // Write the solution to VTK.
    Chimera::VTKOutput vtk_output( comm, partitioner, plist );
    vtk_output.addField( Chimera::VERTEX_FIELD,
			 diffusion_problem->getProblem()->getLHS(),
			 "NEUTRON_FLUX" );
    vtk_output.write();

    return 0;
}

//---------------------------------------------------------------------------//
// end neutron_diffusion.cpp
//---------------------------------------------------------------------------//

