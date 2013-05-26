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
#include <Chimera_OperatorTools.hpp>

#include <MCLS_SolverFactory.hpp>
#include <MCLS_SolverManager.hpp>
#include <MCLS_LinearProblem.hpp>
#include <MCLS_TpetraAdapter.hpp>

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
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Operator.hpp>

#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>

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

    // Read in command line options.
    std::string xml_input_filename;
    Teuchos::CommandLineProcessor clp(false);
    clp.setOption( "xml-in-file",
		   &xml_input_filename,
		   "The XML file to read into a parameter list" );
    clp.parse(argc,argv);

    // Build the parameter list from the xml input.
    Teuchos::RCP<Teuchos::ParameterList> plist =
	Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile(
	xml_input_filename, Teuchos::inoutArg(*plist) );

    // Build and partition the mesh.
    Teuchos::RCP<Chimera::Partitioner> partitioner = Teuchos::rcp(
	new Chimera::Partitioner( comm, plist ) );

    // Build the diffusion problem.
    Teuchos::RCP<Chimera::DiffusionProblem> diffusion_problem = Teuchos::rcp(
	new Chimera::DiffusionProblem( comm, partitioner, plist, true ) );

    // // CHIMERA SOLVE
    // Build the solver.
    Teuchos::RCP<Chimera::LinearSolver<double,int,int> > csolver =
        Chimera::LinearSolverFactory::create( 
            plist, diffusion_problem->getProblem() );

    // Compute the iteration matrix spectral radius.
    double spec_rad = Chimera::OperatorTools::spectralRadius(
        csolver->linearOperatorSplit()->iterationMatrix() );
    if ( comm->getRank() == 0 )
    {
        std::cout << "SPECTRAL RADIUS: " << spec_rad << std::endl; 
    }
    comm->barrier();

    typedef Tpetra::Vector<double,int,int> Vector;
    typedef Tpetra::CrsMatrix<double,int,int> Matrix;
    typedef Tpetra::MultiVector<double,int,int> MV;
    typedef Tpetra::Operator<double,int,int> OP;

    // MCLS SOLVE
    if ( "MCLS" == plist->get<std::string>("Solver Package") )
    {
        Teuchos::RCP<MCLS::LinearProblem<Vector,Matrix> > linear_problem =
            Teuchos::rcp( new MCLS::LinearProblem<Vector,Matrix>(
                              diffusion_problem->getProblem()->getOperator(),
                              diffusion_problem->getProblem()->getLHS(),
                              diffusion_problem->getProblem()->getRHS() ) );

        std::string solver_type = plist->get<std::string>("Solver Type");

        MCLS::SolverFactory<Vector,Matrix> factory;
        Teuchos::RCP<MCLS::SolverManager<Vector,Matrix> > solver_manager =
            factory.create( solver_type, comm, plist );
        solver_manager->setProblem( linear_problem );
        Teuchos::Time timer("");
        timer.start(true);
        solver_manager->solve();
        timer.stop();
        if ( comm->getRank() == 0 )
        {
            std::cout << "MCLS Solve: Complete in " << timer.totalElapsedTime() 
                      << " seconds." << std::endl;
        }
    }

    // Belos Solve
    else if ( "Belos" == plist->get<std::string>("Solver Package") )
    {
        Teuchos::RCP<Teuchos::ParameterList> belosList = Teuchos::parameterList();
        belosList->set( "Maximum Iterations", 1000 );
        belosList->set( "Maximum Restarts", 10 );
        belosList->set( "Convergence Tolerance", 1.0e-8 );
        int verbosity = Belos::Errors + Belos::Warnings + 
                        Belos::TimingDetails + Belos::StatusTestDetails;
        belosList->set( "Verbosity", verbosity );

        Teuchos::RCP<Belos::LinearProblem<double,MV,OP> > linear_problem 
            = Teuchos::rcp( new Belos::LinearProblem<double,MV,OP>(
                                diffusion_problem->getProblem()->getOperator(),
                                diffusion_problem->getProblem()->getLHS(),
                                diffusion_problem->getProblem()->getRHS() ) );

        Teuchos::RCP<Belos::SolverManager<double,MV,OP> > solver_manager;
        
        if ( "GMRES" == plist->get<std::string>("Solver Type") )
        {
            solver_manager = Teuchos::rcp( 
                new Belos::PseudoBlockGmresSolMgr<double,MV,OP>(
                    linear_problem, belosList) );
        }
        else if ( "CG" == plist->get<std::string>("Solver Type") )
        {
            solver_manager = Teuchos::rcp( 
                new Belos::PseudoBlockCGSolMgr<double,MV,OP>(
                    linear_problem, belosList) );
        }

        Teuchos::Time timer("");
        timer.start(true);
        solver_manager->solve();
        timer.stop();
        if ( comm->getRank() == 0 )
        {
            std::cout << "Belos Solve: Complete in " << timer.totalElapsedTime() 
                      << " seconds." << std::endl;
            int numIters = solver_manager->getNumIters();
            std::cout << "Number of iterations performed for this solve: " 
                      << numIters << std::endl;
        }
    }

    // Write the solution to VTK.
    Chimera::VTKOutput vtk_output( comm, partitioner, plist );
    vtk_output.addField( Chimera::VTKOutput::VERTEX_FIELD,
                         diffusion_problem->getProblem()->getLHS(),
                         "NEUTRON_FLUX" );
    vtk_output.write();

    return 0;
}

//---------------------------------------------------------------------------//
// end neutron_diffusion.cpp
//---------------------------------------------------------------------------//

