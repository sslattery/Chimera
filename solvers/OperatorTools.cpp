//---------------------------------------------------------------------------//
// \file OperatorTools.cpp
// \author Stuart R. Slattery
// \brief OperatorTools definition.
//---------------------------------------------------------------------------//

#include "OperatorTools.hpp"

#include <Teuchos_ParameterList.hpp>

#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>

#include <AnasaziTypes.hpp>
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziEpetraAdapter.hpp"

namespace HMCSA
{

/*!
 * \brief Compute the spectral radius of an operator.
 */
double 
OperatorTools::spectralRadius( const Teuchos::RCP<Epetra_Operator>& matrix )
{
    Epetra_Map map = matrix->OperatorDomainMap();

    const int nev = 1;
    const int block_size = 1;
    const int num_blocks = 10;
    const int max_restarts = 100;
    const double tol = 1.0e-4;

    Teuchos::ParameterList krylovschur_params;
    krylovschur_params.set( "Which", "LM" );
    krylovschur_params.set( "Block Size", block_size );
    krylovschur_params.set( "Num Blocks", num_blocks );
    krylovschur_params.set( "Maximum Restarts", max_restarts );
    krylovschur_params.set( "Convergence Tolerance", tol );

    Teuchos::RCP<Epetra_MultiVector> vec 
	= Teuchos::rcp( new Epetra_MultiVector( map, block_size ) );
    vec->Random();

    Teuchos::RCP<Anasazi::BasicEigenproblem< double, 
					     Epetra_MultiVector, 
					     Epetra_Operator> > eigen_problem =
	Teuchos::rcp( 
	    new Anasazi::BasicEigenproblem< double, 
					    Epetra_MultiVector, 
					    Epetra_Operator >(matrix, vec) );
    eigen_problem->setNEV( nev );
    eigen_problem->setProblem();

    Anasazi::BlockKrylovSchurSolMgr< double, 
				     Epetra_MultiVector, 
				     Epetra_Operator> 
	solver_manager(eigen_problem, krylovschur_params);
    solver_manager.solve();

    Anasazi::Eigensolution<double,Epetra_MultiVector> sol = 
	eigen_problem->getSolution();
    std::vector<Anasazi::Value<double> > evals = sol.Evals;
    Teuchos::RCP<Epetra_MultiVector> evecs = sol.Evecs;

    double spectral_radius = pow( evals[0].realpart*evals[0].realpart +
				  evals[0].imagpart*evals[0].imagpart,
				  0.5 );

    return spectral_radius;
}

/*!
 * \brief Compute the stiffness ratio of an operator.
 */
double 
OperatorTools::stiffnessRatio( const Teuchos::RCP<Epetra_Operator>& matrix )
{
    return 0.0;
}

} // end namespace HMCSA

//---------------------------------------------------------------------------//
// end OperatorTools.cpp
//---------------------------------------------------------------------------//

