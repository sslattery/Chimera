//---------------------------------------------------------------------------//
/*!
 * \file DiffusionProblem.cpp
 * \author Stuart R. Slattery
 * \brief Diffusion problem implementation.
 */
//---------------------------------------------------------------------------//

#include "DiffusionProblem.hpp"

#include <Teuchos_Array.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>

namespace Chimera
{

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
DiffusionProblem::DiffusionProblem( const RCP_Comm& comm, 
				    const RCP_Partitioner& partitioner,
				    const RCP_ParameterList& plist )
{
    // Build the linear operator - this is a 2D Transient Diffusion operator.
    int N = 10;
    int problem_size = N*N;
    double dx = 0.01;
    double dy = 0.01;
    double dt = 0.05;
    double alpha = 0.001;
    Teuchos::RCP<const Tpetra::Map<int> > row_map = 
	Tpetra::createUniformContigMap<int,int>( problem_size, comm );
    Teuchos::RCP<Tpetra::CrsMatrix<double,int,int> > A = 
	Tpetra::createCrsMatrix<double,int,int>( row_map, 3 );

    Teuchos::Array<double> diag( 1, 1.0 + 2*dt*alpha*( 1/(dx*dx) + 1/(dy*dy) ) );
    Teuchos::Array<double> i_minus( 1, -dt*alpha/(dx*dx) );
    Teuchos::Array<double> i_plus( 1, -dt*alpha/(dx*dx) );
    Teuchos::Array<double> j_minus( 1, -dt*alpha/(dy*dy) );
    Teuchos::Array<double> j_plus( 1, -dt*alpha/(dy*dy) );

    Teuchos::Array<int> idx(1);
    Teuchos::Array<int> idx_iminus(1);
    Teuchos::Array<int> idx_iplus(1);
    Teuchos::Array<int> idx_jminus(1);
    Teuchos::Array<int> idx_jplus(1);
    Teuchos::Array<double> one(1,1.0);

    if ( comm->getRank() == 0 )
    {
	// Min X boundary Dirichlet.
	for ( int j = 1; j < N-1; ++j )
	{
	    int i = 0;
	    idx[0] = i + j*N;
	    A->insertGlobalValues( idx[0], idx(), one() );
	}

	// Max X boundary Dirichlet.
	for ( int j = 1; j < N-1; ++j )
	{
	    int i = N-1;
	    idx[0] = i + j*N;
	    A->insertGlobalValues( idx[0], idx(), one() );
	}

	// Min Y boundary Dirichlet.
	for ( int i = 0; i < N; ++i )
	{
	    int j = 0;
	    idx[0] = i + j*N;
	    A->insertGlobalValues( idx[0], idx(), one() );
	}

	// Max Y boundary Dirichlet.
	for ( int i = 0; i < N; ++i )
	{
	    int j = N-1;
	    idx[0] = i + j*N;
	    A->insertGlobalValues( idx[0], idx(), one() );
	}

	// Central grid points.
	for ( int i = 1; i < N-1; ++i )
	{
	    for ( int j = 1; j < N-1; ++j )
	    {
		idx[0] = i + j*N;
		idx_iminus[0] = (i-1) + j*N;
		idx_iplus[0] = (i+1) + j*N;
		idx_jminus[0] = i + (j-1)*N;
		idx_jplus[0] = i + (j+1)*N;
	    
		A->insertGlobalValues( idx[0], idx_jminus(), j_minus() );
		A->insertGlobalValues( idx[0], idx_iminus(), i_minus() );
		A->insertGlobalValues( idx[0], idx(),        diag()    );
		A->insertGlobalValues( idx[0], idx_iplus(),  i_plus()  );
		A->insertGlobalValues( idx[0], idx_jplus(),  j_plus()  );
	    }
	}
    }
    comm->barrier();

    A->fillComplete();

    // Build the solution vector.
    double X_val = 0.0;
    Teuchos::RCP<Tpetra::Vector<double,int> > X = 
	Tpetra::createVector<double,int>( row_map );
    X->putScalar( X_val );

    // Build the right hand side.
    double B_val = 0.0;
    Teuchos::RCP<Tpetra::Vector<double,int> > B = 
	Tpetra::createVector<double,int>( row_map );
    B->putScalar( B_val );

    // Set dirichlet boundary conditions.
    int row;
    // left
    for ( int j = 1; j < N-1; ++j )
    {
	int i = 0;
	row = i + j*N;
	if ( row_map->isNodeGlobalElement( row ) )
	{
	    B->replaceGlobalValue( row, 5.0 );
	}
    }
    // right
    for ( int j = 1; j < N-1; ++j )
    {
	int i = N-1;
	row = i + j*N;
	if ( row_map->isNodeGlobalElement( row ) )
	{
	    B->replaceGlobalValue( row, 5.0 );
	}
    }
    // bottom
    for ( int i = 0; i < N; ++i )
    {
	int j = 0;
	row = i + j*N;
	if ( row_map->isNodeGlobalElement( row ) )
	{
	    B->replaceGlobalValue( row, 0.0 );
	}
    }
    // top
    for ( int i = 0; i < N; ++i )
    {
	int j = N-1;
	row = i + j*N;
	if ( row_map->isNodeGlobalElement( row ) )
	{
	    B->replaceGlobalValue( row, 0.0 );
	}
    }

    // Build the linear problem.
    comm->barrier();
    d_linear_problem = Teuchos::rcp( new RCP_LinearProblem( A, X, B ) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
DiffusionProblem::~DiffusionProblem()
{ /* ... */ }

//---------------------------------------------------------------------------//

} // end namespace Chimera

//---------------------------------------------------------------------------//
// end DiffusionProblem.cpp
//---------------------------------------------------------------------------//

