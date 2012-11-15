//---------------------------------------------------------------------------//
/*!
 * \file DiffusionProblem.cpp
 * \author Stuart R. Slattery
 * \brief Diffusion problem implementation.
 */
//---------------------------------------------------------------------------//

#include "DiffusionProblem.hpp"

#include <Chimera_Assertion.hpp>

#include <Teuchos_ArrayView.hpp>

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
    // Get the local rows from the partitioner.
    Teuchos::ArrayView<int> local_rows = partitioner->getLocalRows();

    // Build the map.
    Teuchos::RCP<const Tpetra::Map<int> > row_map = 
	    Tpetra::createNonContigMap<int,int>( local_rows, comm );

    // Build the operator.
    Teuchos::RCP<Tpetra::CrsMatrix<double,int,int> > A = 
	Tpetra::createCrsMatrix<double,int,int>( row_map );

    int N = partitioner->getGlobalEdges().first.size();
    double dx = partitioner->getCellSizes().first;
    double dy = partitioner->getCellSizes().second;
    double xs_s = plist->get<double>("SCATTERING XS");
    double xs_a = plist->get<double>("ABSORPTION XS");
    double D = 1.0 / ( 3.0*(xs_a+xs_s) );

    Teuchos::Array<double> diag( 1, xs_a + D*10.0/(3.0*dx*dx) );
    Teuchos::Array<double> iminus1( 1, -2.0*D/(3.0*dx*dx) );
    Teuchos::Array<double> iplus1( 1, -2.0*D/(3.0*dx*dx) );
    Teuchos::Array<double> jminus1( 1, -2.0*D/(3.0*dy*dy) );
    Teuchos::Array<double> jplus1( 1, -2.0*D/(3.0*dy*dy) );
    Teuchos::Array<double> iminus1jminus1( 1, -D/(6.0*dx*dx) );
    Teuchos::Array<double> iplus1jminus1( 1, -D/(6.0*dx*dx) );
    Teuchos::Array<double> iminus1jplus1( 1, -D/(6.0*dy*dy) );
    Teuchos::Array<double> iplus1jplus1( 1, -D/(6.0*dy*dy) );

    Teuchos::Array<int> idx(1);
    Teuchos::Array<int> idx_iminus1(1);
    Teuchos::Array<int> idx_iplus1(1);
    Teuchos::Array<int> idx_jminus1(1);
    Teuchos::Array<int> idx_jplus1(1);
    Teuchos::Array<int> idx_iminus1jminus1(1);
    Teuchos::Array<int> idx_iplus1jminus1(1);
    Teuchos::Array<int> idx_iminus1jplus1(1);
    Teuchos::Array<int> idx_iplus1jplus1(1);

    // Fill the operator from global data on proc 0.
    if ( comm->getRank() == 0 )
    {
	// Min X boundary vacuum (nonreentrant current).
	for ( int j = 1; j < N-1; ++j )
	{
	    int i = 0;

	    idx[0]                = i + j*N;
	    idx_iplus1[0]         = (i+1) + j*N;
	    idx_jminus1[0]        = i + (j-1)*N;
	    idx_jplus1[0]         = i + (j+1)*N;
	    idx_iplus1jminus1[0]  = (i+1) + (j-1)*N;
	    idx_iplus1jplus1[0]   = (i+1) + (j+1)*N;

	    A->insertGlobalValues( idx[0], idx(),                diag() );
	    A->insertGlobalValues( idx[0], idx_iplus1(),         iplus1  );
	    A->insertGlobalValues( idx[0], idx_jminus1(),        jminus1()  );
	    A->insertGlobalValues( idx[0], idx_jplus1(),         jplus1()  );
	    A->insertGlobalValues( idx[0], idx_iplus1jminus1(),  iplus1jminus1() );
	    A->insertGlobalValues( idx[0], idx_iplus1jplus1(),   iplus1jplus1() );
	}

	// Max X boundary vacuum (nonreentrant current).
	for ( int j = 1; j < N-1; ++j )
	{
	    int i = N-1;

	    idx[0]                = i + j*N;
	    idx_iminus1[0]        = (i-1) + j*N;
	    idx_jminus1[0]        = i + (j-1)*N;
	    idx_jplus1[0]         = i + (j+1)*N;
	    idx_iminus1jminus1[0] = (i-1) + (j-1)*N;
	    idx_iminus1jplus1[0]  = (i-1) + (j+1)*N;

	    A->insertGlobalValues( idx[0], idx(),                diag() );
	    A->insertGlobalValues( idx[0], idx_iminus1(),        iminus1()  );
	    A->insertGlobalValues( idx[0], idx_jminus1(),        jminus1()  );
	    A->insertGlobalValues( idx[0], idx_jplus1(),         jplus1()  );
	    A->insertGlobalValues( idx[0], idx_iminus1jminus1(), iminus1jminus1() );
	    A->insertGlobalValues( idx[0], idx_iminus1jplus1(),  iminus1jplus1() );
	}

	// Min Y boundary vacuum (nonreentrant current).
	for ( int i = 1; i < N-1; ++i )
	{
	    int j = 0;

	    idx[0]                = i + j*N;
	    idx_iminus1[0]        = (i-1) + j*N;
	    idx_iplus1[0]         = (i+1) + j*N;
	    idx_jplus1[0]         = i + (j+1)*N;
	    idx_iminus1jplus1[0]  = (i-1) + (j+1)*N;
	    idx_iplus1jplus1[0]   = (i+1) + (j+1)*N;

	    A->insertGlobalValues( idx[0], idx(),                diag() );
	    A->insertGlobalValues( idx[0], idx_iminus1(),        iminus1()  );
	    A->insertGlobalValues( idx[0], idx_iplus1(),         iplus1  );
	    A->insertGlobalValues( idx[0], idx_jplus1(),         jplus1()  );
	    A->insertGlobalValues( idx[0], idx_iminus1jplus1(),  iminus1jplus1() );
	    A->insertGlobalValues( idx[0], idx_iplus1jplus1(),   iplus1jplus1() );
	}

	// Max Y boundary vacuum (nonreentrant current).
	for ( int i = 1; i < N-1; ++i )
	{
	    int j = N-1;

	    idx[0]                = i + j*N;
	    idx_iminus1[0]        = (i-1) + j*N;
	    idx_iplus1[0]         = (i+1) + j*N;
	    idx_jminus1[0]        = i + (j-1)*N;
	    idx_iminus1jminus1[0] = (i-1) + (j-1)*N;
	    idx_iplus1jminus1[0]  = (i+1) + (j-1)*N;

	    A->insertGlobalValues( idx[0], idx(),                diag() );
	    A->insertGlobalValues( idx[0], idx_iminus1(),        iminus1()  );
	    A->insertGlobalValues( idx[0], idx_jminus1(),        jminus1()  );
	    A->insertGlobalValues( idx[0], idx_jplus1(),         jplus1()  );
	    A->insertGlobalValues( idx[0], idx_iminus1jminus1(), iminus1jminus1() );
	    A->insertGlobalValues( idx[0], idx_iplus1jminus1(),  iplus1jminus1() );
	}

	// Lower left boundary vacuum (nonreentrant current).
	{
	    int i = 0;
	    int j = 0;

	    idx[0]                = i + j*N;
	    idx_iplus1[0]         = (i+1) + j*N;
	    idx_jplus1[0]         = i + (j+1)*N;
	    idx_iplus1jplus1[0]   = (i+1) + (j+1)*N;

	    A->insertGlobalValues( idx[0], idx(),                diag() );
	    A->insertGlobalValues( idx[0], idx_iplus1(),         iplus1  );
	    A->insertGlobalValues( idx[0], idx_jplus1(),         jplus1()  );
	    A->insertGlobalValues( idx[0], idx_iplus1jplus1(),   iplus1jplus1() );
	}

	// Lower right boundary vacuum (nonreentrant current).
	{
	    int i = N-1;
	    int j = 0;

	    idx[0]                = i + j*N;
	    idx_iminus1[0]        = (i-1) + j*N;
	    idx_jplus1[0]         = i + (j+1)*N;
	    idx_iminus1jplus1[0]  = (i-1) + (j+1)*N;

	    A->insertGlobalValues( idx[0], idx(),                diag() );
	    A->insertGlobalValues( idx[0], idx_iminus1(),        iminus1()  );
	    A->insertGlobalValues( idx[0], idx_jplus1(),         jplus1()  );
	    A->insertGlobalValues( idx[0], idx_iminus1jplus1(),  iminus1jplus1() );
	}

	// Upper left boundary vacuum (nonreentrant current).
	{
	    int i = 0;
	    int j = N-1;

	    idx[0]                = i + j*N;
	    idx_iplus1[0]         = (i+1) + j*N;
	    idx_jminus1[0]        = i + (j-1)*N;
	    idx_iplus1jminus1[0]  = (i+1) + (j-1)*N;

	    A->insertGlobalValues( idx[0], idx(),                diag() );
	    A->insertGlobalValues( idx[0], idx_iplus1(),         iplus1  );
	    A->insertGlobalValues( idx[0], idx_jminus1(),        jminus1()  );
	    A->insertGlobalValues( idx[0], idx_iplus1jminus1(),  iplus1jminus1() );
	}

	// Upper right boundary vacuum (nonreentrant current).
	{
	    int i = N-1;
	    int j = N-1;

	    idx[0]                = i + j*N;
	    idx_iminus1[0]        = (i-1) + j*N;
	    idx_jminus1[0]        = i + (j-1)*N;
	    idx_iminus1jminus1[0] = (i-1) + (j-1)*N;

	    A->insertGlobalValues( idx[0], idx(),                diag() );
	    A->insertGlobalValues( idx[0], idx_iminus1(),        iminus1()  );
	    A->insertGlobalValues( idx[0], idx_jminus1(),        jminus1()  );
	    A->insertGlobalValues( idx[0], idx_iminus1jminus1(), iminus1jminus1() );
	}

	// Central grid points
	for ( int i = 1; i < N-1; ++i )
	{
	    for ( int j = 1; j < N-1; ++j )
	    {
		idx[0]                = i + j*N;
		idx_iminus1[0]        = (i-1) + j*N;
		idx_iplus1[0]         = (i+1) + j*N;
		idx_jminus1[0]        = i + (j-1)*N;
		idx_jplus1[0]         = i + (j+1)*N;
		idx_iminus1jminus1[0] = (i-1) + (j-1)*N;
		idx_iplus1jminus1[0]  = (i+1) + (j-1)*N;
		idx_iminus1jplus1[0]  = (i-1) + (j+1)*N;
		idx_iplus1jplus1[0]   = (i+1) + (j+1)*N;

		A->insertGlobalValues( idx[0], idx(),                diag() );
		A->insertGlobalValues( idx[0], idx_iminus1(),        iminus1()  );
		A->insertGlobalValues( idx[0], idx_iplus1(),         iplus1  );
		A->insertGlobalValues( idx[0], idx_jminus1(),        jminus1()  );
		A->insertGlobalValues( idx[0], idx_jplus1(),         jplus1()  );
		A->insertGlobalValues( idx[0], idx_iminus1jminus1(), iminus1jminus1() );
		A->insertGlobalValues( idx[0], idx_iplus1jminus1(),  iplus1jminus1() );
		A->insertGlobalValues( idx[0], idx_iminus1jplus1(),  iminus1jplus1() );
		A->insertGlobalValues( idx[0], idx_iplus1jplus1(),   iplus1jplus1() );
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

    // Build the source.
    double source_strength = plist->get<double>("SOURCE STRENGTH");
    Teuchos::RCP<Tpetra::Vector<double,int> > B = 
	Tpetra::createVector<double,int>( row_map );

    if ( plist->get<std::string>("SOURCE TYPE") == "UNIFORM" )
    {
	B->putScalar( source_strength );
    }
    else if ( plist->get<std::string>("SOURCE TYPE") == "POINT" )
    {
	int source_location = plist->get<int>("SOURCE LOCATION");
	if ( row_map->isNodeGlobalElement( source_location ) )
	{
	    B->replaceGlobalValue( source_location, source_strength );
	}
    }

    // Build the linear problem.
    comm->barrier();
    d_linear_problem = Teuchos::rcp( new LinearProblemType( A, X, B ) );
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

