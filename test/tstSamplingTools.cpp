//---------------------------------------------------------------------------//
/*!
 * \file tstSamplingTools.cpp
 * \author Stuart R. Slattery
 * \brief Boost random number generator unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <Chimera_RNGTraits.hpp>
#include <Chimera_BoostRNG.hpp>
#include <Chimera_SamplingTools.hpp>

#include <boost/random/mersenne_twister.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_as.hpp>

#include <Tpetra_Vector.hpp>

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( SamplingTools, uniform_stratify_sample_test )
{
    using namespace Chimera;

    // Setup parallel distribution.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();
    int comm_size = comm->getSize();

    // Setup linear operator distribution.
    int local_num_rows = 10;
    int global_num_rows = local_num_rows*comm_size;
    Teuchos::RCP<const Tpetra::Map<int> > row_map = 
	Tpetra::createUniformContigMap<int,int>( global_num_rows, comm );

    // Build the PDF in parallel.
    double pdf_val = 1.0;
    Teuchos::RCP<Tpetra::Vector<double,int> > pdf = 
	Tpetra::createVector<double,int>( row_map );
    pdf->putScalar( pdf_val );

    // Stratify sample the PDF.
    int num_histories = global_num_rows;
    Teuchos::ArrayRCP<int> histories_per_local_state =
	SamplingTools::stratifySampleGlobalPDF( num_histories, pdf );

    // Check the sampling.
    TEST_ASSERT( histories_per_local_state.size() == local_num_rows );

    Teuchos::ArrayRCP<int>::const_iterator history_it;
    for ( history_it = histories_per_local_state.begin();
	  history_it != histories_per_local_state.end();
	  ++history_it )
    {
	TEST_ASSERT( *history_it == 1 );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( SamplingTools, nonuniform_stratify_sample_test )
{
    using namespace Chimera;

    // Setup parallel distribution.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    // Setup global PDF distribution.
    int local_num_rows = 10;
    int global_num_rows = local_num_rows*comm_size;
    Teuchos::RCP<const Tpetra::Map<int> > row_map = 
	Tpetra::createUniformContigMap<int,int>( global_num_rows, comm );

    // Build the PDF in parallel.
    double pdf_val = comm_rank + 1.0;
    Teuchos::RCP<Tpetra::Vector<double,int> > pdf = 
	Tpetra::createVector<double,int>( row_map );
    pdf->putScalar( pdf_val );

    // Stratify sample the PDF.
    int num_histories = 0;
    for ( int i = 0; i < comm_size; ++i )
    {
	num_histories += local_num_rows*(i+1);
    }

    Teuchos::ArrayRCP<int> histories_per_local_state =
	SamplingTools::stratifySampleGlobalPDF( num_histories, pdf );

    // Check the sampling.
    TEST_ASSERT( histories_per_local_state.size() == local_num_rows );

    Teuchos::ArrayRCP<int>::const_iterator history_it;
    for ( history_it = histories_per_local_state.begin();
	  history_it != histories_per_local_state.end();
	  ++history_it )
    {
	TEST_ASSERT( *history_it == pdf_val );
    }
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( SamplingTools, local_pdf_sample_test )
{
    using namespace Chimera;

    // Setup a random number generator.
    typedef boost::mt19937 RNG;
    Teuchos::RCP<RNG> rng = RNGTraits<RNG>::create();

    // Setup local PDF.
    int pdf_size = 10;
    Teuchos::Array<double> pdf_values( pdf_size, 0.0 );
    pdf_values[ pdf_size-1 ] = 1.0;
    Teuchos::Array<int> pdf_indices( pdf_size );
    for ( int i = 0; i < pdf_size; ++i )
    {
	pdf_indices[i] = i+1;
    }

    // Sample the local PDF.
    TEST_ASSERT( pdf_size == SamplingTools::sampleLocalDiscretePDF( 
		     Teuchos::as<Teuchos::ArrayView<const double> >(pdf_values()), 
		     Teuchos::as<Teuchos::ArrayView<const int> >(pdf_indices()), rng ) );

    // Reset the PDF to a different state and sample again.
    std::fill( pdf_values.begin(), pdf_values.end(), 0.0 );
    pdf_values[ pdf_size-3 ] = 1.0;
    TEST_ASSERT( pdf_size-2 == SamplingTools::sampleLocalDiscretePDF( 
		     Teuchos::as<Teuchos::ArrayView<const double> >(pdf_values()), 
		     Teuchos::as<Teuchos::ArrayView<const int> >(pdf_indices()), rng ) );
}

//---------------------------------------------------------------------------//
// end tstBoostRNG.cpp
//---------------------------------------------------------------------------//

