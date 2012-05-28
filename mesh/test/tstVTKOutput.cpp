//---------------------------------------------------------------------------//
/*!
 * \file tstVTKOutput.cpp
 * \author Stuart R. Slattery
 * \brief VTKOutput unit tests.
 */
//---------------------------------------------------------------------------//

#include <algorithm>

#include "VTKOutput.hpp"
#include "Partitioner.hpp"
#include "Mesh.hpp"

#include <UnitTestHelpers.hpp>

#include <Teuchos_ParameterList.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//

namespace Chimera
{

//---------------------------------------------------------------------------//
// Uniform mesh test.
TEUCHOS_UNIT_TEST( VTKOutput, vtkoutput_test )
{
    // Get the communicator.
    UnitTestHelpers::RCP_Comm comm = UnitTestHelpers::getDefaultComm();
    unsigned int my_rank = comm->getRank();
    unsigned int my_size = comm->getSize();

    // Block setup.
    int i_blocks = 1;
    int j_blocks = 1;
    if ( my_size == 2 )
    {
	j_blocks = 2;
    }
    else if ( my_size == 4 )
    {
	i_blocks = 2;
	j_blocks = 2;
    }
    double x_min = 0.0;
    double y_min = x_min;
    double x_max = 10.0*my_size;
    double y_max = x_max;
    int global_x_num_cells = 6*my_size;
    int global_y_num_cells = 6*my_size;
    int global_num_cells = global_x_num_cells * global_y_num_cells;

    // Setup a parameter list.
    Teuchos::RCP<Teuchos::ParameterList> plist =
	Teuchos::rcp( new Teuchos::ParameterList );
    plist->set( "I_BLOCKS", i_blocks );
    plist->set( "J_BLOCKS", j_blocks );
    plist->set( "GRID_TYPE", "UNIFORM" );
    plist->set( "X_MIN", x_min );
    plist->set( "X_MAX", x_max );
    plist->set( "Y_MIN", y_min );
    plist->set( "Y_MAX", y_max );
    plist->set( "X_NUM_CELLS", global_x_num_cells );
    plist->set( "Y_NUM_CELLS", global_y_num_cells );
    plist->set( "OUTPUT_FILENAME", "test.vtk" );

    // Partition.
    Teuchos::RCP<Partitioner> partitioner = 
	Teuchos::rcp( new Partitioner( comm, plist ) );

    // Open a VTK file.
    VTKOutput vtk_output( comm, partitioner, plist );

    // Add a parallel vector field to the cells.
    Teuchos::RCP<const Tpetra::Map<int> > map = Teuchos::rcp(
	new Tpetra::Map<int>( global_num_cells, 0, comm ) );

    Teuchos::RCP< Tpetra::MultiVector<double> > multi_vector = Teuchos::rcp(
	new Tpetra::MultiVector<double>( map, 2 ) );

    int num_local_cells = multi_vector->getLocalLength();

    for ( int i = 0; i < num_local_cells; ++i )
    {
	multi_vector->replaceLocalValue( i, 0, 1.0*my_rank );
	multi_vector->replaceLocalValue( i, 1, 1.0*my_rank );
    }

    vtk_output.addField( VTKOutput::CELL_FIELD, multi_vector, "TEST_VECTOR" );

    // Write the VTK file.
    vtk_output.write();

    // Read the file back in and check it.
}

} // end namespace Chimera

//---------------------------------------------------------------------------//
// end tstVTKOutput.cpp
//---------------------------------------------------------------------------//

