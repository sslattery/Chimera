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
    int global_num_cells = 10*my_size;

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
    plist->set( "X_NUM_CELLS", global_num_cells );
    plist->set( "Y_NUM_CELLS", global_num_cells );
    plist->set( "OUTPUT_FILENAME", "test.vtk" );

    // Partition.
    Teuchos::RCP<Partitioner> partitioner = 
	Teuchos::rcp( new Partitioner( comm, plist ) );

    // Write a VTK file.
    VTKOutput vtk_output( comm, partitioner, plist );
    vtk_output.write();
}

} // end namespace Chimera

//---------------------------------------------------------------------------//
// end tstVTKOutput.cpp
//---------------------------------------------------------------------------//

