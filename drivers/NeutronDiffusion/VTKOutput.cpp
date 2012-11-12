//---------------------------------------------------------------------------//
/*!
 * \file VTKOutput.cpp
 * \author Stuart R. Slattery
 * \brief VTK file output definition for mesh and fields.
 */
//---------------------------------------------------------------------------//

#include <vector>

#include "VTKOutput.hpp"

#include <Exception.hpp>

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_CommHelpers.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
VTKOutput::VTKOutput( const RCP_Comm &comm, 
		      const RCP_Partitioner &partitioner,
		      const RCP_ParameterList &plist )
    : d_comm( comm )
    , d_mesh( partitioner->getMesh() )
{
    // Mesh size.
    d_nx = partitioner->getGlobalEdges().first.size();
    d_ny = partitioner->getGlobalEdges().second.size();

    if ( d_comm->getRank() == 0 )
    {
	// Open the file to write.
	d_vtk_file.open( plist->get<std::string>( "OUTPUT_FILENAME" ).c_str() );

	// File header.
	d_vtk_file << "# vtk DataFile Version 3.0\n";
	d_vtk_file << "Chimera version 0.1\n";
	d_vtk_file << "ASCII\n";
	d_vtk_file << "DATASET RECTILINEAR_GRID\n";
	d_vtk_file << "DIMENSIONS " << d_nx << " " << d_ny << " 1\n";

	// Set the mesh.
	std::vector<double>::const_iterator vec_it;

	d_vtk_file << "X_COORDINATES " << d_nx << " double\n";
	for ( vec_it = partitioner->getGlobalEdges().first.begin();
	      vec_it != partitioner->getGlobalEdges().first.end();
	      ++vec_it )
	{
	    d_vtk_file << *vec_it << " ";
	}
	d_vtk_file << "\n";

	d_vtk_file << "Y_COORDINATES " << d_ny << " double\n";
	for ( vec_it = partitioner->getGlobalEdges().second.begin();
	      vec_it != partitioner->getGlobalEdges().second.end();
	      ++vec_it )
	{
	    d_vtk_file << *vec_it << " ";
	}
	d_vtk_file << "\n";

	d_vtk_file << "Z_COORDINATES 1 double\n";
	d_vtk_file << "0\n";
    }

    d_comm->barrier();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
VTKOutput::~VTKOutput()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Add a field to the file.
 */
void VTKOutput::addField( const int field_type,
			  const RCP_MultiVector &field, 
			  const std::string &name )
{
    // Setup.
    int my_rank = d_comm->getRank();
    int global_num_verts = d_nx*d_ny;
    int global_num_cells = (d_nx-1)*(d_ny-1);
    int global_num_values = 0;
    int num_components = field->getNumVectors();
    int global_length = field->getGlobalLength();
    int local_length = field->getLocalLength();
    int local_x_length = d_mesh->getLocalNumCells().first;
    int local_y_length = d_mesh->getLocalNumCells().second;

    if ( field_type == VERTEX_FIELD )
    {
	global_num_values = global_num_verts;
    }
    else if ( field_type == CELL_FIELD )
    {
	global_num_values = global_num_cells;
    }
    
    // Check the input arguments.
    testPrecondition( field_type == VERTEX_FIELD || 
		      field_type == CELL_FIELD,
		      "Field type not supported." );

    testPrecondition( local_length == local_x_length * local_y_length,
		      "Local field_size != local number of vertices or cells." );

    testPrecondition( local_length == global_length / d_comm->getSize(),
		      "Nonuniform partitioning not supported." );

    testPrecondition( global_length % d_comm->getSize() == 0,
		      "Nonuniform partitioning not supported." );

    testPrecondition( global_length == global_num_values,
		      "Global field size != global number of vertices or cells." );

    // Field header.
    if ( my_rank == 0 )
    {
	if ( field_type == VERTEX_FIELD )
	{
	    d_vtk_file << "POINT_DATA " << global_num_verts << "\n";
	}
	else if ( field_type == CELL_FIELD )
	{
	    d_vtk_file << "CELL_DATA " << global_num_cells << "\n";
	}

	d_vtk_file << "SCALARS " << name << " double " << num_components << "\n";
	d_vtk_file << "LOOKUP_TABLE default\n";
    }

    d_comm->barrier();

    // Create temporary global data vector to copy data into.
    std::vector<double> global_data( global_num_values * num_components );

    // Copy the field into the global data vector. It is assumed the given
    // data is globally unique (does not exist on another process). This will
    // also only work if the multi vector has the same number of local
    // components on each process (checked above). This won't work if the
    // partition sizes are not the same on every process.
    Teuchos::ArrayRCP<const double> local_component_data;
    for ( int n = 0; n < num_components; ++n )
    {
	local_component_data = field->getData( n );

	int x_idx, y_idx, local_idx, total_idx;
	int m = 0;
	for ( int j = 0; j < local_y_length; ++j )
	{
	    for ( int i = 0; i < local_x_length; ++i )
	    {
		x_idx = local_x_length*d_mesh->getBlockID().first + i;
		y_idx = local_y_length*d_mesh->getBlockID().second + j;
		local_idx = x_idx + (d_nx-1)*y_idx;
		total_idx = local_idx*num_components + n;

		global_data[ total_idx ] = local_component_data[m];
		++m;
	    }
	}
    }

    d_comm->barrier();

    // Reduce the global data.
    Teuchos::reduceAll<int,double>( *d_comm,
    				    Teuchos::REDUCE_SUM,
    				    (int) global_data.size(),
    				    &global_data[0],
    				    &global_data[0] );

    // Write to file.
    if ( my_rank == 0 )
    {
	for ( int i = 0; i < global_num_values; ++i )
	{
	    for ( int n = 0; n < num_components; ++n )
	    {
		d_vtk_file << global_data[ num_components*i + n ] << " ";
	    }
	    d_vtk_file << "\n";
	}
    }

    d_comm->barrier();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Write the file.
 */
void VTKOutput::write()
{
    if ( d_comm->getRank() == 0 )
    {
	d_vtk_file.close();
    }
    
    d_comm->barrier();
}

} // end namespace Chimera

//---------------------------------------------------------------------------//
// end VTKOutput.cpp
//---------------------------------------------------------------------------//

