//---------------------------------------------------------------------------//
/*!
 * \file VTKWriter.cpp
 * \author Stuart R. Slattery
 * \brief VTK file writer definition for mesh and fields.
 */
//---------------------------------------------------------------------------//

#include <string>
#include <vector>

#include "VTKWriter.hpp"

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
VTKWriter::VTKWriter( const RCP_Comm &comm, 
		      const RCP_Partitioner &partitioner,
		      const RCP_ParameterList &plist )
    : d_comm( comm )
{
    if ( d_comm->getRank() == 0 )
    {
	// Open the file to write.
	d_vtk_file.open( plist->get<std::string>( "OUTPUT_FILENAME" ).c_str() );

	// File header.
	d_vtk_file << "# vtk DataFile Version 2.0\n";
	d_vtk_file << "Chimera 0.1\n";
	d_vtk_file << "ASCII\n";
	d_vtk_file << "DATASET RECTILINEAR_GRID\n";

	// Set the mesh.
	std::vector<double>::const_iterator vec_it;
	int nx = partitioner->getGlobalEdges().first.size();
	int ny = partitioner->getGlobalEdges().second.size();
	int nz = 1;
	d_vtk_file << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";

	d_vtk_file << "X_COORDINATES " << nx << " double\n";
	for ( vec_it = partitioner->getGlobalEdges().first.begin();
	      vec_it != partitioner->getGlobalEdges().first.end();
	      ++vec_it )
	{
	    d_vtk_file << *vec_it << " ";
	}
	d_vtk_file << "\n";

	d_vtk_file << "Y_COORDINATES " << ny << " double\n";
	for ( vec_it = partitioner->getGlobalEdges().second.begin();
	      vec_it != partitioner->getGlobalEdges().second.end();
	      ++vec_it )
	{
	    d_vtk_file << *vec_it << " ";
	}
	d_vtk_file << "\n";

	d_vtk_file << "Z_COORDINATES " << nz << " double\n";
	d_vtk_file << "0\n";
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
VTKWriter::~VTKWriter()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Add a field to the file.
 */
void VTKWriter::addField( const Tpetra::Vector<double>& field )
{

}

//---------------------------------------------------------------------------//
/*!
 * \brief Write the file.
 */
void VTKWriter::write()
{
    if ( d_comm->getRank() == 0 )
    {
	d_vtk_file.close();
    }
}

} // end namespace Chimera

//---------------------------------------------------------------------------//
// end VTKWriter.cpp
//---------------------------------------------------------------------------//

