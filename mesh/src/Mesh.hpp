//---------------------------------------------------------------------------//
/*!
 * \file Mesh.hpp
 * \author Stuart R. Slattery
 * \brief 2D Cartestian mesh declaration.
 */
//---------------------------------------------------------------------------//

#ifndef CHIMERA_MESH_HPP
#define CHIMERA_MESH_HPP

#include <vector>
#include <utility>

namespace Chimera
{

class Mesh
{
  public:

    //@{
    //! typedefs
    typedef unsigned int                                 size_type;
    typedef unsigned long                                long_type;
    typedef std::vector<double>                          VecDbl;
    typedef std::pair<size_type,size_type>               SizePair;
    typedef std::pair<long_type,long_type>               LongPair;
    typedef std::pair<double,double>                     DblPair;
    typedef std::pair<VecDbl, VecDbl>                    VecPair;
    //@}

    //! Direction enum.
    enum MeshDirection { I = 0, J };

  public:

    // Constructor.
    Mesh( long_type global_Ni, long_type global_Nj,
	  double global_i_min, double global_j_min, 
	  double global_i_max, double global_j_max,
	  size_type I_block, size_type J_block,
	  std::vector<double> i_edges, std::vector<double> j_edges );

    // Destructor.
    ~Mesh();

    //! Get the global number of cells.
    LongPair getGlobalNumCells() const
    { return d_global_N; }

    //! Get the global lower bounds.
    DblPair getGlobalLowerBounds() const
    { return d_lower; }

    //! Get the global upper bounds.
    DblPair getGlobalUpperBounds() const
    { return d_upper; }

    //! Get the block ids.
    SizePair getBlockID() const
    { return d_block; }

    //! Get the local number of cells in the specified direction.
    SizePair getLocalNumCells() const
    { return d_local_N; }

    //! Get the edge vector in the specified direction.
    const VecPair& getLocalEdges() const
    { return d_edges; }

    //! Get the width vector in the specified direction.
    const VecPair& getLocalWidths() const
    { return d_widths; }

  private:

    // Compute cell widths.
    void computeCellWidths();

  private:

    // Global number of cells in I and J directions.
    LongPair d_global_N;

    // Global lower bound.
    DblPair d_lower;

    // Global upper bound.
    DblPair d_upper;
    
    // Block IDs.
    SizePair d_block;

    // Local number of cells in I and J directions.
    SizePair d_local_N;

    // Edges.
    VecPair d_edges;

    // Widths.
    VecPair d_widths;
};

} // end namespace Chimera

#endif // end CHIMERA_MESH_HPP

//---------------------------------------------------------------------------//
// end Mesh.hpp
//---------------------------------------------------------------------------//

