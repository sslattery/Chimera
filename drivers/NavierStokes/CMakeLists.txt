INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(SOURCES
  convection_cavity.cpp
  )
  
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ConvectionCavity
  SOURCES ${SOURCES}
  COMM serial mpi
  NUM_MPI_PROCS 4
  )
