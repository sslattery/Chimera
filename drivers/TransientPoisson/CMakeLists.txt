INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TransientPoisson_Aztec
  DEPLIBS chimera panzer
  SOURCES   TransientPoisson_Aztec_main.cpp
  COMM serial mpi
  NUM_MPI_PROCS 4
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TransientPoisson_Chimera
  DEPLIBS chimera panzer
  SOURCES   TransientPoisson_Chimera_main.cpp
  COMM serial mpi
  NUM_MPI_PROCS 4
  )
