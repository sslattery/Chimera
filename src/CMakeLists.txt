INCLUDE(TribitsLibraryMacros)

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.hpp
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  Chimera_Assertion.hpp
  Chimera_BoostRNG.hpp
  Chimera_History.hpp
  Chimera_History_def.hpp
  Chimera_HistoryBank.hpp
  Chimera_HistoryBuffer.hpp
  Chimera_HistoryBuffer_def.hpp
  Chimera_JacobiSplit.hpp
  Chimera_JacobiSplit_def.hpp
  Chimera_LinearOperatorSplit.hpp
  Chimera_LinearOperatorSplitFactory.hpp
  Chimera_LinearOperatorSplitFactory_def.hpp
  Chimera_LinearSolverFactory.hpp
  Chimera_LinearSolverFactory_def.hpp
  Chimera_LinearProblem.hpp
  Chimera_LinearProblem_def.hpp
  Chimera_LinearSolver.hpp
  Chimera_LinearSolverFactory.hpp
  Chimera_MCSA.hpp
  Chimera_MCSA_def.hpp
  Chimera_NeumannUlamSolver.hpp
  Chimera_NeumannUlamSolverFactory.hpp
  Chimera_NeumannUlamSolverFactory_def.hpp
  Chimera_OperatorTools.hpp
  Chimera_OperatorTools_def.hpp
  Chimera_OverlapManager.hpp
  Chimera_OverlapManager_def.hpp
  Chimera_RNGTraits.hpp
  Chimera_SamplingTools.hpp
  Chimera_SamplingTools_def.hpp
  Chimera_SequentialMC.hpp
  Chimera_SequentialMC_def.hpp
  Chimera_StationaryIteration.hpp
  Chimera_StationaryIteration_def.hpp
  Chimera_StationarySolver.hpp
  Chimera_StationarySolver_def.hpp
  ) 

APPEND_SET(SOURCES
  Chimera_Assertion.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  chimera
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

INSTALL(FILES 
  ${HEADERS} 
  DESTINATION 
  ${CMAKE_INSTALL_PREFIX}/include/
)
