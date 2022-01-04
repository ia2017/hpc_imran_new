cmake_minimum_required(VERSION 3.10) 


project (runliddrivencavity) 
add_executable(runliddrivencavity main.cpp LidDrivenCavity.cpp PoissonSolver.cpp)


target_include_directories( runliddrivencavity PUBLIC “${PROJECT_BINARY_DIR}”)


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)

find_package(BLAS)
find_package(LAPACK)
find_package( Boost 1.40 COMPONENTS program_options REQUIRED )
find_package(MPI REQUIRED)

if(LAPACK_FOUND AND BLAS_FOUND)
   set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()


if(Boost_FOUND)
include_directories( ${Boost_INCLUDE_DIR} )
endif()

include_directories(${MPI_INCLUDE_PATH})

target_link_libraries( runliddrivencavity ${lapackblas_libraries} ${Boost_LIBRARIES} ${MPI_LIBRARIES})

set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})


