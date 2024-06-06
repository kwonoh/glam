cmake_minimum_required(VERSION 3.5)
project(glam)

set(CMAKE_VERBOSE_MAKEFILE on)

set(CMAKE_CXX_STANDARD 14)

find_package(TBB REQUIRED)
if(TBB_FOUND)
    include_directories(${TBB_INCLUDE_DIRS})
    link_directories(${TBB_LIBRARY_DIRS})
    add_definitions(-DKW_USE_TBB)
else()
    message(STATUS "Unable to find TBB")
endif()

find_package(CGAL COMPONENTS Core)
if (CGAL_FOUND)
    include(${CGAL_USE_FILE})
    include(CGAL_CreateSingleSourceCGALProgram)
    add_definitions(-DKW_USE_CGAL)
else()
    message(STATUS "Shape based metric requires CGAL")
endif()

add_definitions(-DBOOST_LOG_DYN_LINK)
find_package(Boost 1.62.0 COMPONENTS
    graph
    program_options
    REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "Unable to find boost")
endif()

find_package(OpenCL REQUIRED)
if (OpenCL_FOUND)
    include_directories(${OpenCL_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "Unable to find OpenCL")
endif()

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(glam src/main.cpp)
target_link_libraries(glam ${Boost_LIBRARIES} ${OpenCL_LIBRARIES} TBB::tbb)

