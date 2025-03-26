# Minimum CMake version required
cmake_minimum_required(VERSION 3.10)

# Project name and language
project(CO03plus LANGUAGES CXX CUDA)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add the source files
set(SOURCES
    src/main.cc
    src/DataFile.cc
    src/Track.cc
    src/K0SAnalysis.cc
    src/Vertex.cc
    src/Event.cc
    src/Loop.cc
    src/Helix.cc
    src/CUDATrackKernel.cu
)

# Add the header files
set(HEADERS
    include/DataFile.hh
    include/Event.hh
    include/Track.hh
    include/K0SAnalysis.hh
    include/Vertex.hh
    include/Analysis.hh
    include/Loop.hh
    include/Helix.hh
)

# Include directories
include_directories(include)
find_package(CUDA REQUIRED)
include_directories(${CUDA_INCLUDE_DIRS})
#target_include_directories(CO03plus PRIVATE ${CUDA_INCLUDE_DIRS})

# Lib directories
link_directories(${CUDA_LIBRARY_DIRS})


# Copy the Python script to the build directory
configure_file(${CMAKE_SOURCE_DIR}/scripts/plot_histogram.py ${CMAKE_BINARY_DIR}/plot_histogram.py COPYONLY)

# Create the executable
add_executable(CO03plus ${SOURCES} ${HEADERS})

# Link the CUDA libraries
target_link_libraries(CO03plus cudart)