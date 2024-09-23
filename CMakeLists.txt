# Minimum CMake version required
cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0167 NEW)

project(MCF)

# Set C++ standard (C++17)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Suppress all CMake developer warnings
set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS 1)

# Set the directory of cmake macro file 
list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(FETCHCONTENT_BASE_DIR "${CMAKE_BINARY_DIR}/_deps")

add_custom_target(clean_cache
    COMMAND ${CMAKE_COMMAND} -E rm -rf ${FETCHCONTENT_BASE_DIR}
    COMMENT "Cleaning FetchContent cache"
)

# Include Google Test setup
# include(gtest)

# Include the libigl setup
include(libigl)

# Enable libigl modules
igl_include(stb)
igl_include(glfw)
igl_include(imgui)
igl_include(opengl)

# Find Eigen3
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS filesystem)

# Set Qhull root directory
set(QHULL_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/external/qhull-install")

# Set include and library directories for Qhull manually using QHULL_ROOT
include_directories(
    ${QHULL_ROOT}/include
    ${EIGEN3_INCLUDE_DIR} 
    ${Boost_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/include)
link_directories(${QHULL_ROOT}/lib)
link_directories(${Boost_LIBRARY_DIRS})

# Output current source directory for debugging
message(CMAKE_CURRENT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR})

# Set build type to Debug by default if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

###########
## Build ##
###########

# Collect source files from remeshing and skeletonization
file(GLOB_RECURSE REMESHING_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/remeshing/*.cpp)
file(GLOB_RECURSE SKELETONIZATION_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/skeletonization/*.cpp)

# Combine all sources into a single list
set(SRCFILES    ${REMESHING_SRC} 
                ${SKELETONIZATION_SRC} 
                ${CMAKE_CURRENT_SOURCE_DIR}/src/utils.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/visualization.cpp)
message("Source files found: ${SRCFILES}")

# Create a library from the source files
add_library(MCF_lib ${SRCFILES})

# Link libraries to the shared code library
target_link_libraries(MCF_lib PUBLIC
    igl::glfw
    igl::imgui
    igl::opengl
    igl::stb
    ${QHULL_ROOT}/lib/libqhullcpp.a
    ${QHULL_ROOT}/lib/libqhull_r.dylib  # Ensure this is the reentrant version
    Eigen3::Eigen
)

# Executable animating the algorithm on an example mesh
add_executable(SkeletonExample ${CMAKE_CURRENT_SOURCE_DIR}/src/main_animate.cpp)
target_link_libraries(SkeletonExample PRIVATE MCF_lib)

# Executable for generating a gallery of skeleton and segmentations
add_executable(SkeletonGallery ${CMAKE_CURRENT_SOURCE_DIR}/src/main_gallery.cpp)
target_link_libraries(SkeletonGallery PRIVATE MCF_lib ${Boost_LIBRARIES})

add_executable(PreprocessData ${CMAKE_CURRENT_SOURCE_DIR}/src/preprocess_data.cpp)
target_link_libraries(PreprocessData PRIVATE MCF_lib ${Boost_LIBRARIES})

message("Include directories: ${QHULL_ROOT}/include, ${EIGEN3_INCLUDE_DIR}, ${CMAKE_CURRENT_SOURCE_DIR}/include")

###########
## Testing ##
###########

# Collect test source files
# file(GLOB_RECURSE TEST_SRCFILES ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)

# # Add test executable
# add_executable(${PROJECT_NAME}_tests ${TEST_SRCFILES})

# # Link the test executable with GoogleTest and the shared library
# target_link_libraries(${PROJECT_NAME}_tests PRIVATE gtest gtest_main MCF_lib)

# # Enable CTest framework for test execution
# enable_testing()

# # Add tests to CTest
# add_test(NAME MCFTest COMMAND ${PROJECT_NAME}_tests)
