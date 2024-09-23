# MCFSkeletons
COMP0119 Acquisition and Processing of 3D Geometry

# Project Overview

This project implements the skeleton extraction method described in the following paper:

**Paper Title:** "Mean Curvature Skeletons"  
**Authors:** Andrea Tagliasacchi, Ibraheem Alhashim, Matt Olson, Hao Zhang  
**Journal/Conference:** Computer Graphics Forum  
**Volume:** 31  
**Number:** 5  
**Pages:** 1735--1744  
**Year:** 2012  
**Publisher:** Wiley Online Library

## Features

- **Skeleton Extraction via Mean Curvature Flow:** Implements the algorithm described in the referenced paper for extracting skeletons from 3D meshes.
- **Medial Axis Transform:** Computes the medial axis of a 3D mesh using Voronoi and Delaunay structures.
- **Integration with Qhull:** Utilizes Qhull for geometric computations such as convex hulls and Voronoi diagrams.
- **Automated Dependency Management:** Uses CMake to automatically download and configure necessary libraries, including libigl for geometry processing and Google Test for unit testing.

## Getting Started

### Prerequisites

Before you begin, ensure you have the following tools installed on your system:

- **CMake** (version 3.10 or later)
- **C++ Compiler** that supports C++17
- **Git** (optional, for cloning the repository)

### Included Dependencies

The project includes the following dependencies, which are automatically handled by the build process:

- **libigl:** For geometry processing tasks (included via CMake).
- **Qhull:** For computing Voronoi diagrams and convex hulls (included as an external library in the project).
- **Eigen:** For linear algebra operations (found via CMake).

## Main Executables

The project includes two main executables:

## Main Executables

The project includes two main executables:

### 1. `SkeletonExample`

**Description**:  
`SkeletonExample` is the primary executable that animates the mean curvature skeleton extraction process on a 3D test mesh located in the `data/test_meshes` directory.

**Usage**:  
To run build project
Set qhull external library 'export DYLD_LIBRARY_PATH=/path/to/project/external/qhull-install/lib:$DYLD_LIBRARY_PATH
To run `SkeletonExample`

```bash
./SkeletonExample
