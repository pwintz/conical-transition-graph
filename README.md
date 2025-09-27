# Conical Transition Graphs
This repository contains an implementation of the Conical Transition Graph for analyzing asymptotic stability in conical hybrid systems. 

## Setup 
This project depends on two other MATLAB packages:

* [`MatGeom`](https://github.com/mattools/matGeom)
* [Analyze N -dimensional Convex Polyhedra version 1.9.0.2 by Matt J](https://www.mathworks.com/matlabcentral/fileexchange/30892-analyze-n-dimensional-convex-polyhedra) - The files from this package have been copied to `+polyhedron/` and modified to fix some issues and generally clean up the code (such usage is permitted by the projects [license](+polyhedron/license.txt)).


## Scripts

The following scripts will run CTG--based analysis for example systems.

* `run_conic_abstraction_with_modes.m` - Used in NA:HS submission.
* `run_conic_abstraction_2D.m`
<!-- * `run_conic_abstraction_3D.m` - not yet fully implemented -->
