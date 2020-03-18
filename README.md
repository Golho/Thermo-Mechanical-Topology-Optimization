# General
A general purpose topology optimization code for transient heat transfer. The code takes a mesh in the form of a .msh (GMSH file format), version 4.1 is the only currently supported, where boundary conditions are connected with physical groups in the .msh file. This is parser and input to a FEM model. For each optimization problem there is a corresponding class, which takes the FEM model and uses the optimizer to find the optimal solution. 

The code is written with object oriented paradigms in Matlab, relying on a precompiled MEX-file for the MMA-optimizer. 

# Dependencies
Initially developed in Matlab R2019b. To organize code, I'm wrapping everything as a [Matlab Project](https://blogs.mathworks.com/pick/2019/04/19/matlab-projects/). This is a feature new in Matlab R2019a. It is not however crucial, as the alternative is to add all the folder and subfolders to the MATLAB path at startup.

The mesher tool is called [GMSH](http://gmsh.info) and is a lightweight, open-source mesh tool.

All the Matlab dependencies are located in the folder `/dependencies`, except

* `nlopt_optimize` [NLopt Matlab Reference](https://nlopt.readthedocs.io/en/latest/NLopt_Matlab_Reference/) which holds the MMA solver
* NLOPT helper files (such as `NLOPT_LD_MMA`)

# Installation
To install:
1. `git clone` the repository
2. Make sure to compile the NLOPT Matlab interface ([which is created when NLOPT is built](https://nlopt.readthedocs.io/en/latest/NLopt_Installation/))
3. Done!

# Run
When running the code in Matlab, make sure that the dependencies (`/dependencies` + `nlopt_optimize` and NLOPT helper functions) and all the source code (`/src`) plus their sub-folders are added to the Matlab path. If you have the Matlab Projects feature, all but the NLOPT-folder is automatically added to the path at startup.

# How to use
This section describes how you might use the code to find an optimal solution to a certain problem.

*To be continued* 

# Limitations
Currently, the code has some limitations:
* Handles only 2D geometries
* Handles only linear transient heat transfer
* Handles Robin boundary conditions (convection), but only if the boundary condition is constant over time
* Only handles 1 type of elements in the design domain
* Only handles triangular or quadrilateral 2D elements