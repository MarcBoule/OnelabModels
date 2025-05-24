# Shell and Kelvin transformation benchmarking

This directory contains examples for studying the various shell and Kelvin (inversion) transformations in GetDP. It covers 3D-Spherical, 3D-Cylindrical, 2D-Planar and 2D-Axisymmetric cases on Form0 function spaces. 

More specifically, the electric potential is simulated for a pair of perfectly conducting parallel cylindrical conductors (called wirepair) or a single perfectly conducting sphere. In both cases the electric field energy in the entire unbounded space is computed, from which the capacitance can then be established. For the cylindrical case, it is the capacitance between the two wires, and in the spherical case, it is the capacitance with respect to infinity (i.e. an outer shell at infinity).

## General remarks

### Implementation note

Support for the Kelvin transformation was added in the development version of GetDP in early march 2025, and requires a version of GetDP _greater_ than 3.5.0 (or a development version until such a version is released).

### Meshing

The meshing is kept as simple as possible and is generally uniform. In the interior region it is identical in both the Shell and Kelvin cases. In the exterior region, the mesh becomes coarser near the outer boundary for the Shell case, while it is kept uniform in the Kelvin case. This is done in order to keep the number of degrees of freedom generally comparable between both cases.

### Number of integration points

An option is provided to study the effect of the number of integration points used in the Integration object in the GetDP simulation code. 

When the "Min num of integration pts" option is checked, the theoretical minimum number of points is used. For first order simulations, triangles and tetrahedra can be properly integrated with a single point, while for second order simulations, three points are needed for triangles and four points are needed for tetrahedra. The results shown below show, interestingly, that the choice of points can have a non-negligible effect in some cases. 
When the "Min num of integration pts" option is not checked, the number of points correspond to those given in the Lib_Elasticity_u.pro template that is bundled with GetDP.

### Second order elements

When the "2nd order elements" option is checked, both second order geometrical elements are used, as well as second order interpolation basis functions. The only visible effect of this parameter in the models is in the .geo files, since for Form0 spaces the BF_Node_2E basis functions are automatically added in GetDP (and thus are not explicitly visible in the .pro files' BasisFunction declarations).


## 3D Spherical







