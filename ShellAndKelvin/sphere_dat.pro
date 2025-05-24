// Global constants (gmsh and getdp)
mm  = 1E-3;         // units
rs  = 10*mm;        // radius of conductor
rb  = 20*mm;        // radius of int. boundary
re  = 2*rb;         // radius of ext. shell bound.
xe  = rb + 5*mm;    // x offset for Kelvin region

DefineConstant[
  Flag_3D = {1, Choices{0,1}, Name "Input/1 3D model"}
  Flag_Shell = {0, Choices{0,1}, Name "Input/2 Shell transform"}
  Flag_2ndOrder = {0, Choices{0,1}, Name "Input/3 2nd order elements"}
  Flag_MinPts = {1, Choices{0,1}, Name "Input/4 Min num integration pts"}
];
