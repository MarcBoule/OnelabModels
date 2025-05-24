// Global constants (gmsh and getdp)
mm  = 1E-3;         // units
ell = 2*mm;         // z length of wirepair
rc  = 0.322*mm;     // radius of conductor
rb  = 1.6*mm;       // radius of int. boundary
re  = 2*rb;         // radius of ext. shell bound.
a   = rc + 0.46*mm; // x position of center of wire
xe  = rb + 0.5*mm;  // x offset for Kelvin region

DefineConstant[
  Flag_3D = {1, Choices{0,1}, Name "Input/1 3D model"}
  Flag_Shell = {0, Choices{0,1}, Name "Input/2 Shell transform"}
  Flag_2ndOrder = {0, Choices{0,1}, Name "Input/3 2nd order elements"}
  Flag_MinPts = {1, Choices{0,1}, Name "Input/4 Min num integration pts"}
];
