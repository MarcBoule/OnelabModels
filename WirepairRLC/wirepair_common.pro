// How to run: see wirepair.pro
// Geometrical constants (for Gmsh and GetDP)
mm  = 1E-3;         // units
ell = 4*mm;         // z length of segment
rc  = 0.322*mm;     // radius of conductor
ri  = rc + 0.41*mm; // outer radius of insulator
rb  = 2*mm;         // radius of inner boundary
re  = 2*rb;         // radius of external boundary
a   = ri + 0.05*mm; // x offset for center of wire

hasDefect = 0;      // "V"-shaped cut into wire     
