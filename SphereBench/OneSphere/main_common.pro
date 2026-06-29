// How to run: see main.pro

// Geometrical constants (Gmsh and GetDP)

quarters = 1;    // 1 = quarter-of-domain, 2 = half, 4 = full
order    = 2;    // geometrical element order and basis function interpolation order (1 or 2)
s        = 1.5;  // mesh scaling factor: 1.0=fine, 1.5=coarse

cm  = 1E-2;      // units
rb  = 5*cm;      // radius of interior boundary
re  = 2*rb;      // radius of exterior (infinite) boundary

// sphere:
rs = 1.0*cm;     // radius
xs = 0*cm * (quarters > 1); // x position of center (with symmetry guard)
ys = 0*cm * (quarters > 2); // y position of center (with symmetry guard)
zs = -2.0*cm;    // z position of center
