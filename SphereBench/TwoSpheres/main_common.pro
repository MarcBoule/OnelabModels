// How to run: see main.pro

// Geometrical constants (Gmsh and GetDP)

iabc     = 1;    // 1 = IABC, 0 = shell (for truncation, use 0 and comment out shell Jacobian)
quarters = 1;    // 1 = quarter-of-domain, 2 = half, 4 = full
order    = 2;    // geometrical element order and basis function interpolation order (1 or 2)

cm  = 1E-2;      // units
rb  = 5*cm;      // radius of interior boundary
re  = (iabc ? 1.4 : 2) * rb; // radius of exterior (infinite) boundary

// sphere1:
rs = 1.0*cm;     // radius
xs = 0*cm * (quarters > 1); // x position of center
ys = 0*cm * (quarters > 2); // y position of center
zs = 2*cm;       // z position of center

// sphere2 (same radius):
xs2 = 0*cm * (quarters > 1); // x position of center
ys2 = 0*cm * (quarters > 2); // y position of center
zs2 = -2*cm;     // z position of center
