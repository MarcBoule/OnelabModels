// How to run: see main.pro

// Geometrical constants (Gmsh and GetDP)

iabc     = 1;    // 1 = IABC, 0 = shell
quarters = 1;    // 1 = quarter-of-domain, 2 = half, 4 = full
order    = 2;    // geometrical element order and basis function interpolation order (1 or 2)

cm  = 1E-2;      // units
rb  = 5*cm;      // radius of interior boundary
re  = (iabc ? 1.4 : 2) * rb; // radius of exterior (infinite) boundary

// sphere:
rs = 1.0*cm;     // radius
xs = 0*cm * (quarters > 1); // x position of center
ys = 0*cm * (quarters > 2); // y position of center
zs = -2*cm;      // z position of center

// cylinder:
rc = 1.0*cm;     // radius
hc = 4*cm;       // height
xc = 0*cm * (quarters > 1); // x position of center
yc = 0*cm * (quarters > 2); // y position of center
zc = 2*cm;       // z position of center
