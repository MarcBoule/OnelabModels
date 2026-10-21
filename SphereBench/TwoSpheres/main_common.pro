// How to run: see main.pro

// Geometrical constants (Gmsh and GetDP)

cm  = 1E-2;       // units

DefineConstant[
	quarters = 1, // 1 = quarter-of-domain, 2 = half, 4 = full
	order = 2,    // geometrical element order and basis function interpolation order (1 or 2)
	s   = 1.0,    // mesh scaling factor: 1.0=fine, 1.5=coarse

	rb  = 5*cm,   // radius of interior boundary
	re  = 2*rb,   // radius of exterior (infinite) boundary
	rp  = 1.5*cm, // radius of probe surface for sphere 1

	// sphere1:
	rs = 1*cm,    // radius
	xs = 0*cm * (quarters > 1), // x position of center
	ys = 0*cm * (quarters > 2), // y position of center
	zs = 2*cm,    // z position of center

	// sphere2 (same radius):
	xs2 = 0*cm * (quarters > 1), // x position of center
	ys2 = 0*cm * (quarters > 2), // y position of center
	zs2 = -1*cm   // z position of center
];
