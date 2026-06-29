// How to run: see main.pro

Include "main_common.pro";
SetFactory("OpenCASCADE");

azimut = quarters * Pi/2;

Sphere(1)   = {xs,ys,zs, rs, -Pi/2, Pi/2, azimut};
Sphere(2)   = {0,0,0, rb, -Pi/2, Pi/2, azimut};     // bnd int
Sphere(3)   = {0,0,0, 4*rb/3, -Pi/2, Pi/2, azimut}; // ext1
Sphere(4)   = {0,0,0, 5*rb/3, -Pi/2, Pi/2, azimut}; // ext2
Sphere(5)   = {0,0,0, re, -Pi/2, Pi/2, azimut};     // ext3
BooleanFragments{ Volume{5}; Delete; }{Volume{1:4}; Delete; }

Physical Volume("Vol sphere", 101) = {1};
Physical Surface("Sur sphere", 121) = {1};

If (quarters == 1)

	Physical Volume("Vol vac int", 102) = {5};
	Physical Volume("Vol vac ext1", 103) = {4};
	Physical Volume("Vol vac ext2", 104) = {3};
	Physical Volume("Vol vac ext3", 105) = {2};

	Physical Surface("Sur XZ", 122) = {2,6,9,12,15};
	Physical Surface("Sur YZ", 123) = {3,5,8,11,14};
	Physical Surface("Sur vac ext3", 124) = {4};

ElseIf (quarters == 2)

	Physical Volume("Vol vac int", 102) = {5};
	Physical Volume("Vol vac ext1", 103) = {4};
	Physical Volume("Vol vac ext2", 104) = {3};
	Physical Volume("Vol vac ext3", 105) = {2};

	Physical Surface("Sur XZ", 122) = {2,3,5,6,8,9,11,12,14,15};
	Physical Surface("Sur YZ", 123) = {};
	Physical Surface("Sur vac ext3", 124) = {4};

ElseIf (quarters == 4)

	Physical Volume("Vol vac int", 102) = {3};
	Physical Volume("Vol vac ext1", 103) = {4};
	Physical Volume("Vol vac ext2", 104) = {5};
	Physical Volume("Vol vac ext3", 105) = {2};

	Physical Surface("Sur XZ", 122) = {};
	Physical Surface("Sur YZ", 123) = {};
	Physical Surface("Sur vac ext3", 124) = {2};

EndIf



Mesh.HighOrderOptimize = 2; Mesh.ElementOrder = order;  

// MeshSize{:} = 0.3*cm; // comment out block below when used

// Mesh parameters
Field[1] = MathEval;
Field[1].F = Sprintf("sqrt((x-(%g))^2+(y-(%g))^2+(z-(%g))^2)", xs, ys, zs);
Field[2] = Threshold;
Field[2].InField = 1;
// FINE:
// Field[2].DistMin = rs; Field[2].SizeMin = 0.085*cm; 
// Field[2].DistMax = re; Field[2].SizeMax = 0.9*cm;
// MEDIUM:
Field[2].DistMin = rs; Field[2].SizeMin = 0.12*cm; 
Field[2].DistMax = re; Field[2].SizeMax = 1.0*cm;
// COARSE:
// Field[2].DistMin = rs; Field[2].SizeMin = 0.2*cm; 
// Field[2].DistMax = re; Field[2].SizeMax = 1.3*cm;
Background Field = 2;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
