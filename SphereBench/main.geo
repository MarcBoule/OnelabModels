// How to run: see main.pro

Include "main_common.pro";
SetFactory("OpenCASCADE");

azimut = quarters * Pi/2;           // for sphere
opening = quarters > 1 ? 2*Pi : Pi; // for cylinder
realhc = quarters > 2 ? hc : hc/2;  // for cylinder

Sphere(1)   = {xs,ys,zs, rs, -Pi/2, Pi/2, azimut};
Cylinder(2) = {xc,yc+hc/2,zc, 0,-realhc,0, rc, opening}; 
Sphere(3)   = {0,0,0, rb, -Pi/2, Pi/2, azimut};     // bnd int
Sphere(4)   = {0,0,0, 1.1*rb, -Pi/2, Pi/2, azimut}; // ext1
Sphere(5)   = {0,0,0, 1.2*rb, -Pi/2, Pi/2, azimut}; // ext2
Sphere(6)   = {0,0,0, 1.3*rb, -Pi/2, Pi/2, azimut}; // ext3
Sphere(7)   = {0,0,0, re, -Pi/2, Pi/2, azimut};     // ext4
BooleanFragments{ Volume{7}; Delete; }{Volume{1:6}; Delete; }

Physical Volume("Vol sphere", 101) = {1};
Physical Volume("Vol cylinder", 102) = {2};
Physical Surface("Sur sphere", 121) = {1};

If (quarters == 1)

	Physical Volume("Vol vac int", 103) = {7};
	Physical Volume("Vol vac ext1", 104) = {6};
	Physical Volume("Vol vac ext2", 105) = {5};
	Physical Volume("Vol vac ext3", 106) = {4};
	Physical Volume("Vol vac ext4", 107) = {3};

	Physical Surface("Sur XZ", 122) = {2,5,23,11,14,17,20};
	Physical Surface("Sur YZ", 123) = {3,7,8,22,10,13,16,19};
	Physical Surface("Sur vac ext4", 124) = {9};

ElseIf (quarters == 2)

	Physical Volume("Vol vac int", 103) = {7};
	Physical Volume("Vol vac ext1", 104) = {6};
	Physical Volume("Vol vac ext2", 105) = {5};
	Physical Volume("Vol vac ext3", 106) = {4};
	Physical Volume("Vol vac ext4", 107) = {3};

	Physical Surface("Sur XZ", 122) = {2,3,5,6,8,9,11,12,14,15,17,18,21,22};
	Physical Surface("Sur YZ", 123) = {};
	Physical Surface("Sur vac ext4", 124) = {4};


ElseIf (quarters == 4)

	Physical Volume("Vol vac int", 103) = {4};
	Physical Volume("Vol vac ext1", 104) = {5};
	Physical Volume("Vol vac ext2", 105) = {6};
	Physical Volume("Vol vac ext3", 106) = {7};
	Physical Volume("Vol vac ext4", 107) = {3};

	Physical Surface("Sur XZ", 122) = {};
	Physical Surface("Sur YZ", 123) = {};
	Physical Surface("Sur vac ext4", 124) = {5};
 
EndIf

Mesh.HighOrderOptimize = 2; Mesh.ElementOrder = order;  

// MeshSize{:} = 0.2*cm; // For Feynman sim (comment out block below when used)

// Mesh parameters
Field[1] = MathEval;
Field[1].F = Sprintf("sqrt((x-(%g))^2+(y-(%g))^2+(z-(%g))^2)", xs, ys, zs);
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].DistMin = rs; Field[2].SizeMin = 0.085*cm; 
Field[2].DistMax = rb; Field[2].SizeMax = 0.5*cm;
Background Field = 2;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
