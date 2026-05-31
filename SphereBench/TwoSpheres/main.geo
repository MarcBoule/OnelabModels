// How to run: see main.pro

Include "main_common.pro";
SetFactory("OpenCASCADE");

azimut = quarters * Pi/2;

Sphere(1)   = {xs,ys,zs, rs, -Pi/2, Pi/2, azimut};  // sphere 1 (top)
Sphere(2)   = {xs,ys,zs, rs*2, -Pi/2, Pi/2, azimut}; // sphere 1 (top) probe surface; one of the points is also a reference point for potential when no potential is imposed
Sphere(3)   = {xs2,ys2,zs2, rs, -Pi/2, Pi/2, azimut};  // sphere 2 (bot)

Sphere(4)   = {0,0,0, rb, -Pi/2, Pi/2, azimut};     // bnd int
Sphere(5)   = {0,0,0, 1.1*rb, -Pi/2, Pi/2, azimut}; // ext1
Sphere(6)   = {0,0,0, 1.2*rb, -Pi/2, Pi/2, azimut}; // ext2
Sphere(7)   = {0,0,0, 1.3*rb, -Pi/2, Pi/2, azimut}; // ext3
Sphere(8)   = {0,0,0, re, -Pi/2, Pi/2, azimut};     // ext4
BooleanFragments{ Volume{8}; Delete; }{Volume{1:7}; Delete; }


Physical Volume("Vol sphere 1", 101) = {1};
Physical Volume("Vol sphere 2", 102) = {3};
//Physical Surface("Sur sphere 1", 121) = {1};

If (quarters <= 2)

	Physical Volume("Vol vac int", 103) = {8,9};// 8 is largest volume
	Physical Volume("Vol vac ext1", 104) = {7};
	Physical Volume("Vol vac ext2", 105) = {6};
	Physical Volume("Vol vac ext3", 106) = {5};
	Physical Volume("Vol vac ext4", 107) = {4};

	ReverseMesh Surface {25}; 
	Physical Surface("Sur force probe", 122) = {25};
	Physical Surface("Sur ext4", 123) = {10};
	Physical Point("Pt ref pot", 140) = {17};
	
Else

	Physical Volume("Vol vac int", 103) = {5,6};// 6 is largest volume
	Physical Volume("Vol vac ext1", 104) = {7};
	Physical Volume("Vol vac ext2", 105) = {8};
	Physical Volume("Vol vac ext3", 106) = {9};
	Physical Volume("Vol vac ext4", 107) = {4};

	Physical Surface("Sur force probe", 122) = {6};
	Physical Surface("Sur ext4", 123) = {4};
	Physical Point("Pt ref pot", 140) = {12};
	
EndIf



Mesh.HighOrderOptimize = 2; Mesh.ElementOrder = order;  


// Mesh parameters

smin = 0.085*cm;
smax = (0.5-0.1*iabc)*cm;

Field[1] = MathEval;
Field[1].F = Sprintf("sqrt((x-(%g))^2+(y-(%g))^2+(z-(%g))^2)", xs, ys, zs);
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].DistMin = rs; Field[2].SizeMin = smin; 
Field[2].DistMax = rb; Field[2].SizeMax = smax;

Field[3] = MathEval;
Field[3].F = Sprintf("sqrt((x-(%g))^2+(y-(%g))^2+(z-(%g))^2)", xs2, ys2, zs2);
Field[4] = Threshold;
Field[4].InField = 3;
Field[4].DistMin = rs; Field[4].SizeMin = smin; 
Field[4].DistMax = rb; Field[4].SizeMax = smax;

Field[7] = Min;
Field[7].FieldsList = {2, 4};
Background Field = 7;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
