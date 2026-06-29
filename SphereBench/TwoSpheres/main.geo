// How to run: see main.pro

Include "main_common.pro";
SetFactory("OpenCASCADE");

azimut = quarters * Pi/2;

Sphere(1)   = {xs,ys,zs, rs, -Pi/2, Pi/2, azimut};  // sphere 1 (top)
Sphere(2)   = {xs,ys,zs, rs*2, -Pi/2, Pi/2, azimut}; // sphere 1 probe surface
Sphere(3)   = {xs2,ys2,zs2, rs, -Pi/2, Pi/2, azimut};  // sphere 2 (bot)

Sphere(4)   = {0,0,0, rb, -Pi/2, Pi/2, azimut};     // bnd int
Sphere(5)   = {0,0,0, 4*rb/3, -Pi/2, Pi/2, azimut}; // ext1
Sphere(6)   = {0,0,0, 5*rb/3, -Pi/2, Pi/2, azimut}; // ext2
Sphere(7)   = {0,0,0, re, -Pi/2, Pi/2, azimut};     // ext3
BooleanFragments{ Volume{7}; Delete; }{Volume{1:6}; Delete; }


Physical Volume("Vol sphere 1", 101) = {1};
Physical Volume("Vol sphere 2", 102) = {3};
Physical Surface("Sur sphere 1", 123) = {1};
Physical Point("Pt ref pot", 141) = {2};  // for potential gauging

If (quarters <= 2)

	Physical Volume("Vol vac int", 103) = {7,8};// 7 is largest volume
	Physical Volume("Vol vac ext1", 104) = {6};
	Physical Volume("Vol vac ext2", 105) = {5};
	Physical Volume("Vol vac ext3", 106) = {4};

	ReverseMesh Surface {22}; 
	Physical Surface("Sur force probe", 121) = {22};
	Physical Surface("Sur vac ext3", 122) = {10};
	Physical Surface("Sur sphere 2", 124) = {7};
	
Else

	Physical Volume("Vol vac int", 103) = {5,6};// 6 is largest volume
	Physical Volume("Vol vac ext1", 104) = {7};
	Physical Volume("Vol vac ext2", 105) = {8};
	Physical Volume("Vol vac ext3", 106) = {4};

	Physical Surface("Sur force probe", 121) = {6};
	Physical Surface("Sur vac ext3", 122) = {4};
	Physical Surface("Sur sphere 2", 124) = {3};
	
EndIf


// Mesh parameters

Mesh.HighOrderOptimize = 2; 
Mesh.ElementOrder = order;  

smin = s*0.085*cm;
smax = s*0.9*cm;

Field[1] = MathEval;
Field[1].F = Sprintf("sqrt((x-(%g))^2+(y-(%g))^2+(z-(%g))^2)", xs, ys, zs);
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].DistMin = rs; Field[2].SizeMin = smin; 
Field[2].DistMax = re; Field[2].SizeMax = smax;

Field[3] = MathEval;
Field[3].F = Sprintf("sqrt((x-(%g))^2+(y-(%g))^2+(z-(%g))^2)", xs2, ys2, zs2);
Field[4] = Threshold;
Field[4].InField = 3;
Field[4].DistMin = rs; Field[4].SizeMin = smin; 
Field[4].DistMax = re; Field[4].SizeMax = smax;

Field[5] = Min;
Field[5].FieldsList = {2, 4};
Background Field = 5;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
