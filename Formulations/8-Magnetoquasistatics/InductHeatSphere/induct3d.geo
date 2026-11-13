Include "induct_common.pro";
SetFactory("OpenCASCADE");

azimut = 2*Pi/coef;

Sphere(1) = {0,0,0, rs, 0, Pi/2, azimut};
Torus(2) = {0,0,0, rt, rc, azimut};
Box (3) = {0,0,0, rb, rb, -rb};
BooleanDifference(4) = { Volume{2}; Delete; }{ Volume{3}; Delete; };

Sphere(5) = {0,0,0, rb, 0, Pi/2, azimut};
Sphere(6) = {0,0,0, re, 0, Pi/2, azimut};
BooleanFragments{Volume{6}; Delete;}{Volume{1,4,5}; Delete;}

Physical Volume("Vol sphere",  100) = {1};
Physical Volume("Vol coil",    101) = {4};
Physical Volume("Vol vac int", 102) = {6};
Physical Volume("Vol vac ext", 103) = {5};

Physical Surface("Sur front ex sph", 121) = {7,12,15}; // coil, ext, int
Physical Surface("Sur back ex sph",  122) = {8,10,14}; // coil, ext, int
Physical Surface("Sur ext",   123) = {9};
Physical Surface("Sur front sph", 124) = {20}; 
Physical Surface("Sur back sph",  125) = {21}; 

Mesh.HighOrderOptimize = 2; 
Mesh.ElementOrder = order;  

Field[1] = MathEval;
Field[1].F = Sprintf("sqrt((x-(%g))^2+(y-(%g))^2+(z-(%g))^2)", 0, 0, 0);
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].DistMin = rs; Field[2].SizeMin = meshSize; 
Field[2].DistMax = re; Field[2].SizeMax = meshSize*25;

Field[3] = MathEval;
Field[3].F = Sprintf("sqrt((x-(%g))^2+(y-(%g))^2+(z-(%g))^2)", rt, 0, 0);
Field[4] = Threshold;
Field[4].InField = 3;
Field[4].DistMin = rc; Field[4].SizeMin = meshSize*1.5; 
Field[4].DistMax = re; Field[4].SizeMax = meshSize*25;

Field[5] = Min;
Field[5].FieldsList = {2, 4};
Background Field = 5;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;

Color Black{ Surface{:};}
