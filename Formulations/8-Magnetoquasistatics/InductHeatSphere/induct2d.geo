Include "induct_common.pro";
SetFactory("OpenCASCADE");

Point(1) = {0,0,0};
Point(2) = {rs,0,0};
Point(3) = {rt-rc,0,0};
Point(4) = {rt+rc,0,0};
Point(5) = {rb,0,0};
Point(6) = {re,0,0};
Point(7) = {0,rs,0};
Point(8) = {0,rb,0};
Point(9) = {0,re,0};
Point(10) = {rt,0,0}; // center of torus

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Circle(4) = {3,10,4};
Line(5) = {4,5};
Circle(6) = {5,1,8};
Line(7) = {5,6};
Circle(8) = {6,1,9};
Line(9) = {1,7};
Circle(10) = {2,1,7};
Line(11) = {7,8};
Line(12) = {8,9};    

Curve Loop(1) = {1,10,-9}; Plane Surface(2) = {1};// sphere
Curve Loop(3) = {3,-4}; Plane Surface(4) = {3};// coil
Curve Loop(5) = {2,4,5,6,-11,-10}; Plane Surface(6) = {5};// vac int
Curve Loop(7) = {7,8,-12,-6}; Plane Surface(8) = {7};// vac ext

Physical Surface("Vol sphere",  100) = {2};
Physical Surface("Vol coil",    101) = {4};
Physical Surface("Vol vac int", 102) = {6};
Physical Surface("Vol vac ext", 103) = {8};

Physical Curve("Sur left", 121) = {9,11,12}; 
Physical Curve("Sur ext",   122) = {8};

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
Mesh.MeshSizeFromCurvature = 0;

Color Black{ Surface{:};}
