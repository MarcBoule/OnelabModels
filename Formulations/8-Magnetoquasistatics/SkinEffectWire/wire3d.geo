Include "wire_common.pro";
SetFactory("OpenCASCADE");

Cylinder(1) = {0,0,0, 0,0,d, r};
Cylinder(2) = {0,0,0, 0,0,d, rb};
BooleanFragments{Volume{2}; Delete;}{Volume{1}; Delete;}

Physical Volume("Vol wire", 100) = {1};
Physical Volume("Vol vacuum", 101) = {2};
Physical Surface("Sur ext", 120) = {4};
Physical Surface("Sur wire back", 121) = {3}; 
Physical Surface("Sur wire front", 122) = {2}; //ReverseMesh Surface {2};
Physical Surface("Sur vacuum frontback", 123) = {5,6}; 

Mesh.HighOrderOptimize = 2; Mesh.ElementOrder = order;

Field[1] = Distance;
Field[1].SurfacesList = {1};
Field[1].Sampling = 200;

Field[2] = MathEval;
Field[2].F = Sprintf("((F1/%g)^0.75) * %g",rb-r,rb-r);

Field[3] = Threshold;
Field[3].InField = 2;
Field[3].DistMin = 0;    Field[3].SizeMin = meshSize; 
Field[3].DistMax = rb-r; Field[3].SizeMax = 4*meshSize; 

Background Field = 3;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

Color Black{ Surface{:};}
//Geometry.Surfaces = 1;
//Geometry.SurfaceLabels = 1;
//Geometry.Normals = 50;
