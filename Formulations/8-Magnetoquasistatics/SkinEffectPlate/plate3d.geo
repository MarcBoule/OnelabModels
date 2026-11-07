Include "plate_common.pro";
SetFactory("OpenCASCADE");

Box(1) = {-w/2,-hb/2,-d/2, w,hb,d};
Box(2) = {-w/2,-h/2,-d/2, w,h,d};
BooleanFragments{Volume{2}; Delete;}{Volume{1}; Delete;}

Physical Volume("Vol plate", 100) = {2};
Physical Volume("Vol vacuum", 101) = {3,4};
Physical Surface("Sur vacuum topbot", 120) = {16,20};
Physical Surface("Sur vacuum frontback", 121) = {14,15,19,21};
Physical Surface("Sur leftright", 122) = {7,8,13,17,18,22};
Physical Surface("Sur plate back", 123) = {11};
Physical Surface("Sur plate front", 124) = {12};

Mesh.ElementOrder = 1;

MeshSize {:} = 3*meshSize;
MeshSize { PointsOf{ Volume{2}; } } = meshSize;

Mesh.MeshSizeExtendFromBoundary = 1;
Mesh.MeshSizeFromPoints = 1;
Mesh.MeshSizeFromCurvature = 0;

Color Black{ Surface{:};}
