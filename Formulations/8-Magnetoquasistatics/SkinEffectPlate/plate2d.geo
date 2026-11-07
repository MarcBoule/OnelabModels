Include "plate_common.pro";
SetFactory("OpenCASCADE");

Rectangle(1) = {-w/2,-hb/2,0, w,hb};
Rectangle(2) = {-w/2,-h/2,0, w,h};
BooleanFragments{Surface{2}; Delete;}{Surface{1}; Delete;}

Physical Surface("Vol plate", 100) = {2};
Physical Surface("Vol vacuum", 101) = {3,4};
Physical Curve("Sur top bot", 120) = {6, 4};

Mesh.ElementOrder = 1;

MeshSize {:} = 3*meshSize;
MeshSize { PointsOf{ Surface{2}; } } = meshSize;

Mesh.MeshSizeExtendFromBoundary = 1;
Mesh.MeshSizeFromPoints = 1;
Mesh.MeshSizeFromCurvature = 0;

Color Black{ Surface{:};}
