Include "wire_common.pro";
SetFactory("OpenCASCADE");

Disk(1) = {0,0,0, r};
Disk(2) = {0,0,0, rb};
BooleanFragments{Surface{2}; Delete;}{Surface{1}; Delete;}

Physical Surface("Vol wire", 100) = {1};
Physical Surface("Vol vacuum", 101) = {2};
Physical Curve("Sur ext", 120) = {1};

Mesh.HighOrderOptimize = 2; Mesh.ElementOrder = order;

MeshSize {:} = 4*meshSize;
MeshSize { PointsOf{ Surface{1}; } } = meshSize;

Mesh.MeshSizeExtendFromBoundary = 1;
Mesh.MeshSizeFromPoints = 1;
Mesh.MeshSizeFromCurvature = 0;

Color Black{ Surface{:};}
