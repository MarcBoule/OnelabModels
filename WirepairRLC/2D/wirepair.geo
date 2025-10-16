// How to run: see wirepair.pro
Include "wirepair_common.pro";
SetFactory("OpenCASCADE");

Cylinder(1) = {a,0,0, 0,0,ell, rc,Pi};   // cond.
Cylinder(2) = {a,0,0, 0,0,ell, ri,Pi};   // insul.
Cylinder(3) = {0,0,0, 0,0,ell, rb,Pi/2}; // bnd int
Cylinder(4) = {0,0,0, 0,0,ell, re,Pi/2}; // bnd ext
BooleanFragments{ Volume{4}; Delete; }{ Volume{1:3}; Delete; }

Rectangle (30) = {0,0,0, re, re};
BooleanIntersection(21) = { Volume{1}; Delete;}{ Surface{30}; };
BooleanIntersection(22) = { Volume{2}; Delete;}{ Surface{30}; };
BooleanIntersection(23) = { Volume{3}; Delete;}{ Surface{30}; };
BooleanIntersection(24) = { Volume{4}; Delete;}{ Surface{30}; Delete;};


Physical Surface("Wire cond", 100) = {21};
Physical Surface("Wire insul", 101) = {24};
Physical Surface("Vacuum int", 102) = {23};
Physical Surface("Vacuum ext", 103) = {22};
Physical Curve("Left", 120) = {45, 48};
Physical Curve("Ext", 121) = {43};

// Mesh parameters
MeshSize{:}                      = 0.04*mm;
MeshSize{PointsOf{Curve{43};}}   = 0.08*mm;//ext boundary

Mesh.HighOrderOptimize = 2; 
Mesh.ElementOrder = 2;  

Color Black{ Surface{:};}
