// How to run: see wirecurr.pro
Include "wirecurr_common.pro";
SetFactory("OpenCASCADE");

order = 2;     // choices: 1 or 2

cz = 0.07*mm;  // half of z-length of cut (now in y)
cd = 0.1*mm;   // depth of cut into conductor
rf = 0.02*mm;  // radius of inner fillet of cut

// Single wire
yf = (cd*cz - rf*Sqrt(cd^2+cz^2))/cz;
xi = cd*rf/Sqrt(cd^2+cz^2);
yi = cz*rf/Sqrt(cd^2+cz^2);

Point(1)  = {rc-yf, 0, 0};
Point(2)  = {rc-yf-yi, -xi, 0};
Point(13) = {rc-yf-yi,  xi, 0};
Point(3)  = {rc, -cz, 0};
Point(12) = {rc,  cz, 0};
Point(6)  = {rc, -ell/2, 0};
Point(9)  = {rc,  ell/2, 0};
Point(7)  = {0, -ell/2, 0};
Point(8)  = {0,  ell/2, 0};

Line(13) = {6,3};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};
Line(14) = {12,9};

// Curve loops and Surfaces
If (hasDefect)
	Line(12) = {12,13};
	Circle(1) = {13,1,2};
	Line(2) = {2,3};
	Curve Loop(1) = {12,1,2,-13,6,7,8,-14};
Else
	Line(15) = {12,3};
	Curve Loop(1) = {-13,6,7,8,-14,15};
EndIf
Plane Surface(1) = {1};

// Physical entities
Physical Curve("Top cond", 120) = {6};
Physical Curve("Bot cond", 121) = {8};
Physical Surface("Wire cond", 100) = {1};

// Mesh parameters
MeshSize{ : } = 0.005*mm;

Mesh.HighOrderOptimize = 2; 
Mesh.ElementOrder = order;  

Color Black{ Surface{:};}
