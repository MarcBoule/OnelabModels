// How to run: see wirepair.pro
Include "wirepair_common.pro";
SetFactory("OpenCASCADE");

// See end of file for mesh scaling

cz = 0.07*mm; // half of z-length of cut (defect)
cd = 0.1*mm;  // depth of cut into conductor
rf = 0.02*mm; // radius of inner fillet of cut

If (!hasDefect)

	Cylinder(1) = {a,0,0, 0,0,ell, rc,Pi};  // cond.
	Cylinder(2) = {a,0,0, 0,0,ell, ri,Pi};  // insul.
	Cylinder(3) = {0,0,0, 0,0,ell, rb,Pi/2};// bnd int
	Cylinder(4) = {0,0,0, 0,0,ell, re,Pi/2};// bnd ext
	BooleanFragments{ Volume{4}; Delete; }{ Volume{1:3}; Delete; }

	Physical Volume("Wire cond", 100) = {1};
	Physical Volume("Wire insul", 101) = {4};
	Physical Volume("Vacuum int", 102) = {3};
	Physical Volume("Vacuum ext", 103) = {2}; // shell
	Physical Surface("Left", 120) = {3, 11};
	Physical Surface("Front", 121) = {2, 7, 13, 18};
	Physical Surface("Back", 122) = {4, 12, 17, 19};
	Physical Surface("Ext", 123) = {1}; // curved r=re

Else

	// conductor and insulator (extrusion)
	yf = (cd*cz - rf*Sqrt(cd^2+cz^2))/cz;
	xi = cd*rf/Sqrt(cd^2+cz^2);
	yi = cz*rf/Sqrt(cd^2+cz^2);
	zi = cz*(ri-rc)/cd;

	Point(1)  = {a+rc-yf, 0, 0};
	Point(2)  = {a+rc-yf-yi, 0, -xi};
	Point(13) = {a+rc-yf-yi, 0,  xi};
	Point(3)  = {a+rc, 0, -cz};
	Point(12) = {a+rc, 0,  cz};
	Point(4)  = {a+ri, 0, -cz-zi};
	Point(11) = {a+ri, 0,  cz+zi};
	Point(5)  = {a+ri, 0, -ell/2};
	Point(10) = {a+ri, 0,  ell/2};
	Point(6)  = {a+rc, 0, -ell/2};
	Point(9)  = {a+rc, 0,  ell/2};
	Point(7)  = {a, 0, -ell/2};
	Point(8)  = {a, 0,  ell/2};

	Circle(1) = {13,1,2};
	Line(2) = {2,3};
	Line(3) = {3,4};
	Line(4) = {4,5};
	Line(5) = {5,6};
	Line(6) = {6,7};
	Line(7) = {7,8};
	Line(8) = {8,9};
	Line(9) = {9,10};
	Line(10) = {10,11};
	Line(11) = {11,12};
	Line(12) = {12,13};
	Line(13) = {6,3};
	Line(14) = {12,9};

	Curve Loop(1) = {3,4,5,13};// insul back
	Curve Loop(2) = {10,11,14,9};// insul front
	Curve Loop(3) = {1,2,-13,6,7,8,-14,12};// cond

	Plane Surface(1) = {1};
	Plane Surface(2) = {2};
	Plane Surface(3) = {3};

	Extrude { {0,0,1}, {a,0,0}, Pi}{ Surface{1,2,3}; }

	// empty space 
	Cylinder(12) = {0,0,-ell/2, 0,0,ell, rb, Pi/2};
	Cylinder(13) = {0,0,-ell/2, 0,0,ell, re, Pi/2};

	// fragment all relevant volumes
	BooleanFragments{ Volume{13}; Delete; }{ Volume{1,2,3,12}; Delete; }

	Physical Volume("Wire cond", 100) = {3};
	Physical Volume("Wire insul", 101) = {1,2};
	Physical Volume("Vacuum int", 102) = {5};
	Physical Volume("Vacuum ext", 103) = {4};
	Physical Surface("Left", 120) = {22,29};
	Physical Surface("Front", 121) = {21,26,12,17};
	Physical Surface("Back", 122) = {23,30,6,16};
	Physical Surface("Ext", 123) = {20};

EndIf 


// Mesh parameters
s = 1.5;// scaling factor: 1.0=fine, 1.5=coarse

If (!hasDefect)

	MeshSize{:}                    = s*0.08*mm;
	MeshSize{PointsOf{Volume{1};}} = s*0.04*mm;// cond.
	MeshSize{PointsOf{Volume{2};}} = s*0.12*mm;// shell
	MeshSize{PointsOf{Surface{1};}}= s*0.20*mm;// ext boundary

Else

	MeshSize{:}                     = s*0.08*mm;
	MeshSize{PointsOf{Volume{3};}}  = s*0.04*mm; // cond.
	MeshSize{PointsOf{Volume{4};}}  = s*0.12*mm; // shell
	MeshSize{PointsOf{Surface{20};}}= s*0.20*mm; // ext boundary
	MeshSize{PointsOf{Surface{14};}}= 1*0.008*mm;// fillet

EndIf

Mesh.HighOrderOptimize = 2; 
Mesh.ElementOrder = 2;  

Color Black{ Surface{:};}
