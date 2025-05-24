Include "wirepair_dat.pro";
SetFactory("OpenCASCADE");


Cylinder(1) = {a,0,0, 0,0,ell, rc,Pi};  // cond
Cylinder(2) = {0,0,0, 0,0,ell, rb,Pi/2};// bnd int
surMeshDouble = -1;


If(Flag_Shell == 1)


	Cylinder(3) = {0,0,0, 0,0,ell, re,Pi/2};// bnd ext
	BooleanFragments{ Volume{3}; Delete; }{ Volume{2}; Delete; }
	BooleanDifference(4) = { Volume{2}; Delete; }{ Volume{1}; Delete; };

	If (Flag_3D == 1)
	
		Physical Volume("Vol int", 100) = {4};
		Physical Volume("Vol ext", 101) = {3};

		Physical Surface("Cond", 120) = {16};
		Physical Surface("Lefts", 121) = {8,13};
		Physical Surface("Bnd int", 122) = {11};// unused
		surMeshDouble = 6;
		Physical Surface("Bnd ext", 123) = {surMeshDouble};
	
	Else // 2D
	
		Rectangle (20) = {0,0,0, 3*rb, 3*rb};
		BooleanIntersection(21) = { Volume{4}; Delete;}{ Surface{20}; };
		BooleanIntersection(22) = { Volume{3}; Delete;}{ Surface{20}; Delete;};
		
		Physical Surface("Vol int", 100) = {21};
		Physical Surface("Vol ext", 101) = {22};

		Physical Curve("Cond", 120) = {38};
		Physical Curve("Lefts", 121) = {40,43};
		Physical Curve("Bnd int", 122) = {20};// unused
		surMeshDouble = 41;
		Physical Curve("Bnd ext", 123) = {surMeshDouble};
	
	EndIf


Else // Kelvin


	BooleanDifference(3) = { Volume{2}; Delete; }{ Volume{1}; Delete; };
	Cylinder(4) = {xe,0,0, 0,0,ell, rb,Pi/2};// bnd ext

	If (Flag_3D == 1)
	
		Physical Volume("Vol int", 100) = {3};
		Physical Volume("Vol ext", 101) = {4};

		Physical Surface("Cond", 120) = {6};
		Physical Surface("Lefts", 121) = {3,12};
		Physical Surface("Bnd int", 122) = {1};
		Physical Surface("Bnd ext", 123) = {8};

		Periodic Surface {8} = {1} Translate {xe, 0, 0};

	Else // 2D
	
		Rectangle (20) = {0,0,0, 3*rb, 3*rb};
		BooleanIntersection(21) = { Volume{4}; Delete;}{ Surface{20}; };
		BooleanIntersection(22) = { Volume{3}; Delete;}{ Surface{20}; Delete;};
	
		Physical Surface("Vol int", 100) = {22};
		Physical Surface("Vol ext", 101) = {21};

		Physical Curve("Cond", 120) = {34};
		Physical Curve("Lefts", 121) = {30,36};
		Physical Curve("Bnd int", 122) = {32};
		Physical Curve("Bnd ext", 123) = {29};

		Periodic Curve {29} = {32} Translate {xe, 0, 0};
	
	EndIf
	

EndIf


Mesh.ElementOrder = 1 + Flag_2ndOrder; 
mSize = 0.04*mm * (1 + Flag_2ndOrder);
MeshSize{:} = mSize;
If (surMeshDouble != -1) 
	If (Flag_3D)
		MeshSize{ PointsOf{ Surface{surMeshDouble}; } } = 2*mSize;
	Else
		MeshSize{ PointsOf{ Curve{surMeshDouble}; } } = 2*mSize;
	EndIf
EndIf


Color Black{ Surface{:};}
