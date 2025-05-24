Include "sphere_dat.pro";
SetFactory("OpenCASCADE");


Sphere(1) = {0,0,0, rs, 0, Pi/2, Pi/2};// cond
Sphere(2) = {0,0,0, rb, 0, Pi/2, Pi/2};// bnd int
surMeshDouble = -1;


If(Flag_Shell == 1)

	
	Sphere(3) = {0,0,0, re, 0, Pi/2, Pi/2};// bnd ext
	BooleanFragments{ Volume{3}; Delete; }{ Volume{2}; Delete; }
	BooleanDifference(4) = { Volume{2}; Delete; }{ Volume{1}; Delete; };

	If (Flag_3D == 1)

		Physical Volume("Vol int", 100) = {4};
		Physical Volume("Vol ext", 101) = {3};

		Physical Surface("Cond", 120) = {13};
		Physical Surface("Bnd int", 121) = {9};// unused
		surMeshDouble = 5;
		Physical Surface("Bnd ext", 122) = {surMeshDouble};

		Physical Point("Pnt inf", 140) = {5};// unused
		
	Else // 2D
	
		Rectangle (20) = {0,0,0, 3*rb, 3*rb};
		BooleanIntersection(21) = { Volume{4}; Delete;}{ Surface{20}; };
		BooleanIntersection(22) = { Volume{3}; Delete;}{ Surface{20}; Delete;};

		Physical Surface("Vol int", 100) = {21};
		Physical Surface("Vol ext", 101) = {22};

		Physical Curve("Cond", 120) = {31};
		Physical Curve("Bnd int", 121) = {16};// unused
		surMeshDouble = 33;
		Physical Curve("Bnd ext", 122) = {surMeshDouble};

		Physical Point("Pnt inf", 140) = {21};// unused
		
	EndIf


Else // Kelvin


	BooleanDifference(3) = { Volume{2}; Delete; }{ Volume{1}; Delete; };
	Sphere(4) = {xe,0,0, rb, 0, Pi/2, Pi/2};// bnd ext

	If (Flag_3D == 1)
	
		Physical Volume("Vol int", 100) = {3};
		Physical Volume("Vol ext", 101) = {4};

		Physical Surface("Cond", 120) = {5};
		Physical Surface("Bnd int", 121) = {1};
		Physical Surface("Bnd ext", 122) = {6};

		Periodic Surface {6} = {1} Translate {xe, 0, 0};

		Physical Point("Pnt inf", 140) = {10};

	Else // 2D
	
		Rectangle (20) = {0,0,0, 3*rb, 3*rb};
		BooleanIntersection(21) = { Volume{4}; Delete;}{ Surface{20}; };
		BooleanIntersection(22) = { Volume{3}; Delete;}{ Surface{20}; Delete;};
		
		Physical Surface("Vol int", 100) = {22};
		Physical Surface("Vol ext", 101) = {21};

		Physical Curve("Cond", 120) = {28};
		Physical Curve("Bnd int", 121) = {26};
		Physical Curve("Bnd ext", 122) = {23};

		Periodic Curve {23} = {26} Translate {xe, 0, 0};

		Physical Point("Pnt inf", 140) = {17};
	
	EndIf
	

EndIf


Mesh.ElementOrder = 1 + Flag_2ndOrder; 
mSize = 0.47*mm * (1 + Flag_2ndOrder);
MeshSize{:} = mSize;
If (surMeshDouble != -1) 
	If (Flag_3D)
		MeshSize{ PointsOf{ Surface{surMeshDouble}; } } = 2*mSize;
	Else
		MeshSize{ PointsOf{ Curve{surMeshDouble}; } } = 2*mSize;
	EndIf
EndIf


Color Black{ Surface{:};}
