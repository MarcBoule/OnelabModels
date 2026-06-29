// How to run: 
// * adjust the desired flags and constants here and in main_common.pro
// * open this file (sphere.pro) in Gmsh
// * press the Run button
// * results are appended to output.txt

Include "main_common.pro";


// Constants
eps0 = 8.8541878188E-12; // ref permittivity (F/m)
mu0  = 1.25663706127E-6; // ref permeability (H/m)
coef = 4.0 / quarters;   // for post-processing integrals

// Problem type
PROB_SPHERE2_COND  = 1;
PROB_SPHERE2_INSUL = 2;
PROB_SPHERE2_MAG   = 3;  // requires quarters >= 2
DefineConstant[prob = {PROB_SPHERE2_MAG, Name "Input/0Problem type",
	Choices{PROB_SPHERE2_COND  = "Sphere2 Cond",
			PROB_SPHERE2_INSUL = "Sphere2 Insul",
			PROB_SPHERE2_MAG   = "Sphere2 Mag"}}];

// Boundary type
BOUND_TRUNC = 1;
BOUND_ABC   = 2;
BOUND_IABC  = 3;
BOUND_SHELL = 4;
DefineConstant[bound = {BOUND_SHELL, Name "Input/1Boundary type",
	Choices{BOUND_TRUNC = "Truncation",
			BOUND_ABC   = "ABC 1st order",
			BOUND_IABC  = "IABC 3rd order",
			BOUND_SHELL = "Shell"}}];

// Simulation parameters
Mp      = 3.2E5; // permanent magnetic dipole moment per unit volume (A/m)
rho_f   = 8E-5;  // volume free charge density (C/m^3)
V       = 100;   // electric potential (V)


Group {
	VolSphere1 = #{101}; 
	VolSphere2 = #{102}; 
	VolVacInt  = #{103}; // interior vacuum
	VolExt1    = #{104}; // exterior shell 1 (inner)
	VolExt2    = #{105}; // exterior shell 2
	VolExt3    = #{106}; // exterior shell 3 (outer)
	SurProbe   = #{121};
	SurExt     = #{122};
	SurSphere1 = #{123};
	SurSphere2 = #{124};
	PtRefPot   = #{141}; // interior reference point for scalar potential gauging

	VolSpheres = #{VolSphere1,VolSphere2};
	VolExts    = #{VolExt1,VolExt2,VolExt3};
	VolAll     = #{VolSpheres,VolVacInt,VolExts};
	VolProbeLayer = ElementsOf[ VolAll, OnPositiveSideOf SurProbe ];
}


Function {
	// IABC coefficients
	If (bound == BOUND_IABC) 
		mu[VolExt1]  = mu0 * 1.04128553077778; 
		mu[VolExt2]  = mu0 * 0.566491614524866;
		mu[VolExt3]  = mu0 * 9.53571368003040;
		eps[VolExt1] = eps0 * 1.04128553077778;
		eps[VolExt2] = eps0 * 0.566491614524866;
		eps[VolExt3] = eps0 * 9.53571368003040;
	EndIf
	eps[All]     = eps0; // All others
	mu[All]      = mu0; // All others

	// Vectors
	// u1 and u2 are magnetization directions for sphere1 and sphere2
	//   they should not have a Y component when quarters < 4, and
	//   they should not have an X component when quarters < 2, and
	u1[] = Unit[Vector[0, 0, 1]]; // for Mp1 vector
	u2[] = Unit[Vector[0, 0, 1]]; // for Mp2 vector
	vd[] = Vector[xs,ys,zs]-Vector[xs2,ys2,zs2]; // distance vector 2 to 1
	d[]  = Norm[vd[]];
	ud[] = vd[] / d[]; 

	c1[]  = Vector[xs,  ys,  zs];  // Center of sphere 1
	c2[]  = Vector[xs2, ys2, zs2]; // Center of sphere 2
	r1[]  = XYZ[] - c1[]; // Position from center of sphere 1
	r2[]  = XYZ[] - c2[]; // Position from center of sphere 2
	nr1[] = Norm[r1[]];
	nr2[] = Norm[r2[]];
}


Constraint {
	{ Name CstUn; Case { 
		{ Region SurProbe;  Value 1; }
		{ Region NodesOf[ VolProbeLayer, Not SurProbe ]; Value 0; }
	}}
}


FunctionSpace {
	{ Name H_Un; Type Form0;
		BasisFunction {
			{ Name sn; NameOfCoef un; Support VolProbeLayer; 
				Function BF_Node; Entity NodesOf[All]; }
		}
		Constraint {
			{ NameOfCoef un; EntityType NodesOf; NameOfConstraint CstUn; }
		}
	}
}


Jacobian {
	{ Name J1; Case {
		{ Region SurExt; Jacobian Sur; }
		If (bound == BOUND_SHELL) 
			{ Region VolExts; Jacobian VolSphShell{rb, re}; }
		EndIf
		{ Region All; Jacobian Vol; }// All that are unassigned
	}}      
}


Integration {
	{ Name I1; Case {{ Type Gauss; Case { 
		{ GeoElement Triangle;     NumberOfPoints 1; }
		{ GeoElement Triangle2;    NumberOfPoints 3; }
		{ GeoElement Tetrahedron;  NumberOfPoints 1; }
		{ GeoElement Tetrahedron2; NumberOfPoints 4; }
	}}}}
	{ Name I2; Case {{ Type Gauss; Case { // for L2 error norm 
		{ GeoElement Tetrahedron;  NumberOfPoints 4; }
		{ GeoElement Tetrahedron2; NumberOfPoints 16; }
	}}}}
}


If (prob == PROB_SPHERE2_COND)
	Include "sphere2_Cond.pro";
ElseIf (prob == PROB_SPHERE2_INSUL)
	Include "sphere2_Insul.pro";
Else //(prob == PROB_SPHERE2_MAG)
	Include "sphere2_Mag.pro";
EndIf
