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

PROB_SPHERE2_COND  = 1;
PROB_SPHERE2_INSUL = 2;
PROB_SPHERE2_MAG   = 3;  // requires quarters >= 2

DefineConstant[prob = {PROB_SPHERE2_MAG, Name "Input/0Problem type",
	Choices{PROB_SPHERE2_COND  = "Sphere2 Cond",
			PROB_SPHERE2_INSUL = "Sphere2 Insul",
			PROB_SPHERE2_MAG   = "Sphere2 Mag"}}];

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
	VolExt3    = #{106}; // exterior shell 3
	VolExt4    = #{107}; // exterior shell 4 (outer)
	SurProbe   = #{122};
	SurExt4    = #{123};
	PtRefPot   = #{140}; // interior reference point for potential 

	VolSpheres = #{VolSphere1,VolSphere2};
	VolExts    = #{VolExt1,VolExt2,VolExt3,VolExt4};
	VolAll     = #{VolSpheres,VolVacInt,VolExts};
	VolProbeLayer = ElementsOf[ VolAll, OnPositiveSideOf SurProbe ];
}


Function {
	// IABC coefficients
	If (iabc == 1) 
		eps[VolExt1] = eps0 * 0.801713;
		eps[VolExt2] = eps0 * 2.88849;
		eps[VolExt3] = eps0 * 0.163862;
		eps[VolExt4] = eps0 * 37.8756;
		mu[VolExt1]  = mu0 * 0.801713; 
		mu[VolExt2]  = mu0 * 2.88849;
		mu[VolExt3]  = mu0 * 0.163862;
		mu[VolExt4]  = mu0 * 37.8756;
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
		If (iabc == 0) 
			{ Region VolExts; Jacobian VolSphShell{rb, re}; }// Shell transformation when not IABC
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
}


If (prob == PROB_SPHERE2_COND)
	Include "sphere2_Cond.pro";
ElseIf (prob == PROB_SPHERE2_INSUL)
	Include "sphere2_Insul.pro";
Else //(prob == PROB_SPHERE2_MAG)
	Include "sphere2_Mag.pro";
EndIf
