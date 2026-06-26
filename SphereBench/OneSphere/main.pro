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
PROB_SPHERE_P     = 1;
PROB_SPHERE_RHO   = 2;
PROB_SPHERE_V     = 3;
PROB_SPHERE_Q     = 4;
PROB_SPHERE_SIGMA = 5;
PROB_SPHERE_M     = 6;
PROB_SPHERE_J     = 7;
PROB_SPHERE_K     = 8;
DefineConstant[prob = {PROB_SPHERE_J, Name "Input/0Problem type",
	Choices{PROB_SPHERE_P       = "Sphere P",
			PROB_SPHERE_RHO     = "Sphere Rho",
			PROB_SPHERE_V       = "Sphere V",
			PROB_SPHERE_Q       = "Sphere Q",
			PROB_SPHERE_SIGMA   = "Sphere Sigma",
			PROB_SPHERE_M       = "Sphere M",
			PROB_SPHERE_J       = "Sphere J",
			PROB_SPHERE_K       = "Sphere K"}}];

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
epsR    = 3.1;   // relative permittivity (unitless)
muR     = 1.6;   // relative permeability (unitless)
Pp      = 5E-5;  // permanent electric dipole moment per unit volume (C/m^2)
Mp      = 3.2E5; // permanent magnetic dipole moment per unit volume (A/m)
rho_f   = 8E-5;  // volume free charge density (C/m^3)
V       = 100;   // electric potential (V)
sigma_f = V * eps0 / rs; // surface free charge density (C/m^2)
omega   = 1E3;   // angular speed (rad/s)
axis    = 3;     // Mp and omega axis (1=X, 2=Y, 3=Z)
axisP   = 3;     // Pp axis (1=X, 2=Y, 3=Z)


Group {
	VolSphere = #{101}; 
	VolVacInt = #{102}; // interior vacuum
	VolExt1   = #{103}; // exterior shell 1 (inner)
	VolExt2   = #{104}; // exterior shell 2
	VolExt3   = #{105}; // exterior shell 3 (outer)
	SurSphere = #{121}; 
	SurXZ     = #{122}; // cut of XZ plane (Y=0), empty when quarters > 2
	SurYZ     = #{123}; // cut of YZ plane (X=0), empty when quarters > 1
	SurExt    = #{124}; // curved outer bound (of VolExt3)

	VolExts   = #{VolExt1,VolExt2,VolExt3};
	VolAll    = #{VolSphere,VolVacInt,VolExts};

	// Electric Dirichlet (can be dipole or monopole),
	//  subproblems can rework this as needed
	SurDiriV   = #{SurExt}; 

	// Magnetic Dirichlet (dipole only, so do here)
	SurDiriPhi = #{}; // B normal to Dirichlet=0 surfaces
	SurDiriA   = #{SurExt}; // B parallel to Dirichlet=0 surfaces
	If (axis == 1)
		SurDiriPhi += #{SurYZ};
		SurDiriA   += #{SurXZ};
	ElseIf (axis == 2)
		SurDiriPhi += #{SurXZ};
		SurDiriA   += #{SurYZ};
	Else // (axis == 3)
		SurDiriA   += #{SurXZ,SurYZ};
	EndIf
}


Function {
	// IABC 3rd order coefficients, see IABC3_Axi.sage
	Macro MuIabcDiriVectAndNeumScal
		mu[VolExt1]  = mu0 * 1.04128553077778; 
		mu[VolExt2]  = mu0 * 0.566491614524866;
		mu[VolExt3]  = mu0 * 9.53571368003040;
	Return
	Macro EpsIabcNeumScal
		eps[VolExt1] = eps0 * 1.04128553077778;
		eps[VolExt2] = eps0 * 0.566491614524866;
		eps[VolExt3] = eps0 * 9.53571368003040;
	Return
	Macro EpsIabcDiriScal
		eps[VolExt1] = eps0 * 0.863724546405317;
		eps[VolExt2] = eps0 * 3.13613101072653;
		eps[VolExt3] = eps0 * 0.150891868561234;
	Return


	// Vectors
	u[]  = Vector[axis  == 1, axis  == 2, axis  == 3];// for Mp and omega vectors 
	uP[] = Vector[axisP == 1, axisP == 2, axisP == 3];// for Pp vector
	c[]  = Vector[xs, ys, zs]; // Center of sphere
	r[]  = XYZ[] - c[]; // Position from center of sphere
	nr[] = Norm[r[]];
}


Jacobian {
	{ Name J1; Case { 
		{ Region #{SurSphere, SurExt}; Jacobian Sur; }
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


If (prob == PROB_SPHERE_P)
	Include "sphere_P.pro";
ElseIf (prob == PROB_SPHERE_RHO)
	Include "sphere_Rho.pro";
ElseIf (prob == PROB_SPHERE_V)
	Include "sphere_V.pro";
ElseIf (prob == PROB_SPHERE_Q)
	Include "sphere_Q.pro";
ElseIf (prob == PROB_SPHERE_SIGMA)
	Include "sphere_Sigma.pro";
ElseIf (prob == PROB_SPHERE_M)
	Include "sphere_M.pro";
ElseIf (prob == PROB_SPHERE_J)
	Include "sphere_J.pro";
Else //(prob == PROB_SPHERE_K)
	Include "sphere_K.pro";
EndIf
