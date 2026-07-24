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
PROB_SPHERE_RHO_M   = 1;
PROB_SPHERE_V_M     = 2;
PROB_SPHERE_SIGMA_M = 3;
PROB_SPHERE_Q_M     = 4;
PROB_SPHERE_P_M     = 5;
PROB_SPHERE_RHO_J   = 6;
PROB_SPHERE_SIGMA_K = 7;
PROB_FEYNMAN        = 8;
DefineConstant[prob = {PROB_SPHERE_RHO_M, Name "Input/0Problem type",
	Choices{PROB_SPHERE_RHO_M   = "Sphere Rho,M",
			PROB_SPHERE_V_M     = "Sphere V,M",
			PROB_SPHERE_SIGMA_M = "Sphere Sigma,M",
			PROB_SPHERE_Q_M     = "Sphere Q,M",
			PROB_SPHERE_P_M     = "Sphere P,M",
			PROB_SPHERE_RHO_J   = "Sphere Rho,J",
			PROB_SPHERE_SIGMA_K = "Sphere Sigma,K",
			PROB_FEYNMAN        = "Feynman"}}];

// Boundary type
BOUND_IABC  = 1;
BOUND_SHELL = 0;
DefineConstant[iabc = {BOUND_IABC, Name "Input/1Boundary type",
	Choices{BOUND_IABC  = "IABC 3rd order",
			BOUND_SHELL = "Shell"}}];
// (for truncation, use 0 and comment out shell Jacobian)


// Simulation parameters
epsR    = 3.1;   // relative permittivity (unitless)
muR     = 1.6;   // relative permeability (unitless)
Pp      = 5E-5;  // permanent electric dipole moment per unit volume (C/m^2)
Mp      = 3.2E5; // permanent magnetic dipole moment per unit volume (A/m)
rho_f   = 8E-5;  // volume free charge density (C/m^3)
V0      = 100;   // electric potential (V)
sigma_f = V0 * eps0 / rs; // surface free charge density (C/m^2)
omega   = 1E3;   // angular speed (rad/s)
axis    = 2;     // Mp and omega axis (1=X, 2=Y, 3=Z)
axisP   = (axis == 3 ? 1 : axis + 1); // Pp axis (1=X, 2=Y, 3=Z)


Group {
	VolSphere = #{101}; 
	VolCylinder = #{102}; 
	VolVacInt = #{103}; // interior vacuum
	VolExt1   = #{104}; // exterior shell 1 (inner)
	VolExt2   = #{105}; // exterior shell 2
	VolExt3   = #{106}; // exterior shell 3 (outer)
	SurSphere = #{121}; 
	SurXZ     = #{122}; // cut of XZ plane (Y=0), empty when quarters > 2
	SurYZ     = #{123}; // cut of YZ plane (X=0), empty when quarters > 1
	SurExt    = #{124}; // curved outer bound (of VolExt3)

	VolExts   = #{VolExt1,VolExt2,VolExt3};
	VolAll    = #{VolSphere,VolCylinder,VolVacInt,VolExts};

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
	// IABC coefficients
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

	// Position vector for angular momentum calculation (when shell transform)
	//  angular momentum is calculated about the center of the sphere
	nXYZ[] = Norm[XYZ[]];
	tXYZ[] = (nXYZ[] > rb && nXYZ[] < re && iabc == 0) ? 
		XYZ[] * rb*(re-rb) / (nXYZ[]*(re-nXYZ[])) // shell transform
		:
		XYZ[]; // no transform
	tr[] = tXYZ[] - c[]; // Transformed position vector from center of sphere
}


Jacobian {
	{ Name J1; Case { 
		{ Region SurSphere; Jacobian Sur; }
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


If (prob == PROB_SPHERE_RHO_M)
	Include "sphere_Rho_M.pro";
ElseIf (prob == PROB_SPHERE_V_M)
	Include "sphere_V_M.pro";
ElseIf (prob == PROB_SPHERE_SIGMA_M)
	Include "sphere_Sigma_M.pro";
ElseIf (prob == PROB_SPHERE_Q_M)
	Include "sphere_Q_M.pro";
ElseIf (prob == PROB_SPHERE_P_M)
	Include "sphere_P_M.pro";
ElseIf (prob == PROB_SPHERE_RHO_J)
	Include "sphere_Rho_J.pro";
ElseIf (prob == PROB_SPHERE_SIGMA_K)
	Include "sphere_Sigma_K.pro";
Else //(prob == PROB_FEYNMAN)
	Include "feynman.pro";
EndIf
