// How to run: 
// * adjust the desired flags and constants here and in main_common.pro
// * open this file (sphere.pro) in Gmsh
// * press the Run button
// * results are appended in output.txt

Include "main_common.pro";

// Constants
eps0 = 8.8541878188E-12; // ref permittivity (F/m)
mu0  = 1.25663706127E-6; // ref permeability (H/m)
coef = 4.0 / quarters;   // for post-processing integrals

PROB_SPHERE_RHO_M   = 1;
PROB_SPHERE_V_M     = 2;
PROB_SPHERE_SIGMA_M = 3;
PROB_SPHERE_P_M     = 4;
PROB_SPHERE_RHO_J   = 5;
PROB_SPHERE_SIGMA_K = 6;
PROB_FEYNMAN        = 7;

DefineConstant[prob = {PROB_SPHERE_RHO_M, Name "Input/0Problem type",
	Choices{PROB_SPHERE_RHO_M   = "Sphere Rho,M",
			PROB_SPHERE_V_M     = "Sphere V,M",
			PROB_SPHERE_SIGMA_M = "Sphere Sigma,M",
			PROB_SPHERE_P_M     = "Sphere P,M",
			PROB_SPHERE_RHO_J   = "Sphere Rho,J",
			PROB_SPHERE_SIGMA_K = "Sphere Sigma,K",
			PROB_FEYNMAN        = "Feynman"}}];

// Simulation parameters
epsR  = 3.1;   // relative permittivity (unitless)
muR   = 1.6;   // relative permeability (unitless)
Pp    = 5E-5;  // electric dipole moment per unit volume (C/m^2)
Mp    = 3.2E5; // magnetic dipole moment per unit volume (A/m)
rho   = 8E-5;  // volume charge density (C/m^3)
V     = 100;   // electric potential (V)
sigma = V * eps0 / rs; // surface charge density (C/m^2)
omega = 1E3;   // angular speed (rad/s)
axis  = 2;     // Mp and omega axis (1=X, 2=Y, 3=Z)
axisP = (axis == 3 ? 1 : axis + 1); // Pp axis (1=X, 2=Y, 3=Z)


Group {
	VolSphere = #{101}; 
	VolCylinder = #{102}; 
	VolVacInt = #{103}; // interior vacuum
	VolExt1   = #{104}; // exterior shell 1 (inner)
	VolExt2   = #{105}; // exterior shell 2
	VolExt3   = #{106}; // exterior shell 3
	VolExt4   = #{107}; // exterior shell 4 (outer)
	SurSphere = #{121}; 
	SurXZ     = #{122}; // cut of XZ plane (Y=0), empty when quarters > 2
	SurYZ     = #{123}; // cut of YZ plane (X=0), empty when quarters > 1
	SurExt4   = #{124}; // curved outer bound (of VolExt4)

	VolExts   = #{VolExt1,VolExt2,VolExt3,VolExt4};
	VolAll    = #{VolSphere,VolCylinder,VolVacInt,VolExts};

	// Electric Dirichlet (can be dipole or monopole),
	//  subproblems can rework this as needed
	SurDiriV   = #{SurExt4}; 

	// Magnetic Dirichlet (dipole only, so do here)
	SurDiriPhi = #{}; // B normal to Dirichlet=0 surfaces
	SurDiriA   = #{SurExt4}; // B parallel to Dirichlet=0 surfaces
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
	mu[VolSphere] = mu0*muR;
	mu[VolCylinder] = mu0;
	mu[VolVacInt] = mu0;
	// Dirichlet-vector and Neumann-scalar:
	mu[VolExt1]  = (iabc ? 0.801713 : 1) * mu0; 
	mu[VolExt2]  = (iabc ? 2.88849 : 1) * mu0;
	mu[VolExt3]  = (iabc ? 0.163862 : 1) * mu0;
	mu[VolExt4]  = (iabc ? 37.8756 : 1) * mu0;

	eps[VolSphere] = eps0*epsR;
	eps[VolCylinder] = eps0;
	eps[VolVacInt] = eps0;
	// VolExt1-VolExt4 done separately in subproblems via calls of:
	Macro EpsNeumannScalar
		eps[VolExt1] = (iabc ? 0.801713 : 1) * eps0;
		eps[VolExt2] = (iabc ? 2.88849 : 1) * eps0;
		eps[VolExt3] = (iabc ? 0.163862 : 1) * eps0;
		eps[VolExt4] = (iabc ? 37.8756 : 1) * eps0;
	Return
	Macro EpsDirichletScalar
		eps[VolExt1] = (iabc ? 1.5363 : 1) * eps0;
		eps[VolExt2] = (iabc ? 0.2472 : 1) * eps0;
		eps[VolExt3] = (iabc ? 14.892 : 1) * eps0;
		eps[VolExt4] = (iabc ? 0.0872062 : 1) * eps0;
	Return

	u[]  = Vector[axis  == 1, axis  == 2, axis  == 3];// for Mp and omega vectors 
	uP[] = Vector[axisP == 1, axisP == 2, axisP == 3];// for Pp vector
	c[]  = Vector[xs, ys, zs]; // Center of sphere
	r[]  = XYZ[] - c[];

	// Position vector for angular momentum calculation (when shell transform)
	//  angular momentum is calculated about the center of the sphere
	nXYZ[] = Norm[XYZ[]];
	tXYZ[] = (nXYZ[] > rb && nXYZ[] < re && iabc == 0) ? 
		XYZ[] * rb*(re-rb) / (nXYZ[]*(re-nXYZ[])) // shell transform
		:
		XYZ[]; // no transform
	tr[] = tXYZ[] - c[];
}


Jacobian {
	{ Name J1; Case { 
		{ Region SurSphere; Jacobian Sur; }
		{ Region #{VolAll,-VolExts}; Jacobian Vol; }
		If (iabc)
			{ Region VolExts; Jacobian Vol; }
		Else
			{ Region VolExts; Jacobian VolSphShell{rb, re}; }
		EndIf
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
ElseIf (prob == PROB_SPHERE_P_M)
	Include "sphere_P_M.pro";
ElseIf (prob == PROB_SPHERE_RHO_J)
	Include "sphere_Rho_J.pro";
ElseIf (prob == PROB_SPHERE_SIGMA_K)
	Include "sphere_Sigma_K.pro";
Else //(prob == PROB_FEYNMAN)
	Include "feynman.pro";
EndIf
