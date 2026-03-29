// How to run: see main.pro

// Spherical charge near a cylindrical magnet (figure 27.6 in  https://www.feynmanlectures.caltech.edu/II_27.html)

// Notes for generating the Feynman image

// parameters: iabc = 0, quarters = 2, order = 1
// in .geo: comment out the scaled mesh and use "MeshSize{:} = 0.2*cm;"
// After main.pro is loaded, select the Feynman problem type
// Press the run button
// In the new view's options, General tab, choose Logarithmic
// In the new view's options, Map tab, press 0 to 9 for different col maps
// Show mesh (Alt+M), then Alt+Shift+V to remove volume mesh (leaves only surface), then Ctrl+Shift+V and select only sphere and cylinder physical groups.
// Tools > Options > Geometry > Visibility: un-check all
// Tools > Options > Mesh > Visibility: check 2D edges and 2D surfaces, 
// Tools > Options > Mesh > Color: then color triangles as needed rgb=0.7,0.7,0.7; Enable lighting; Coloring Mode "By element type"


Group {
	ElemXZ = ElementsOf[ VolVacInt, OnNegativeSideOf SurXZ ]; // For post-processing
}


Function {
	M[VolCylinder] = Mp * u[]; 
	M[All] = Vector[0,0,0]; // All others

	If (iabc == 1) 
		Call EpsIabcDiriScal;
		Call MuIabcDiriVectAndNeumScal;
	EndIf
	
	eps[All] = eps0; // All that are unassigned
	mu[All]  = mu0; // All that are unassigned
}


Constraint {
	{ Name CstV; Case { // Neumann: any SurXZ, SurYZ or SurExt4 not in SurDiriV
		{ Region SurDiriV; Value 0; }
	}}
	{ Name CstPhi; Case { // Neumann: SurExt4 and any SurXZ or SurYZ not in SurDiriPhi
		{ Region SurDiriPhi; Value 0; } 
	}} 
}


FunctionSpace {
	{ Name HgradV; Type Form0; // electric scalar potential
		BasisFunction {
			{ Name sn; NameOfCoef vn; Support VolAll; 
				Function BF_Node; Entity NodesOf[All]; }
		}
		Constraint {
			{ NameOfCoef vn; EntityType NodesOf; NameOfConstraint CstV; }
		}
	}
	{ Name HgradPhi; Type Form0; // magnetic scalar potential
		BasisFunction {
			{ Name sn; NameOfCoef pn;  Support VolAll; 
				Function BF_Node; Entity NodesOf[All]; }
		}
		Constraint {
			{ NameOfCoef pn; EntityType NodesOf; NameOfConstraint CstPhi; }
		}
	}
}


Formulation {
	{ Name FrmV; Type FemEquation;
		Quantity {
			{ Name v; Type Local; NameOfSpace HgradV; }
		}
		Equation {
			Integral { [ eps[] * Dof{d v}, {d v} ]; 
			Integration I1; Jacobian J1; In VolAll; }
			Integral { [ -rho_f, {v} ]; 
			Integration I1; Jacobian J1; In VolSphere; }
		}
	}  
	{ Name FrmPhi; Type FemEquation;
		Quantity {
			{ Name v; Type Local; NameOfSpace HgradV; }
			{ Name p; Type Local; NameOfSpace HgradPhi; }
		}
		Equation {
			Integral { [ mu[] * Dof{d p}, {d p} ];
			Integration I1; Jacobian J1; In VolAll; }
			Integral { [ -M[]*mu0, {d p} ]; // -Br
			Integration I1; Jacobian J1; In VolCylinder; }
		}
	}
}


Resolution {
	{ Name ResMain;
		System {
			{ Name SV; NameOfFormulation FrmV; }
			{ Name SP; NameOfFormulation FrmPhi; }
		}
		Operation {
			Generate[SV]; Solve[SV]; SaveSolution[SV];
			Generate[SP]; Solve[SP]; SaveSolution[SP];
		}
	}
}


PostProcessing {
	{ Name PostMain; NameOfFormulation FrmPhi;
		Quantity {
			{ Name E; Value {Local {
				[ -{d v} ]; In VolVacInt; Jacobian J1; }}
			}
			{ Name B; Value {Local {
				[ mu0*M[]-mu[]*{d p} ]; In VolVacInt; Jacobian J1; }}
			}
			{ Name S; Value {Local {
				[ 1/mu0 * Cross[-{d v}, mu0*M[]-mu[]*{d p}] ]; In ElemXZ; Jacobian J1; }}
			}
			{ Name Lc; Value {Integral {Type Global; // needs iabc=0
				[ coef* eps0* Cross[ tr[], Cross[-{d v}, mu0*M[]-mu[]*{d p}] ] ];
				Integration I1; Jacobian J1; In VolAll;}}
			}
		}
	}
}


PostOperation {
	{ Name PostFields; NameOfPostProcessing PostMain; 
		Operation {
			// Print[ E, OnElementsOf VolVacInt, File "sphere_E.pos" ];    
			// Print[ H, OnElementsOf VolVacInt, File "sphere_H.pos" ];    
			Print[ S, OnElementsOf ElemXZ, File "sphere_S.pos" ];    
		}
	}
	{ Name PostMain; NameOfPostProcessing PostMain; 
		Format Table;
		Operation {
			Print[{prob, quarters, axis, iabc, epsR, muR}, Format "Prob=%g, Quarters=%g, Axis=%g, IABC=%g, epsR=%g, muR=%g:", File > "output.txt"]; 

			If (iabc == 0)
				Print[ Lc, OnGlobal, StoreInVariable $Lc ]; 
				Print[ {CompY[$Lc]}, Format " LcY = %.8g [kg*m^2/s]", File > "output.txt" ];
			Else
				Echo[ " LcY = [needs iabc = 0]", File > "output.txt" ];
			EndIf

			// Echo[ "Lc = ", Format Table, File > "output.txt"] ;
			// Print[ Lc, OnGlobal, Format Table, File > "output.txt" ];
		}
	} 
}
