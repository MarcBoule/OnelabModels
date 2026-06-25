// How to run: see main.pro

// Volumetrically uniformly charged full sphere


Group { 
	If (bound == BOUND_ABC)
		SurDiriV -= #{SurExt};
	EndIf
}


Function {
	If (bound == BOUND_IABC) 
		Call EpsIabcDiriScal;
	EndIf
	eps[VolSphere] = eps0 * epsR;
	eps[All] = eps0; // All that are unassigned

	// Exact results (for post analysis):
	We[] = 2*Pi * rs^5 * rho_f^2 / (9*eps0) * (1/(5*epsR)+1);
}


Constraint {
	{ Name CstV; Case { // Neumann: any SurExt, SurXZ or SurYZ not in SurDiriV
		{ Region SurDiriV; Value 0; }
	}}
}


FunctionSpace {
	{ Name HgradV; Type Form0; // electric scalar potential
		BasisFunction {
			{ Name sn; NameOfCoef vn; Support #{VolAll,SurExt}; 
				Function BF_Node; Entity NodesOf[All]; }
		}
		Constraint {
			{ NameOfCoef vn; EntityType NodesOf; NameOfConstraint CstV; }
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
			
			If (bound == BOUND_ABC)
				// 1st order ABC (n = 1 since monopole is leading harmonic)
				Integral{ [ eps0 * 1 / re * Dof{v}, {v}]; 
				Integration I1; Jacobian J1; In SurExt; }		
			EndIf
		}
	} 
}


Resolution {
	{ Name ResMain;
		System {
			{ Name SV; NameOfFormulation FrmV; }
		}
		Operation {
			Generate[SV]; Solve[SV]; SaveSolution[SV];
		}
	}
}


PostProcessing {
	{ Name PostMain; NameOfFormulation FrmV;
		Quantity {
			{ Name E; Value {Local {
				[ -{d v} ]; In VolAll; Jacobian J1; }}
			}
			{ Name We; Value {Integral {Type Global; 
				[ coef*(eps[]/2) * SquNorm[-{d v}] ]; 
				Integration I1; Jacobian J1; In VolAll;}}
			}
			{ Name We2; Value {Integral {Type Global; 
				[ coef/2 * rho_f * {v} ]; // rho_f for linear dielectrics (Griffiths 5ed p.198)
				Integration I1; Jacobian J1; In VolSphere;}}
			}
		}
	}
}


PostOperation {
	{ Name PostMain; NameOfPostProcessing PostMain; 
		Format Table;
		Operation {
			Print[{prob, quarters, axis, bound, epsR}, Format "Prob=%g, Quarters=%g, Axis=%g, Bound=%g, epsR=%g:", File > "output.txt"]; 

			If (bound != BOUND_ABC)
				Print[ We, OnGlobal, StoreInVariable $We ];
				Print[ {$We, We[], ($We-We[])/We[]*10^6}, Format 
				" We  = %.8g [J] (analyt %.8g, %.3g ppm)", File > "output.txt" ];
			EndIf

			Print[ We2, OnGlobal, StoreInVariable $We2 ];
			Print[ {$We2, We[], ($We2-We[])/We[]*10^6}, Format 
			" We2 = %.8g [J] (analyt %.8g, %.3g ppm)", File > "output.txt" ];
		}
	} 
	{ Name PostFields; NameOfPostProcessing PostMain; 
		Operation {
			//Print[ E, OnElementsOf VolAll, File "sphere_E.pos" ];
		}
	}
}
