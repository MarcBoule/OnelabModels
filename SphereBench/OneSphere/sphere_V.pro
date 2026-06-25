// How to run: see main.pro

// Conducting (V) full sphere


Group {
	VolExSphere = #{VolAll,-VolSphere};
	If (bound == BOUND_ABC)
		SurDiriV -= #{SurExt};
	EndIf
}


Function {
	If (bound == BOUND_IABC) 
		Call EpsIabcDiriScal;
	EndIf
	eps[All] = eps0; // All that are unassigned

	// Exact results (for post analysis):
	We[] = 2*Pi * rs * V^2 * eps0; 
}


Constraint {
	{ Name CstV; Case { // Neumann: any SurExt, SurXZ or SurYZ not in SurDiriV
		{ Region SurDiriV; Value 0; }
	}}
}


FunctionSpace {
	{ Name HgradV; Type Form0; // electric scalar potential
		BasisFunction {
			{ Name sn; NameOfCoef vn; Support #{VolExSphere,SurSphere,SurExt}; 
				Function BF_Node; Entity NodesOf[All]; }
		}
		Constraint {
			{ NameOfCoef vn; EntityType NodesOf; NameOfConstraint CstV; }
		}
	}
	{ Name HhalfLam; Type Form0;// Lagrange multiplier
		BasisFunction {
			{ Name sn; NameOfCoef vn; Support SurSphere; 
				Function BF_Node; Entity NodesOf[All];}
		}
	}
}


Formulation {
	{ Name FrmV; Type FemEquation;
		Quantity {
			{ Name v; Type Local; NameOfSpace HgradV; }
			{ Name lam; Type Local; NameOfSpace HhalfLam;}
		}
		Equation {
			Integral { [ eps[] * Dof{d v}, {d v} ]; 
			Integration I1; Jacobian J1; In VolExSphere; }
			Integral{ [ -Dof{lam}, {v}];
			Integration I1; Jacobian J1; In SurSphere; }
			Integral{ [eps0 * Dof{v}, {lam}];
			Integration I1; Jacobian J1; In SurSphere; }
			Integral{ [-eps0 * V, {lam}];
			Integration I1; Jacobian J1; In SurSphere; }

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
				[ -{d v} ]; In VolExSphere; Jacobian J1; }}
			}
			{ Name We; Value {Integral {Type Global; 
				[ coef*(eps[]/2) * SquNorm[-{d v}] ]; 
				Integration I1; Jacobian J1; In VolExSphere;}}
			}
			{ Name We2; Value {Integral {Type Global; 
				[ coef/2 * {lam} * {v} ]; // sigma = lam (Lagrange multiplier)
				Integration I1; Jacobian J1; In SurSphere;}}
			}
		}
	}
}


PostOperation {
	{ Name PostMain; NameOfPostProcessing PostMain; 
		Format Table;
		Operation {
			Print[{prob, quarters, axis, bound}, Format "Prob=%g, Quarters=%g, Axis=%g, Bound=%g:", File > "output.txt"]; 

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
			//Print[ E, OnElementsOf VolExSphere, File "sphere_E.pos" ];
		}
	}
}
