// How to run: see main.pro

// Conducting (V) full sphere
// Note: this is imposed Q on conductor face only, since imposing surface Q or volume Q on insulator is same as (and is done by) imposing equivalent uniform sigma or rho.


Qf = sigma_f * 4*Pi*rs^2 / coef;


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
	{ Name CstQs; Case { 
		{ Region SurSphere; Value Qf; }
	}}
}


FunctionSpace {
	{ Name HgradV; Type Form0; // electric scalar potential
		BasisFunction {
			{ Name sn; NameOfCoef vn; Support #{VolExSphere/*,SurSphere*/,SurExt}; 
				Function BF_Node; Entity NodesOf[All, Not SurSphere]; }
			{ Name sf; NameOfCoef vf; Support #{VolExSphere,SurSphere};
				Function BF_GroupOfNodes; Entity GroupsOfNodesOf[ SurSphere ]; }
		}
		GlobalQuantity {
			{ Name GlobalVs; Type AliasOf       ; NameOfCoef vf; }
			{ Name GlobalQs; Type AssociatedWith; NameOfCoef vf; }
		}
		Constraint {
			{ NameOfCoef vn; EntityType NodesOf; NameOfConstraint CstV; }
			{ NameOfCoef GlobalQs; EntityType GroupsOfNodesOf; NameOfConstraint CstQs; }
		}
	}
}


Formulation {
	{ Name FrmV; Type FemEquation;
		Quantity {
			{ Name v; Type Local; NameOfSpace HgradV; }
			{ Name U; Type Global; NameOfSpace HgradV[GlobalVs]; }
			{ Name Q; Type Global; NameOfSpace HgradV[GlobalQs]; }
		}
		Equation {
			Integral { [ eps[] * Dof{d v}, {d v} ]; 
				Integration I1; Jacobian J1; In VolExSphere; }
			GlobalTerm { [ -Dof{Q} , {U} ]; In SurSphere; }

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
			{ Name We2; Value { Term {
				[ coef/2 * {Q} * {U} ]; In SurSphere;}}
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

			Print[ We2, OnRegion SurSphere, StoreInVariable $We2 ];
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
