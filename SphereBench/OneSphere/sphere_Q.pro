// How to run: see main.pro

// Conducting (V) full sphere
// This is imposed Q on conductor only, since imposing uniform Q on insulator can be done by imposing an equivalent sigma or rho.
// Also has alternate solution to sphere_V using global quantities when "ImposeV = 1"


ImposeV = 0; // 0 = impose Qf, 1 = impose potential

Qf = sigma_f * 4*Pi*rs^2 / coef;


Group {
	VolExSphere = #{VolAll,-VolSphere};
	If (bound == BOUND_ABC)
		SurDiriV -= #{SurExt};
	EndIf
	DomHGradV = #{VolExSphere,SurSphere,SurExt};
}


Function {
	If (bound == BOUND_IABC) 
		Call EpsIabcDiriScal;
	EndIf
	eps[All] = eps0; // All that are unassigned

	// Exact results (for post analysis):
	If (ImposeV == 0)
		Eex[VolVacInt] = Qf*coef*r[]/(4*Pi*eps0*nr[]^3);
	Else
		Eex[VolVacInt] = V*rs*r[]/nr[]^3;
	EndIf
	We[] = 2*Pi * rs * V^2 * eps0; 
}


Constraint {
	{ Name CstV; Case { // Neumann: any SurExt, SurXZ or SurYZ not in SurDiriV
		{ Region SurDiriV; Value 0; }
	}}
	{ Name CstQs; Case { 
		{ Region SurSphere; Value Qf; }
	}}
	{ Name CstVs; Case { 
		{ Region SurSphere; Value V; }
	}}
}


FunctionSpace {
	{ Name HgradV; Type Form0; // electric scalar potential
		BasisFunction {
			{ Name sn; NameOfCoef vn; Support DomHGradV; 
				Function BF_Node; Entity NodesOf[All, Not SurSphere]; }
			{ Name sf; NameOfCoef vf; Support DomHGradV;
				Function BF_GroupOfNodes; Entity GroupsOfNodesOf[ SurSphere ]; }
		}
		GlobalQuantity {
			{ Name GlobalVs; Type AliasOf       ; NameOfCoef vf; }
			{ Name GlobalQs; Type AssociatedWith; NameOfCoef vf; }
		}
		Constraint {
			{ NameOfCoef vn; EntityType NodesOf; NameOfConstraint CstV; }
			If (ImposeV == 1)
			{ NameOfCoef GlobalVs; EntityType GroupsOfNodesOf; NameOfConstraint CstVs; }
			Else
			{ NameOfCoef GlobalQs; EntityType GroupsOfNodesOf; NameOfConstraint CstQs; }
			EndIf
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
			GlobalTerm { [ -Dof{Q}, {U} ]; In SurSphere; }

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
			// { Name E; Value {Local {
				// [ -{d v} ]; In VolVacInt; Jacobian J1; }}
			// }
			// { Name Eex; Value {Local {
				// [ Eex[] ]; In VolVacInt; Jacobian J1; }}
			// }
			{ Name L2error; Value {Integral {Type Global; 
				[ coef* SquNorm[Eex[]-(-{d v})] ]; // square root in PostOperation
				Integration I2; Jacobian J1; In VolVacInt;}}
			}
			{ Name E2; Value {Integral {Type Global; // E^2 exact integral
				[ coef* SquNorm[Eex[]] ];
				Integration I2; Jacobian J1; In VolVacInt;}}
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
			Print[{prob, quarters, axis, bound, ImposeV}, Format "Prob=%g, Quarters=%g, Axis=%g, Bound=%g, ImposeV=%g:", File > "output.txt"]; 

			Print[ L2error, OnGlobal, StoreInVariable $L2error ];
			Print[ E2, OnGlobal, StoreInVariable $E2 ];
			Print[ {Sqrt[$L2error/$E2]}, Format 
			" RelL2e = %.8g [1]", File > "output.txt" ];

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
			// Print[ E, OnElementsOf VolVacInt, File "sphere_E.pos" ];
			// Print[ Eex, OnElementsOf VolVacInt, File "sphere_Eex.pos" ];
		}
	}
}
