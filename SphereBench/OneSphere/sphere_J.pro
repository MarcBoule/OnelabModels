// How to run: see main.pro

// Spinning uniformly volume-charged nonconducting full sphere
// epsilon_r = mu_r = 1 (or else there will be spinning bound charges and currents which would complicate things)


Group {
	VolExt3Shell = ElementsOf[ VolExt3, OnNegativeSideOf SurExt ];
	VolCoul = #{VolSphere,VolExt3Shell};
	If (bound == BOUND_ABC)
		SurDiriA -= #{SurExt};
	EndIf
}


Function {
	v[] = Cross[ omega*u[], r[] ];
	J[] = rho_f * v[];

	If (bound == BOUND_IABC) 
		Call MuIabcDiriVectAndNeumScal;
	EndIf
	mu[All]  = mu0; // All that are unassigned

	// Exact results (for post analysis):
	Aex[VolSphere] = (5*rs^2 - 3*nr[]^2) * J[] * (mu0/30);
	Aex[VolVacInt] = rs^5 * J[] * (mu0/(15*nr[]^3));
	Wb[] = 8*Pi  * rs^7 * rho_f^2 * omega^2 * (mu0/315);
}


Constraint {
	{ Name CstA; Case { // Neumann: any SurXZ or SurYZ not in SurDiriA
		{ Region SurDiriA; Value 0; }
	}} 
	{ Name CstGaugeA; Case {
		{ Region VolAll; SubRegion SurDiriA; Value 0; }
	}}
	// Boundary condition for the Coulomb gauge multiplier "xi" (only used when
	// "bound == BOUND_ABC" or when L2 norm is wanted):
	{ Name xi_Mag; Case {
		{ Region #{SurDiriA,SurExt,SurSphere}; Value 0; }
	}}
}


FunctionSpace {
	{ Name HcurlA; Type Form1; // magnetic vector potential
		BasisFunction {
			{ Name se; NameOfCoef ae; Support #{VolAll,SurExt}; 
				Function BF_Edge; Entity EdgesOf[All]; }
			If (order >= 2)
			{ Name s3a; NameOfCoef a3a; Support #{VolAll,SurExt};
				Function BF_Edge_3F_a; Entity FacetsOf[All]; }
			{ Name s3b; NameOfCoef a3b; Support #{VolAll,SurExt};
				Function BF_Edge_3F_b; Entity FacetsOf[All]; }
			EndIf
		}
		Constraint {
			{ NameOfCoef ae; EntityType EdgesOf; NameOfConstraint CstA; }
			If (order >= 2)
			{ NameOfCoef a3a; EntityType FacetsOf; NameOfConstraint CstA; }
			{ NameOfCoef a3b; EntityType FacetsOf; NameOfConstraint CstA; }
			EndIf
			If (bound != BOUND_ABC)
				{ NameOfCoef ae; EntityType EdgesOfTreeIn; 
					EntitySubType StartingOn; NameOfConstraint CstGaugeA; }
			EndIf
		}
	}
	// Scalar Lagrange multiplier "xi" in H1_0 to enforce Coulomb gauge constraint
	{ Name H1_xi_Mag; Type Form0;
		BasisFunction {
			{ Name sn; NameOfCoef xin; Support VolCoul;
				 Function BF_Node; Entity NodesOf[All]; }
		}
		Constraint {
			{ NameOfCoef xin; EntityType NodesOf; NameOfConstraint xi_Mag; }
		}
	}
}


Formulation {
	{ Name FrmA; Type FemEquation;
		Quantity {
			{ Name a; Type Local; NameOfSpace HcurlA; }
			{ Name xi; Type Local; NameOfSpace H1_xi_Mag; }
		}
		Equation {
			Integral { [ 1 / mu[] * Dof{d a}, {d a} ];
			Integration I1; Jacobian J1; In VolAll; }
			Integral { [ -J[], {a} ];
			Integration I1; Jacobian J1; In VolSphere; }
			
			If (bound == BOUND_ABC)
				Integral{ [ 1 / (mu0*re) * Dof{a}, {a}]; 
				Integration I1; Jacobian J1; In SurExt; }

				// Coulomb Gauge
				Integral { [ Dof{a}, {d xi} ] ;
				In VolCoul; Jacobian J1; Integration I1; }
				Integral { [ Dof{d xi}, {a} ] ;
				In VolCoul; Jacobian J1; Integration I1; }
				If (order > 1)
					// Penalty factor for stabilization when higher order
					Integral { [ -1E-6 * Dof{d xi}, {d xi} ] ; 
					In VolCoul; Jacobian J1; Integration I1; }
				EndIf
			EndIf
		}
	} 
}


Resolution {
	{ Name ResMain;
		System {
			{ Name SA; NameOfFormulation FrmA; }
		}
		Operation {
			Generate[SA]; Solve[SA]; SaveSolution[SA];
		}
	}
}


PostProcessing {
	{ Name PostMain; NameOfFormulation FrmA;
		Quantity {
			// { Name J; Value {Local {
				// [ J[] ]; In VolSphere; Jacobian J1; }}
			// }
			/*{ Name B; Value {Local {
				[ Norm[{d a}] ]; In VolAll; Jacobian J1; }}
			}
			{ Name H; Value {Local {
				[ {d a}/mu[] ]; In VolAll; Jacobian J1; }}
			}*/
			/*{ Name Aex; Value {Local {
				[ {a}*0 + Aex[] ]; In #{VolSphere}; Jacobian J1; }}
			}
			{ Name A; Value {Local {
				[ {a} ]; In #{VolSphere}; Jacobian J1; }}
			}
			{ Name L2error; Value {Integral {Type Global;
				[ coef* SquNorm[Aex[]-{a}] ]; // square root done at the end in PostOperation
				Integration I2; Jacobian J1; In #{VolSphere};}}
			}*/

			{ Name Wb; Value {Integral {Type Global; 
				[ coef/(2*mu[]) * SquNorm[{d a}] ]; 
				Integration I1; Jacobian J1; In VolAll; }}
			}
			{ Name Wb2; Value {Integral {Type Global;
				[ coef/2 * J[] * {a} ]; 
				Integration I1; Jacobian J1; In VolSphere; }}
			}
		}
	}
}


PostOperation {
	{ Name PostMain; NameOfPostProcessing PostMain; 
		Format Table;
		Operation {
			Print[{prob, quarters, axis, bound, muR}, Format "Prob=%g, Quarters=%g, Axis=%g, Bound=%g, muR=%g:", File > "output.txt"]; 

			/*Print[ L2error, OnGlobal, StoreInVariable $L2error ];
			Print[ {Sqrt[$L2error]}, Format 
			" L2e = %.8g [T*m^2.5]", File > "output.txt" ];*/

			If (bound != BOUND_ABC)
				Print[ Wb, OnGlobal, StoreInVariable $Wb ];
				Print[ {$Wb, Wb[], ($Wb-Wb[])/Wb[]*10^6}, Format 
				" Wb  = %.8g [J] (analyt %.8g, %.3g ppm)", File > "output.txt" ];
			EndIf

			Print[ Wb2, OnGlobal, StoreInVariable $Wb2 ];
			Print[ {$Wb2, Wb[], ($Wb2-Wb[])/Wb[]*10^6}, Format 
			" Wb2 = %.8g [J] (analyt %.8g, %.3g ppm)", File > "output.txt" ];
		}
	} 
	{ Name PostFields; NameOfPostProcessing PostMain; 
		Operation {
			// Print[ J, OnElementsOf VolSphere, File "sphere_J.pos" ];
			Print[ Aex, OnElementsOf #{VolVacInt,VolSphere}, File "sphere_Aex.pos" ];
			Print[ A, OnElementsOf #{VolVacInt,VolSphere}, File "sphere_A.pos" ];
			// Print[ B, OnElementsOf VolAll, File "sphere_B.pos", Smoothing ];    
			// Print[ B, OnLine {{0,0,0}{rb,0,0}} {1000},
			// Format TimeTable, File "sphere_B.txt" ];
			// Print[ H, OnLine {{0,0,0}{rb,0,0}} {1000},
			// Format TimeTable, File "sphere_H.txt" ];
		}
	}
}
