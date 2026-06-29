// How to run: see main.pro

// Spinning uniformly surface-charged nonconducting full sphere
// epsilon_r = mu_r = 1 (or else there will be spinning bound charges and currents which would complicate things)


Group { 
	SurCoul = #{SurSphere,SurDiriA};
	VolExt3Shell = ElementsOf[ VolExt3, OnNegativeSideOf SurExt ];
	VolCoul = #{VolSphere,VolExt3Shell}; // don't Coulomb gauge everything since slow
	
	If (bound == BOUND_ABC)
		SurDiriA -= #{SurExt};
	EndIf
}


Function {
	v[] = Cross[ omega*u[], r[] ];
	K[] = sigma_f * v[];

	If (bound == BOUND_IABC) 
		Call MuIabcDiriVectAndNeumScal;
	EndIf
	mu[All]  = mu0; // All that are unassigned

	// Exact results (for post analysis):
	Bex[VolSphere] = mu0*sigma_f*omega*rs*u[]*2/3;
	Bex[All] = mu0*sigma_f*omega*rs^4/3*( 3/nr[]^5*(u[]*r[])*r[] - u[]/nr[]^3 );
	Wb[] = 4*Pi  * rs^5 * sigma_f^2 * omega^2 * (mu0/9); 
}


Constraint {
	{ Name CstA; Case { // Neumann: any SurXZ or SurYZ not in SurDiriA
		{ Region SurDiriA; Value 0; }
	}} 
	{ Name CstGaugeA; Case {
		{ Region VolAll; SubRegion SurDiriA; Value 0; }
	}}
	// Boundary condition for the Coulomb gauge multiplier "xi" (only used when
	// "bound == BOUND_ABC"):
	{ Name xi_Mag; Case {
		{ Region SurCoul; Value 0; }
	}}
}


FunctionSpace {
	{ Name HcurlA; Type Form1; // magnetic vector potential
		BasisFunction {
			{ Name se; NameOfCoef ae; Support #{VolAll,SurSphere,SurExt}; 
				Function BF_Edge; Entity EdgesOf[All]; }
			If (order >= 2)
			{ Name s3a; NameOfCoef a3a; Support #{VolAll,SurSphere,SurExt};
				Function BF_Edge_3F_a; Entity FacetsOf[All]; }
			{ Name s3b; NameOfCoef a3b; Support #{VolAll,SurSphere,SurExt};
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
			Integral { [ -K[], {a} ]; // this works because when we add the weak equations for each region considered separately (VolSphere and VolExSphere), which are separated by SurSphere, we can consider SurSphere as a Neumann term for each side; when added, the two Neumann terms on SurSphere end up representing precisely the matching condition "n x (H1-H2) = K"
			Integration I1; Jacobian J1; In SurSphere; }

			If (bound == BOUND_ABC)
				// 1st order ABC (n = 1 for vect potential when dipole is leading harmonic)
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
			// { Name B; Value {Local {
				// [ {d a} ]; In VolAll; Jacobian J1; }}
			// }

			{ Name L2error; Value {Integral {Type Global; 
				[ coef* SquNorm[Bex[]-{d a}] ]; // square root in PostOperation
				Integration I2; Jacobian J1; In #{VolVacInt,VolSphere};}}
			}
			{ Name B2; Value {Integral {Type Global; // B^2 exact integral
				[ coef* SquNorm[Bex[]] ];
				Integration I2; Jacobian J1; In #{VolVacInt,VolSphere};}}
			}

			{ Name Wb; Value {Integral {Type Global; 
				[ coef/(2*mu[]) * SquNorm[{d a}] ];
				Integration I1; Jacobian J1; In VolAll; }}
			}
			{ Name Wb2; Value {Integral {Type Global; 
				[ coef/2 * K[] * {a} ]; 
				Integration I1; Jacobian J1; In SurSphere;}}
			}
		}
	}
}


PostOperation {
	{ Name PostMain; NameOfPostProcessing PostMain; 
		Format Table;
		Operation {
			Print[{prob, quarters, bound, s, axis}, Format "Prob=%g, Quarters=%g, Bound=%g, s=%g, Axis=%g:", File > "output.txt"]; 

			Print[ L2error, OnGlobal, StoreInVariable $L2error ];
			Print[ B2, OnGlobal, StoreInVariable $B2 ];
			Print[ {Sqrt[$L2error/$B2]}, Format 
			" RelL2e = %.8g [1]", File > "output.txt" ];

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
			// Print[ B, OnElementsOf VolAll, File "sphere_B.pos" ];
		}
	}
}
