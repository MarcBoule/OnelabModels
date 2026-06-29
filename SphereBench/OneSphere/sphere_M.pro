// How to run: see main.pro

// Uniformly magnetized full sphere
// Note: magnetic scalar potential is ungauged when SurDiriPhi is empty


ScalarMagPotential = 1; // 1 = scalar mag potential, 0 = vector mag potential


Group {
	SurCoul = #{SurSphere,SurDiriA};
	VolExt3Shell = ElementsOf[ VolExt3, OnNegativeSideOf SurExt ];
	VolCoul = #{VolSphere,VolExt3Shell}; // don't Coulomb gauge everything since slow
		
	If (bound == BOUND_ABC)
		SurDiriA -= #{SurExt};
	EndIf
}


Function {
	M[VolSphere] = Mp * u[]; 
	M[All] = Vector[0,0,0]; // All that are unassigned

	If (bound == BOUND_IABC) 
		Call MuIabcDiriVectAndNeumScal;
	EndIf
	mu[VolSphere]  = mu0 * muR;
	mu[All]  = mu0; // All that are unassigned

	// Exact results (for post analysis):
	Bex[VolSphere] = mu0*Mp*u[]*2/(muR+2);
	Bex[VolVacInt] = mu0*Mp*rs^3/(muR+2)*( 3/nr[]^5*(u[]*r[])*r[] - u[]/nr[]^3 );
	Hex[VolSphere] = -Mp*u[]/(muR+2);
	Hex[VolVacInt] = Bex[]/mu0;
	Wb[] = 4*Pi * rs^3 * Mp^2  * mu0 / (3*muR*(muR+2));
	Wh[] = 2*Pi * rs^3 * Mp^2  * mu0 / (3*(muR+2));
}


Constraint {
	{ Name CstPhi; Case { // Neumann: SurExt and any SurXZ or SurYZ not in SurDiriPhi
		{ Region SurDiriPhi; Value 0; } 
	}} 
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
	{ Name HgradPhi; Type Form0; // magnetic scalar potential
		BasisFunction {
			{ Name sn; NameOfCoef pn; Support #{VolAll,SurExt}; 
				Function BF_Node; Entity NodesOf[All]; }
		}
		Constraint {
			{ NameOfCoef pn; EntityType NodesOf; NameOfConstraint CstPhi; }
		}
	}
}


Formulation {
	{ Name FrmPhi; Type FemEquation;
		Quantity {
			{ Name p; Type Local; NameOfSpace HgradPhi; }
		}
		Equation {
			Integral { [ mu[] * Dof{d p}, {d p} ];
			Integration I1; Jacobian J1; In VolAll; }
			Integral { [ -M[]*mu0, {d p} ]; // -Br
			Integration I1; Jacobian J1; In VolSphere; }

			If (bound == BOUND_ABC)
				// 1st order ABC (n = 2 since dipole is leading harmonic)
				Integral{ [ mu0 * 2 / (re) * Dof{p}, {p}];
				Integration I1; Jacobian J1; In SurExt; }
			EndIf
		}
	}
	{ Name FrmA; Type FemEquation;
		Quantity {
			{ Name a; Type Local; NameOfSpace HcurlA; }
			{ Name xi; Type Local; NameOfSpace H1_xi_Mag; }
		}
		Equation {
			Integral { [ 1 / mu[] * Dof{d a}, {d a} ];
			Integration I1; Jacobian J1; In VolAll; }
			Integral { [ -M[]*mu0/mu[], {d a} ]; // -Br/mu (=-Hc)
			Integration I1; Jacobian J1; In VolSphere; }

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
			If (ScalarMagPotential) 
				{ Name SP; NameOfFormulation FrmPhi; }
			Else
				{ Name SA; NameOfFormulation FrmA; }
			EndIf
		}
		Operation {
			If (ScalarMagPotential) 
				Generate[SP]; Solve[SP]; SaveSolution[SP];
			Else
				Generate[SA]; Solve[SA]; SaveSolution[SA];
			EndIf
		}
	}
}


PostProcessing {
	If (ScalarMagPotential)
		{ Name PostMain; NameOfFormulation FrmPhi; }
	Else
		{ Name PostMain; NameOfFormulation FrmA; }
	EndIf

	{ Append; Name PostMain;
		Quantity {
			{ Name B2; Value {Integral {Type Global; // B^2 exact integral
				[ coef* SquNorm[Bex[]] ];
				Integration I2; Jacobian J1; In #{VolVacInt,VolSphere};}}
			}

			If (ScalarMagPotential)

			{ Name Phi; Value {Local {
				[ {p} ]; In VolAll; Jacobian J1; }}
			}
			// { Name B; Value {Local {
				// [ mu0*M[]-mu[]*{d p} ]; In VolAll; Jacobian J1; }}
			// }
			// { Name H; Value {Local {
				// [ -{d p} ]; In VolAll; Jacobian J1; }}
			// }
			{ Name L2error; Value {Integral {Type Global; 
				[ coef* SquNorm[Bex[]-(mu0*M[]-mu[]*{d p})] ]; // square root in PostOperation
				Integration I2; Jacobian J1; In #{VolVacInt,VolSphere};}}
			}

			{ Name Wb; Value {Integral {Type Global; 
				[ coef/(2*mu[]) * SquNorm[mu0*M[]-mu[]*{d p}] ]; 
				Integration I1; Jacobian J1; In VolAll; }}
			}
			{ Name Wb2; Value {Integral {Type Global; 
				[ coef* mu0/(2*mu[]) * M[] * (mu0*M[]-mu[]*{d p}) ]; 
				Integration I1; Jacobian J1; In VolSphere; }}
			}
			{ Name Wh; Value {Integral {Type Global; 
				[ coef*mu[]/2 * SquNorm[-{d p}] ]; 
				Integration I1; Jacobian J1; In VolAll; }}
			}
			{ Name Wh2; Value {Integral {Type Global; 
				[ coef* -mu0/2 * M[] * (-{d p}) ]; 
				Integration I1; Jacobian J1; In VolSphere;}}
			}

			Else // (ScalarMagPotential == 0)

			// { Name B; Value {Local {
				// [ {d a} ]; In VolAll; Jacobian J1; }}
			// }
			// { Name H; Value {Local {
				// [ ({d a}-mu0*M[])/mu[] ]; In VolAll; Jacobian J1; }}
			// }

			{ Name L2error; Value {Integral {Type Global; 
				[ coef* SquNorm[(Bex[]-{d a})] ]; // square root in PostOperation
				Integration I2; Jacobian J1; In #{VolVacInt,VolSphere};}}
			}
			{ Name Wb; Value {Integral {Type Global; 
				[ coef/(2*mu[]) * SquNorm[{d a}] ];
				Integration I1; Jacobian J1; In VolAll; }}
			}
			{ Name Wb2; Value {Integral {Type Global;
				[ coef* mu0/(2*mu[]) * M[] * {d a} ];
				Integration I1; Jacobian J1; In VolSphere; }}
			}
			{ Name Wh; Value {Integral {Type Global; 
				[ coef*mu[]/2 * SquNorm[({d a}-mu0*M[])/mu[] ] ]; 
				Integration I1; Jacobian J1; In VolAll; }}
			}
			{ Name Wh2; Value {Integral {Type Global; 
				[ coef* -mu0/2 * M[] * (({d a}-mu0*M[])/mu[]) ];
				Integration I1; Jacobian J1; In VolSphere;}}
			}

			EndIf   

		}
	}
}


PostOperation {
	{ Name PostMain; NameOfPostProcessing PostMain; 
		Format Table;
		Operation {
			Print[{prob, quarters, axis, bound, muR, ScalarMagPotential}, Format "Prob=%g, Quarters=%g, Axis=%g, Bound=%g, muR=%g, Phi=%g:", File > "output.txt"]; 

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

			If (bound != BOUND_ABC)
				Print[ Wh, OnGlobal, StoreInVariable $Wh ];
				Print[ {$Wh, Wh[], ($Wh-Wh[])/Wh[]*10^6}, Format 
				" Wh  = %.8g [J] (analyt %.8g, %.3g ppm)", File > "output.txt" ];
			EndIf

			Print[ Wh2, OnGlobal, StoreInVariable $Wh2 ];
			Print[ {$Wh2, Wh[], ($Wh2-Wh[])/Wh[]*10^6}, Format 
			" Wh2 = %.8g [J] (analyt %.8g, %.3g ppm)", File > "output.txt" ];
		}
	} 
	{ Name PostFields; NameOfPostProcessing PostMain; 
		Operation {
			Print[ Phi, OnElementsOf VolAll, File "sphere_Phi.pos" ];
			// Print[ B, OnElementsOf VolSphere, File "sphere_B.pos" ];
			// Print[ H, OnElementsOf VolSphere, File "sphere_H.pos" ];
		}
	}
}
