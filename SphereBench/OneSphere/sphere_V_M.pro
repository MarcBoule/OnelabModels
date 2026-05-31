// How to run: see main.pro

// Conducting (V) and magnetized full sphere
// Lc2 not applicable in phi formulation since no A. Ok in A formulation, and no Lagrange mult nor mixed needed since A happens to be parallel to the surface of the sphere
// Lc2 in A formulation has high ppm though


ScalarMagPotential = 0; // 1 = scalar mag potential, 0 = vector mag potential


Group {
	VolExSphere = #{VolAll,-VolSphere};
}


Function {
	M[VolSphere] = Mp * u[]; 
	M[All] = Vector[0,0,0]; // All that are unassigned

	If (iabc == 1) 
		Call EpsIabcDiriScal;
		Call MuIabcDiriVectAndNeumScal;
	EndIf

	eps[VolSphere] = eps0 * epsR;
	mu[VolSphere]  = mu0  * muR;
	
	eps[All] = eps0; // All that are unassigned
	mu[All]  = mu0; // All that are unassigned

	// Exact results (for post analysis):
	We[] = 2*Pi * rs   * V^2    * eps0; 
	Wb[] = 4*Pi * rs^3 * Mp^2   * mu0 / (3*muR*(muR+2)); 
	Wh[] = 2*Pi * rs^3 * Mp^2   * mu0 / (3*(muR+2)); 
	Lc[] = 8*Pi * rs^3 * V * Mp * eps0 * mu0 / (3*(muR+2));
}


Constraint {
	{ Name CstV; Case { // Neumann: any SurExt4, SurXZ or SurYZ not in SurDiriV
		{ Region SurDiriV; Value 0; }
	}}
	{ Name CstPhi; Case { // Neumann: SurExt4 and any SurXZ or SurYZ not in SurDiriPhi
		{ Region SurDiriPhi; Value 0; } 
	}} 
	{ Name CstA; Case { // Neumann: any SurXZ or SurYZ not in SurDiriA
		{ Region SurDiriA; Value 0; }
	}} 
	{ Name CstGaugeA; Case {
		{ Region VolAll; SubRegion SurDiriA; Value 0; }
	}}
}


FunctionSpace {
	{ Name HgradV; Type Form0; // electric scalar potential
		BasisFunction {
			{ Name sn; NameOfCoef vn; Support #{VolExSphere,SurSphere}; 
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
	{ Name HcurlA; Type Form1; // magnetic vector potential
		BasisFunction {
			{ Name se; NameOfCoef ae; Support #{VolAll,SurSphere}; 
				Function BF_Edge; Entity EdgesOf[All]; }
			If (order >= 2)
			{ Name s3a; NameOfCoef a3a; Support #{VolAll,SurSphere};
				Function BF_Edge_3F_a; Entity FacetsOf[All]; }
			{ Name s3b; NameOfCoef a3b; Support #{VolAll,SurSphere};
				Function BF_Edge_3F_b; Entity FacetsOf[All]; }
			EndIf
		}
		Constraint {
			{ NameOfCoef ae; EntityType EdgesOf; NameOfConstraint CstA; }
			If (order >= 2)
			{ NameOfCoef a3a; EntityType FacetsOf; NameOfConstraint CstA; }
			{ NameOfCoef a3b; EntityType FacetsOf; NameOfConstraint CstA; }
			EndIf
			{ NameOfCoef ae; EntityType EdgesOfTreeIn; 
				EntitySubType StartingOn; NameOfConstraint CstGaugeA; }
		}
	}
	{ Name HgradPhi; Type Form0; // magnetic scalar potential
		BasisFunction {
			{ Name sn; NameOfCoef pn; Support VolAll; 
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
			{ Name lam; Type Local; NameOfSpace HhalfLam;}
		}
		Equation {
			Integral { [ eps[] * Dof{d v}, {d v} ]; 
			Integration I1; Jacobian J1; In VolExSphere; }
			Integral{ [ -/*eps0 */  Dof{lam}, {v}];
			Integration I1; Jacobian J1; In SurSphere; }
			Integral{ [eps0 * Dof{v}, {lam}];
			Integration I1; Jacobian J1; In SurSphere; }
			Integral{ [-eps0 * V, {lam}];
			Integration I1; Jacobian J1; In SurSphere; }
		}
	} 
	{ Name FrmPhi; Type FemEquation;
		Quantity {
			{ Name v; Type Local; NameOfSpace HgradV; }
			{ Name lam; Type Local; NameOfSpace HhalfLam;}
			{ Name p; Type Local; NameOfSpace HgradPhi; }
		}
		Equation {
			Integral { [ mu[] * Dof{d p}, {d p} ];
			Integration I1; Jacobian J1; In VolAll; }
			Integral { [ -M[]*mu0, {d p} ]; // -Br
			Integration I1; Jacobian J1; In VolSphere; }
		}
	}
	{ Name FrmA; Type FemEquation;
		Quantity {
			{ Name v; Type Local; NameOfSpace HgradV; }
			{ Name lam; Type Local; NameOfSpace HhalfLam;}
			{ Name a; Type Local; NameOfSpace HcurlA; }
		}
		Equation {
			Integral { [ 1 / mu[] * Dof{d a}, {d a} ];
			Integration I1; Jacobian J1; In VolAll; }
			Integral { [ -M[]*mu0/mu[], {d a} ]; // -Br/mu
			Integration I1; Jacobian J1; In VolSphere; }
		}
	}
}


Resolution {
	{ Name ResMain;
		System {
			{ Name SV; NameOfFormulation FrmV; }
			If (ScalarMagPotential) 
				{ Name SP; NameOfFormulation FrmPhi; }
			Else
				{ Name SA; NameOfFormulation FrmA; }
			EndIf
		}
		Operation {
			Generate[SV]; Solve[SV]; SaveSolution[SV];
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

			If (ScalarMagPotential)

			{ Name B; Value {Local {
				[ mu0*M[]-mu[]*{d p} ]; In VolAll; Jacobian J1; }}
			}
			{ Name Phi; Value {Local {
				[ {p} ]; In VolAll; Jacobian J1; }}
			}
			{ Name H; Value {Local {
				[ -{d p} ]; In VolAll; Jacobian J1; }}
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
			{ Name Lc; Value {Integral {Type Global; // needs iabc=0
				[ coef*Cross[ tr[], eps0*Cross[-{d v},mu0*M[]-mu[]*{d p}] ] ];
				Integration I1; Jacobian J1; In VolExSphere;}}
			}
			{ Name Lc2; Value {Integral {Type Global; 
			// [ coef*Cross[ r[], {lam} * {a} ] ]; // needs gauging; assumes sigma=sigma_free
				[ Vector[0,0,0] ]; // can't use above line here since no A
				Integration I1; Jacobian J1; In SurSphere;}}
			}

			Else // (ScalarMagPotential == 0)

			{ Name B; Value {Local {
				[ {d a} ]; In VolAll; Jacobian J1; }}
			}
			{ Name H; Value {Local {
				[ ({d a}-mu0*M[])/mu[] ]; In VolAll; Jacobian J1; }}
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
			{ Name Lc; Value {Integral {Type Global; // needs iabc=0
				[ coef*Cross[ tr[], eps0*Cross[-{d v},{d a}] ] ];
				Integration I1; Jacobian J1; In VolExSphere;}}
			}
			{ Name Lc2; Value {Integral {Type Global;
				[ coef*Cross[ r[], {lam} * {a} ] ]; // needs gauging; assumes sigma=sigma_free; Lagrange mult or mixed formulation not needed since vecA happens to be parallel to surface, so we can use vecA directly in surface integral
				Integration I1; Jacobian J1; In SurSphere;}}
			}

			EndIf   

		}
	}
}


PostOperation {
	{ Name PostMain; NameOfPostProcessing PostMain; 
		Format Table;
		Operation {
			Print[{prob, quarters, axis, iabc, epsR, muR, ScalarMagPotential}, Format "Prob=%g, Quarters=%g, Axis=%g, IABC=%g, epsR=%g, muR=%g, Phi=%g:", File > "output.txt"]; 

			Print[ We, OnGlobal, StoreInVariable $We ];
			Print[ {$We, We[], ($We-We[])/We[]*10^6}, Format 
			" We  = %.8g [J] (analyt %.8g, %.3g ppm)", File > "output.txt" ];

			Print[ We2, OnGlobal, StoreInVariable $We2 ];
			Print[ {$We2, We[], ($We2-We[])/We[]*10^6}, Format 
			" We2 = %.8g [J] (analyt %.8g, %.3g ppm)", File > "output.txt" ];

			Print[ Wb, OnGlobal, StoreInVariable $Wb ];
			Print[ {$Wb, Wb[], ($Wb-Wb[])/Wb[]*10^6}, Format 
			" Wb  = %.8g [J] (analyt %.8g, %.3g ppm)", File > "output.txt" ];

			Print[ Wb2, OnGlobal, StoreInVariable $Wb2 ];
			Print[ {$Wb2, Wb[], ($Wb2-Wb[])/Wb[]*10^6}, Format 
			" Wb2 = %.8g [J] (analyt %.8g, %.3g ppm)", File > "output.txt" ];

			Print[ Wh, OnGlobal, StoreInVariable $Wh ];
			Print[ {$Wh, Wh[], ($Wh-Wh[])/Wh[]*10^6}, Format 
			" Wh  = %.8g [J] (analyt %.8g, %.3g ppm)", File > "output.txt" ];

			Print[ Wh2, OnGlobal, StoreInVariable $Wh2 ];
			Print[ {$Wh2, Wh[], ($Wh2-Wh[])/Wh[]*10^6}, Format 
			" Wh2 = %.8g [J] (analyt %.8g, %.3g ppm)", File > "output.txt" ];

			If (iabc == 0)
				Print[ Lc, OnGlobal, StoreInVariable $Lc ]; 
				Print[ {CompY[$Lc], Lc[], (CompY[$Lc]-Lc[])/Lc[]*10^6}, Format 
				" LcY = %.8g [kg*m^2/s] (analyt %.8g, %.3g ppm)", File > "output.txt" ];
			Else
				Echo[ " Lcy = [needs iabc = 0]", File > "output.txt" ];
			EndIf

			If (ScalarMagPotential == 0)
				Print[ Lc2, OnGlobal, StoreInVariable $Lc2 ]; 
				Print[ {CompY[$Lc2], Lc[], (CompY[$Lc2]-Lc[])/Lc[]*10^6}, Format 
				" LcY2 = %.8g [kg*m^2/s] (analyt %.8g, %.3g ppm)", File > "output.txt" ];
			Else
				Echo[ " Lcy2 = [needs ScalarMagPotential = 0]", File > "output.txt" ];
			EndIf

			//Print[ Lc, OnGlobal, File > "output.txt" ];
		}
	} 
	{ Name PostFields; NameOfPostProcessing PostMain; 
		Operation {
			//Print[ Phi, OnElementsOf VolAll, File "sphere_Phi.pos" ];
			//Print[ B, OnElementsOf VolSphere, File "sphere_B.pos" ];
			//Print[ H, OnElementsOf VolSphere, File "sphere_H.pos" ];
		}
	}
}
