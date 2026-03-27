// How to run: see main.pro

// Spinning uniformly surfacically charged nonconducting full sphere
// epsilon_r = mu_r = 1 only for mag and momentum (or else there will be spinning bound charges and currents which would complexify things)


Group { 
}


Function {
	v[] = Cross[ omega*u[], r[] ];
	K[] = sigma * v[];

	Call EpsDirichletScalar;

	// Exact results (for post analysis):
	We[] = 2*Pi  * rs^3 * sigma^2    / eps0; 
	Wb[] = 4*Pi  * rs^5 * sigma^2 * omega^2 * (mu0/9); 
	Lc[] = 8*Pi  * rs^5 * sigma^2 * omega   * (mu0/9);
}


Constraint {
	{ Name CstV; Case { // Neumann: any SurXZ, SurYZ or SurExt4 not in SurDiriV
		{ Region SurDiriV; Value 0; }
	}}
	{ Name CstAout; Case { // Neumann: any SurXZ or SurYZ not in SurDiriA
		{ Region SurDiriA; Value 0; }
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
			{ Name sn; NameOfCoef vn; Support #{VolAll,SurSphere}; 
				Function BF_Node; Entity NodesOf[All]; }
		}
		Constraint {
			{ NameOfCoef vn; EntityType NodesOf; NameOfConstraint CstV; }
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
}


Formulation {
	{ Name FrmV; Type FemEquation;
		Quantity {
			{ Name v; Type Local; NameOfSpace HgradV; }
		}
		Equation {
			Integral { [ eps[] * Dof{d v}, {d v} ]; 
			Integration I1; Jacobian J1; In VolAll; }
			Integral{ [ -sigma, {v}]; // this works because when we add the weak equations for each region considered separately (VolSphere and VolExSphere), which are separated by SurSphere, we can consider SurSphere as a Neumann term for each side; when added, the two Neumann terms on SurSphere end up representing precisely the matching condition "n dot (D1-D2) = sigma"
			Integration I1; Jacobian J1; In SurSphere; }
		}
	} 
	{ Name FrmA; Type FemEquation;
		Quantity {
			{ Name v; Type Local; NameOfSpace HgradV; }
			{ Name a; Type Local; NameOfSpace HcurlA; }
		}
		Equation {
			Integral { [ 1 / mu[] * Dof{d a}, {d a} ];
			Integration I1; Jacobian J1; In VolAll; }
			Integral { [ -K[], {a} ]; // this works because when we add the weak equations for each region considered separately (VolSphere and VolExSphere), which are separated by SurSphere, we can consider SurSphere as a Neumann term for each side; when added, the two Neumann terms on SurSphere end up representing precisely the matching condition "n x (H1-H2) = K"
			Integration I1; Jacobian J1; In SurSphere; }
		}
	}
}


Resolution {
	{ Name ResMain;
		System {
			{ Name SV; NameOfFormulation FrmV; }
			{ Name SA; NameOfFormulation FrmA; }
		}
		Operation {
			Generate[SV]; Solve[SV]; SaveSolution[SV];
			Generate[SA]; Solve[SA]; SaveSolution[SA];
		}
	}
}


PostProcessing {
	{ Name PostMain; NameOfFormulation FrmA;
		Quantity {
			{ Name E; Value {Local {
				[ -{d v} ]; In VolAll; Jacobian J1; }}
			}
			{ Name B; Value {Local {
				[ {d a} ]; In VolAll; Jacobian J1; }}
			}

			{ Name We; Value {Integral {Type Global; 
				[ coef*(eps[]/2) * SquNorm[-{d v}] ];
				Integration I1; Jacobian J1; In VolAll;}}
			}
			{ Name We2; Value {Integral {Type Global; 
				[ coef/2 * sigma * {v} ]; 
				Integration I1; Jacobian J1; In SurSphere;}}
			}
			{ Name Wb; Value {Integral {Type Global; 
				[ coef/(2*mu[]) * SquNorm[{d a}] ];
				Integration I1; Jacobian J1; In VolAll; }}
			}
			{ Name Wb2; Value {Integral {Type Global; 
				[ coef/2 * K[] * {a} ]; 
				Integration I1; Jacobian J1; In SurSphere;}}
			}
			{ Name Lc; Value {Integral {Type Global; // needs iabc=0
				[ coef*Cross[ tr[], eps0*Cross[-{d v},{d a}] ] ];
				Integration I1; Jacobian J1; In VolAll;}}
			}
			{ Name Lc2; Value {Integral {Type Global;
				[ coef*Cross[ r[], sigma * {a} ] ]; // needs gauging. Lagrange mult or mixed  formulation not needed since vecA happens to be parallel to surface, so we can use vecA directly in surface integral
				Integration I1; Jacobian J1; In SurSphere;}}
			}
		}
	}
}


PostOperation {
	{ Name PostMain; NameOfPostProcessing PostMain; 
		Format Table;
		Operation {
			Print[{prob, quarters, axis, iabc, epsR, muR}, Format "Prob=%g, Quarters=%g, Axis=%g, IABC=%g, epsR=%g, muR=%g:", File > "output.txt"]; 

			Print[ We, OnGlobal, StoreInVariable $We ];
			Print[ {$We, We[], ($We-We[])/We[]*10^6}, Format 
			" We  = %.8g [J] (analyt %.8g, %.3g ppm)", File > "output.txt" ];

			Print[ We2, OnGlobal, StoreInVariable $We2 ];
			Print[ {$We2, We[], ($We2-We[])/We[]*10^6}, Format 
			" We2 = %.8g [J] (analyt %.8g, %.3g ppm)", File > "output.txt" ];

			If (epsR != 1 || muR != 1)
				Echo[ " [needs epsR = 1 and muR = 1 for other results]", File > "output.txt" ];
			Else
				
				Print[ Wb, OnGlobal, StoreInVariable $Wb ];
				Print[ {$Wb, Wb[], ($Wb-Wb[])/Wb[]*10^6}, Format 
				" Wb  = %.8g [J] (analyt %.8g, %.3g ppm)", File > "output.txt" ];

				Print[ Wb2, OnGlobal, StoreInVariable $Wb2 ];
				Print[ {$Wb2, Wb[], ($Wb2-Wb[])/Wb[]*10^6}, Format 
				" Wb2 = %.8g [J] (analyt %.8g, %.3g ppm)", File > "output.txt" ];

				If (iabc == 0)
				Print[ Lc, OnGlobal, StoreInVariable $Lc ]; 
				Print[ {CompY[$Lc], Lc[], (CompY[$Lc]-Lc[])/Lc[]*10^6}, Format 
				" LcY = %.8g [kg*m^2/s] (analyt %.8g, %.3g ppm)", File > "output.txt" ];
				Else
				Echo[ " LcY = [needs iabc = 0]", File > "output.txt" ];
				EndIf

				Print[ Lc2, OnGlobal, StoreInVariable $Lc2 ]; 
				Print[ {CompY[$Lc2], Lc[], (CompY[$Lc2]-Lc[])/Lc[]*10^6}, Format 
				" LcY2 = %.8g [kg*m^2/s] (analyt %.8g, %.3g ppm)", File > "output.txt" ];
				
			EndIf
		}
	} 
	{ Name PostFields; NameOfPostProcessing PostMain; 
		Operation {
			// Print[ V, OnElementsOf VolExSphere, File "sphere_V.pos" ];
			//Print[ V, OnLine {{rs,1e-6,0}{rb,1e-6,0}} {1000}, 
			// Format TimeTable, File "sphere_Vs.txt" ];
			//Print[ E, OnElementsOf VolAll, File "sphere_E.pos" ];
			//Print[ Bout, OnElementsOf VolExSphere, File "sphere_Bout.pos" ];    
			//Print[ B,  OnElementsOf VolAll,   File "sphere_B.pos" ];    
		}
	}
}
