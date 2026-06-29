// How to run: see main.pro

// Two non-conducting spheres
// epsR = 1 to avoid cross polarizations


Function {
	Rho[VolSphere1] =  rho_f;
	Rho[VolSphere2] = -rho_f;
	Rho[All] = 0; // All that are unassigned	
	
	// Maxwell stress tensor (electric components only)
	TM[] = (SquDyadicProduct[$1] - SquNorm[$1] * TensorDiag[0.5, 0.5, 0.5]) * eps0;

	// Exact results (for post analysis):
	// E field (superposition)
	E1ex[VolSphere1] = rho_f*r1[]/(3*eps0);
	E1ex[All] = rho_f*rs^3/nr1[]^3*r1[]/(3*eps0);
	E2ex[VolSphere2] = -rho_f*r2[]/(3*eps0);
	E2ex[All] = -rho_f*rs^3/nr2[]^3*r2[]/(3*eps0);
	// Energy
	// Q = rho_f * 4/3*Pi*rs^3;
	// We[] = Q*Q*(6*d[]-5*rs)/(20*Pi*eps0*rs*d[]);
	We[] = 4*Pi * rs^5 * rho_f^2 * (6/5 - rs/d[]) / (9*eps0);
	// Force on top sphere (Sphere1)
	// F[] = -1/(4*Pi*eps0)*Q*Q/(d[]*d[]);
	F[] = -4*Pi * rs^6 * rho_f^2 / (9 * d[]*d[] * eps0);
}


Constraint {
	{ Name CstV; Case { 
		If (bound != BOUND_ABC)
		{ Region PtRefPot; Value 0; }// helps improve We2; value and location are arbitrary
		EndIf
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
			{ Name un; Type Local; NameOfSpace H_Un; }
		}
		Equation {
			Integral { [ eps[] * Dof{d v}, {d v} ]; 
			Integration I1; Jacobian J1; In VolAll; }
			Integral { [ -Rho[], {v} ]; 
			Integration I1; Jacobian J1; In VolSpheres; }
			Integral { [ Dof{d un}, {d un} ]; // no actual DoFs here
			Integration I1; Jacobian J1; In VolProbeLayer; }
			
			If (bound == BOUND_ABC)
				// 1st order ABC (n = 2 since dipole is leading harmonic)
				Integral{ [ eps0 * 2 / re * Dof{v}, {v}]; 
				Integration I1; Jacobian J1; In SurExt; }		
			EndIf
		}
	} 
}

Resolution {
	{ Name ResMain;
		System {
			{ Name SV; NameOfFormulation FrmV;}
		}
		Operation {
			Generate[SV]; Solve[SV]; SaveSolution[SV];
		}
	} 
}

PostProcessing {
	{ Name PostMain; NameOfFormulation FrmV;
		Quantity {
			// { Name un; Value {Local {
				// [ {un} ]; In VolProbeLayer; Jacobian J1; }}
			// }
			// { Name V; Value {Local {
				// [ {v} ]; In VolAll; Jacobian J1; }}
			// }
			{ Name L2error; Value {Integral {Type Global; 
				[ coef* SquNorm[(E1ex[]+E2ex[])-(-{d v})] ]; // square root in PostOperation
				Integration I2; Jacobian J1; In #{VolVacInt,VolSpheres};}}
			}
			{ Name E2; Value {Integral {Type Global; // E^2 exact integral
				[ coef* SquNorm[E1ex[]+E2ex[]] ];
				Integration I2; Jacobian J1; In #{VolVacInt,VolSpheres};}}
			}
			{ Name We; Value {Integral {Type Global; 
				[ coef* eps[] / (2.0) * SquNorm[-{d v}] ];
				Integration I1; Jacobian J1; In VolAll;}}
			}
			{ Name We2; Value {Integral {Type Global; 
				[ coef/2 * Rho[] * {v} ]; // rho_f for linear dielectrics (Griffiths 5ed p.198)
				Integration I1; Jacobian J1; In VolSpheres;}}
			}
			{ Name F; Value {Integral {Type Global; 
				[ coef* TM[-{d v}] * -{d un} ];
				Integration I1; Jacobian J1; In VolProbeLayer;}}
			}
		}
	}
}

PostOperation {
	{ Name PostMain; NameOfPostProcessing PostMain;
		Format Table;
		Operation {
			Print[{prob, quarters, bound, s}, Format "Prob=%g, Quarters=%g, Bound=%g, s=%g:", File > "output.txt"]; 

			Print[ L2error, OnGlobal, StoreInVariable $L2error ];
			Print[ E2, OnGlobal, StoreInVariable $E2 ];
			Print[ {Sqrt[$L2error/$E2]}, Format 
			" RelL2e = %.8g [1]", File > "output.txt" ];

			If (bound != BOUND_ABC)
				Print[ We, OnGlobal, StoreInVariable $We ];
				Print[ {$We, We[], ($We-We[])/We[]*10^6}, Format " We  = %.8g [J] (analyt %.8g, %.3g ppm)", File > "output.txt" ];
			EndIf

			Print[ We2, OnGlobal, StoreInVariable $We2 ];
			Print[ {$We2, We[], ($We2-We[])/We[]*10^6}, Format " We2 = %.8g [J] (analyt %.8g, %.3g ppm)", File > "output.txt" ];

			Print[ F, OnGlobal, StoreInVariable $F ]; 
			Print[ {$Fz = CompZ[$F], $FzA = F[], ($Fz-$FzA)/$FzA*10^6}, Format " Fz = %.8g [N] (analyt %.8g, %.3g ppm)", File > "output.txt" ];
		}
	} 
	{ Name PostFields; NameOfPostProcessing PostMain; 
		Operation {
			// Print[ un, OnElementsOf VolProbeLayer, File "sphere_un.pos" ];
			// Print[ V, OnElementsOf VolExSpheres, File "sphere_v.pos" ];
		}
	}
}
