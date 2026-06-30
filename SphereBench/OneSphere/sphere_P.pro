// How to run: see main.pro

// Uniformly polarized full sphere
// Note: electric scalar potential is ungauged when SurDiriV is empty


Group { 
	SurDiriV -= #{SurExt}; // remove since Neumann V
	If (axisP == 1)
		SurDiriV += #{SurYZ};
	ElseIf (axisP == 2)
		SurDiriV += #{SurXZ};
	Else // (axisP == 3)
	EndIf
}


Function {
	P[VolSphere] = Pp * uP[];
	P[All] = Vector[0,0,0]; // All that are unassigned

	If (bound == BOUND_IABC) 
		Call EpsIabcNeumScal;
	EndIf
	eps[VolSphere] = eps0 * epsR;
	eps[All] = eps0; // All that are unassigned

	// Exact results (for post analysis):
	Dex[VolSphere] = Pp*uP[]*2/(epsR+2);
	Dex[All] = Pp*rs^3/((epsR+2)/nr[]^3)*( 3/nr[]^2*(uP[]*r[])*r[] - uP[] );
	Eex[VolSphere] = -Pp*uP[]/(eps0*(epsR+2));
	Eex[All] = Dex[]/eps0;
	We[] = 2*Pi * rs^3 * Pp^2 / (3*eps0*(epsR+2));
	Wd[] = 4*Pi * rs^3 * Pp^2 / (3*eps0*epsR*(epsR+2)); 
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
			Integral { [ -P[], {d v} ]; 
			Integration I1; Jacobian J1; In VolSphere; }
			
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
			// { Name V; Value {Local {
				// [ {v} ]; In #{VolVacInt,VolSphere}; Jacobian J1; }}
			// }
			// { Name E; Value {Local {
				// [ -{d v} ]; In VolAll; Jacobian J1; }}
			// }
			// { Name D; Value {Local {
				// [ eps[]*(-{d v}) + P[] ]; In VolAll; Jacobian J1; }}
			// }
			{ Name L2error; Value {Integral {Type Global; 
				[ coef* SquNorm[Eex[]-(-{d v})] ]; // square root in PostOperation
				Integration I2; Jacobian J1; In #{VolVacInt,VolSphere};}}
			}
			{ Name E2; Value {Integral {Type Global; // E^2 exact integral
				[ coef* SquNorm[Eex[]] ];
				Integration I2; Jacobian J1; In #{VolVacInt,VolSphere};}}
			}

			{ Name We; Value {Integral {Type Global; 
				[ coef*(eps[]/2) * SquNorm[-{d v}] ]; 
				Integration I1; Jacobian J1; In VolAll;}}
			}
			{ Name We2; Value {Integral {Type Global; 
				[ coef* -1/2 * P[] * -{d v} ];
				Integration I1; Jacobian J1; In VolSphere;}}
			}
			{ Name Wd; Value {Integral {Type Global; 
				[ coef/(eps[]*2) * SquNorm[eps[]*(-{d v}) + P[]] ]; 
				Integration I1; Jacobian J1; In VolAll;}}
			}
			{ Name Wd2; Value {Integral {Type Global; 
				[ coef/(eps[]*2) * P[] * (eps[]*(-{d v}) + P[]) ];
				Integration I1; Jacobian J1; In VolSphere;}}
			}
		}
	}
}


PostOperation {
	{ Name PostMain; NameOfPostProcessing PostMain; 
		Format Table;
		Operation {
			Print[{prob, quarters, bound, s, axisP, epsR}, Format "Prob=%g, Quarters=%g, Bound=%g, s=%g, Axis=%g, epsR=%g:", File > "output.txt"]; 

			Print[ L2error, OnGlobal, StoreInVariable $L2error ];
			Print[ E2, OnGlobal, StoreInVariable $E2 ];
			Print[ {Sqrt[$L2error/$E2]}, Format 
			" RelL2e = %.8g [1]", File > "output.txt" ];

			If (bound != BOUND_ABC)
				Print[ We, OnGlobal, StoreInVariable $We ];
				Print[ {$We, We[], ($We-We[])/We[]*10^6}, Format 
				" We  = %.8g [J] (analyt %.8g, %.3g ppm)", File > "output.txt" ];
			EndIf

			Print[ We2, OnGlobal, StoreInVariable $We2 ];
			Print[ {$We2, We[], ($We2-We[])/We[]*10^6}, Format 
			" We2 = %.8g [J] (analyt %.8g, %.3g ppm)", File > "output.txt" ];

			If (bound != BOUND_ABC)
				Print[ Wd, OnGlobal, StoreInVariable $Wd ];
				Print[ {$Wd, Wd[], ($Wd-Wd[])/Wd[]*10^6}, Format 
				" Wd  = %.8g [J] (analyt %.8g, %.3g ppm)", File > "output.txt" ];
			EndIf

			Print[ Wd2, OnGlobal, StoreInVariable $Wd2 ];
			Print[ {$Wd2, Wd[], ($Wd2-Wd[])/Wd[]*10^6}, Format 
			" Wd2 = %.8g [J] (analyt %.8g, %.3g ppm)", File > "output.txt" ];
		}
	} 
	{ Name PostFields; NameOfPostProcessing PostMain; 
		Operation {
			// Print[ V, OnElementsOf #{VolVacInt,VolSphere}, File "sphere_V.pos" ];
			// Print[ dV, OnElementsOf #{VolVacInt,VolSphere}, File "sphere_dV.pos" ];
			// Print[ E, OnPoint {xs+rs*1.5,ys,zs},
			// Format TimeTable, File > "sphere_pts.txt" ];
			// Print[ D, OnPoint {xs+rs*1.5,ys,zs},
			// Format TimeTable, File > "sphere_pts.txt" ];
		}
	}
}
