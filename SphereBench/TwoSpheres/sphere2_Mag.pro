// How to run: see main.pro

// Two magnetized spheres
// see constraints on u1[] and u2[] in main.pro when quarters < 4


Function {
	M[VolSphere1] = Mp * u1[];
	M[VolSphere2] = Mp * u2[];
	M[All] = Vector[0,0,0]; // All that are unassigned


	// Maxwell stress tensor (magnetic components only)
	TM[] = (SquDyadicProduct[$1] - SquNorm[$1] * TensorDiag[0.5, 0.5, 0.5]) / mu0;

	// Analytical energies
	Wb[] = 4*Pi * rs^3 * Mp^2 * mu0 / 9 * (2 + (rs/d[])^3 * (3 * (u1[]*ud[])*(u2[]*ud[])-u1[]*u2[]) );
	Wh[] = 4*Pi * rs^3 * Mp^2 * mu0 / 9 * (1 - (rs/d[])^3 * (3 * (u1[]*ud[])*(u2[]*ud[])-u1[]*u2[]) );

	// Analytical force on top sphere
	F[] = 4*Pi * rs^6 * Mp^2 * mu0 / (3 * d[]^4) * (
		ud[] * (u1[] *u2[]) +
		u1[] * (ud[] *u2[]) +
		u2[] * (ud[] *u1[]) -
		5 * ud[] * (ud[] *u1[]) * (ud[] *u2[])
	);

	// Analytical torque on top sphere
	Tau[] = 4*Pi * rs^6 * Cross[u1[], 3*(u2[]*ud[])*ud[] - u2[]] * Mp^2 * mu0 / (9 * d[]^3);
}


Constraint {
	{ Name CstPhi; Case { 
		{ Region PtRefPot;  Value 0; }// not needed, but helps to anchor the scalar potential to a reasonable value
	}}
}


FunctionSpace {
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
			{ Name un; Type Local; NameOfSpace H_Un; }
		}
		Equation {
			Integral { [ mu[] * Dof{d p}, {d p} ];
			Integration I1; Jacobian J1; In VolAll; }
			Integral { [ -M[]*mu0, {d p} ]; // -Br
			Integration I1; Jacobian J1; In VolSpheres; }
			Integral { [ Dof{d un}, {d un} ]; // no actual DoFs here
			Integration I1; Jacobian J1; In VolProbeLayer; }
			
			If (bound == BOUND_ABC)
				// 1st order ABC (n = 2 since dipole is leading harmonic)
				Integral{ [ mu0 * 2 / (re) * Dof{p}, {p}];
				Integration I1; Jacobian J1; In SurExt; }
			EndIf
		}
	}

}

Resolution {
	{ Name ResMain;
		System {
			{ Name SP; NameOfFormulation FrmPhi;}
		}
		Operation {
			Generate[SP]; Solve[SP]; SaveSolution[SP];
		}
	} 
}

PostProcessing {
	{ Name PostMain; NameOfFormulation FrmPhi;
		Quantity {
			// { Name un; Value {Local {
				// [ -{d un} ]; In VolProbeLayer; Jacobian J1; }}
			// }
			{ Name Wb; Value {Integral {Type Global; 
				[ coef/(2*mu[]) * SquNorm[mu0*M[]-mu[]*{d p}] ]; 
				Integration I1; Jacobian J1; In VolAll; }}
			}
			{ Name Wb2; Value {Integral {Type Global; 
				[ coef/2 * M[] * (mu0*M[]-mu[]*{d p}) ]; 
				Integration I1; Jacobian J1; In VolSpheres; }}
			}
			{ Name Wh; Value {Integral {Type Global; 
				[ coef* mu[]/2 * SquNorm[-{d p}] ]; 
				Integration I1; Jacobian J1; In VolAll; }}
			}
			{ Name Wh2; Value {Integral {Type Global; 
				[ coef* -mu[]/2 * M[] * (-{d p}) ]; 
				Integration I1; Jacobian J1; In VolSpheres; }}
			}
			{ Name F; Value {Integral {Type Global; 
				[ coef* TM[mu0*M[] - mu[]*{d p}] * -{d un} ];
				In VolProbeLayer; Integration I1; Jacobian J1;}}
			}
			{ Name Tau; Value {Integral {Type Global;
				[ coef* Cross[(XYZ[] - Vector[xs,ys,zs]), (TM[mu0*M[] - mu[]*{d p}] * -{d un})] ];
				In VolProbeLayer; Integration I1; Jacobian J1;}}
			} 
		}
	}
}

PostOperation {
	{ Name PostMain; NameOfPostProcessing PostMain;
		Format Table;
		Operation {
			// Print[ un, OnElementsOf VolProbeLayer, File "sphere_un.pos" ];

			Print[{prob, quarters, bound, CompX[u1[]], CompZ[u1[]], CompX[u2[]], CompZ[u2[]]}, Format "Prob=%g, Quarters=%g, Bound=%g, u1xz=%.3g,%.3g, u2xz=%.3g,%.3g:", File > "output.txt"];
				
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

			// Divisions by 0 will return inf, so no need to test
			Print[ F, OnGlobal, StoreInVariable $F ]; 
			Print[ {CompX[$F], CompX[F[]], (CompX[$F]-CompX[F[]])/CompX[F[]]*10^6}, Format " Fx = %.8g [N] (analyt %.8g, %.3g ppm)", File > "output.txt" ];
			//Print[ {CompY[$F], CompY[F[]], (CompY[$F]-CompY[F[]])/CompY[F[]]*10^6}, Format " Fy = %.8g [N] (analyt %.8g, %.3g ppm)", File > "output.txt" ];
			Print[ {CompZ[$F], CompZ[F[]], (CompZ[$F]-CompZ[F[]])/CompZ[F[]]*10^6}, Format " Fz = %.8g [N] (analyt %.8g, %.3g ppm)", File > "output.txt" ];

			Print[ Tau, OnGlobal, StoreInVariable $Tau ]; 
			Print[ {$Ty = CompY[$Tau], $TyA = CompY[Tau[]], ($Ty-$TyA)/$TyA*10^6}, Format " Ty = %.8g [N] (analyt %.8g, %.3g ppm)", File > "output.txt" ];
		}
	} 
}
