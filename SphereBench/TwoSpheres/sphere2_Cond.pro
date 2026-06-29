// How to run: see main.pro

// Two conducting spheres


Group {
	VolExSpheres = #{VolAll,-VolSpheres};
	SurSpheres = #{SurSphere1,SurSphere2};
	DomHGradV = #{VolExSpheres,SurExt,SurSpheres};
}


Function {
	// Maxwell stress tensor (electric components only)
	TM[] = (SquDyadicProduct[$1] - SquNorm[$1] * TensorDiag[0.5, 0.5, 0.5]) * eps0;

	// Exact results (for post analysis):
	// Energy
	k[] = d[]/(2*rs);
	x[] = Sqrt[k[]*k[]-1]; // beta
	y[] = k[] + x[];       // alpha
	argw[] = 1 / (y[]^(2*$1+1) - 1); // w_i
	We[] = 2*Pi*eps0*(2*V)^2*rs*x[] * (argw[0]+argw[1]+argw[2]+argw[3]+argw[4]+argw[5]+argw[6]+argw[7]+argw[8]+argw[9]+argw[10]);
	// Force on top sphere (Sphere1)
	//argf[] = k[]/(x[]*(y[]^(2*$1+1)-1)) - (2*$1+1)*y[]^(2*$1+1) / (y[]^(2*$1+1)-1)^2;
	argf[] = k[]/(x[]*(( $z=(y[]^(2*$1+1)) )-1)) - (2*$1+1)*$z / ($z-1)^2;
	F[] = Pi*eps0*(2*V)^2*(argf[0]+argf[1]+argf[2]+argf[3]+argf[4]+argf[5]+argf[6]+argf[7]+argf[8]+argf[9]+argf[10]);
}


Constraint {
	{ Name CstV; Case { 
		{ Region SurSphere1; Value V; }
		{ Region SurSphere2; Value -V; }
	}}
}


FunctionSpace {
	{ Name HgradV; Type Form0; // electric scalar potential
		BasisFunction {
			{ Name sn; NameOfCoef vn; Support DomHGradV; 
				Function BF_Node; Entity NodesOf[All, Not SurSpheres]; }
			{ Name sf; NameOfCoef vf; Support DomHGradV;
				Function BF_GroupOfNodes; Entity GroupsOfNodesOf[ SurSpheres ]; }
		}
		GlobalQuantity {
			{ Name GlobalVs; Type AliasOf       ; NameOfCoef vf; }
			{ Name GlobalQs; Type AssociatedWith; NameOfCoef vf; }
		}
		Constraint {
			{ NameOfCoef GlobalVs; EntityType NodesOf; NameOfConstraint CstV; }
		}
	}
}


Formulation {
	{ Name FrmV; Type FemEquation;
		Quantity {
			{ Name v; Type Local; NameOfSpace HgradV; }
			{ Name un; Type Local; NameOfSpace H_Un; }
			{ Name U; Type Global; NameOfSpace HgradV[GlobalVs]; }
			{ Name Q; Type Global; NameOfSpace HgradV[GlobalQs]; }
		}
		Equation {
			Integral { [ eps[] * Dof{d v}, {d v} ]; 
			Integration I1; Jacobian J1; In VolExSpheres; }
			Integral { [ Dof{d un}, {d un} ]; // no actual DoFs here
			Integration I1; Jacobian J1; In VolProbeLayer; }
			GlobalTerm { [ -Dof{Q}, {U} ]; In SurSpheres; }
			
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
				// [ {v} ]; In VolExSpheres; Jacobian J1; }}
			// }
			{ Name We; Value {Integral {Type Global; 
				[ coef* eps[] / (2.0) * SquNorm[-{d v}] ];
				Integration I1; Jacobian J1; In VolExSpheres;}}
			}
			{ Name We21; Value { Term {
				[ coef/2 * {Q} * {U} ]; In SurSphere1;}} 
			}
			{ Name We22; Value { Term {
				[ coef/2 * {Q} * {U} ]; In SurSphere2;}} 
			}
			{ Name We3; Value {Integral {Type Global;
				[ coef * V * (eps0*-{d v}  )* -{d un} ];
				// Q_inside = vol_int(div(D)), with same weighted integral method as used in weighted stress tensor force calculation
				// Pp = 0, and result simply doubled since energy of sphere2
				//   is identical (and it doesn't have a probe surface)
				// this calculation is 2x the Q*V/2 of sphere1, using weighted 
				Integration I1; Jacobian J1; In VolProbeLayer;}}
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

			If (bound != BOUND_ABC)
				Print[ We, OnGlobal, StoreInVariable $We ];
				Print[ {$We, We[], ($We-We[])/We[]*10^6}, Format " We  = %.8g [J] (analyt %.8g, %.3g ppm)", File > "output.txt" ];
			EndIf

			Print[ We21, OnRegion SurSphere1, StoreInVariable $We21 ];
			Print[ We22, OnRegion SurSphere2, StoreInVariable $We22 ];
			Print[ {$We21+$We22, We[], ($We21+$We22-We[])/We[]*10^6}, Format " We2 = %.8g [J] (analyt %.8g, %.3g ppm)", File > "output.txt" ];

			Print[ We3, OnGlobal, StoreInVariable $We3 ];
			Print[ {$We3, We[], ($We3-We[])/We[]*10^6}, Format " We3 = %.8g [J] (analyt %.8g, %.3g ppm)", File > "output.txt" ];

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
