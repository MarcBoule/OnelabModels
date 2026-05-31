// Two conducting spheres


Function {
	// Maxwell stress tensor (electric components only)
	TM[] = (SquDyadicProduct[$1] - SquNorm[$1] * TensorDiag[0.5, 0.5, 0.5]) * eps0;

	// Analytical energy
	k[] = d[]/(2*rs);
	x[] = Sqrt[k[]*k[]-1];
	y[] = k[] + x[];
	argw[] = 1 / (y[]^(2*$1+1) - 1);
	We[] = 2*Pi*eps0*(2*V)^2*rs*x[] * (argw[0]+argw[1]+argw[2]+argw[3]+argw[4]+argw[5]+argw[6]+argw[7]+argw[8]+argw[9]+argw[10]);
	
	// Analytical force 
	//argf[] = k[]/(x[]*(y[]^(2*$1+1)-1)) - (2*$1+1)*y[]^(2*$1+1) / (y[]^(2*$1+1)-1)^2;
	argf[] = k[]/(x[]*(( $z=(y[]^(2*$1+1)) )-1)) - (2*$1+1)*$z / ($z-1)^2;
	F[] = Pi*eps0*(2*V)^2*(argf[0]+argf[1]+argf[2]+argf[3]+argf[4]+argf[5]+argf[6]+argf[7]+argf[8]+argf[9]+argf[10]);
}

Constraint {
	{ Name CstV; Case { 
		{ Region VolSphere1; Value V; }
		{ Region VolSphere2; Value -V; }
	}}
}

FunctionSpace {
	{ Name HgradV; Type Form0; // electric scalar potential
		BasisFunction {
			{ Name sn; NameOfCoef vn; Support #{VolAll,SurProbe}; 
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
			Integral { [ Dof{d un}, {d un} ]; // no actual DoFs here
				Integration I1; Jacobian J1; In VolProbeLayer; }
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
			{ Name We; Value {Integral {Type Global; 
				[ coef* eps[] / (2.0) * SquNorm[-{d v}] ];
				Integration I1; Jacobian J1; In VolAll;}}
			}
			{ Name We2; Value {Integral {Type Global;
				[ coef * V * (eps0*-{d v}  )* -{d un} ];
				// Q_inside = vol_int(div(D)), with same weighted integral method as used in weighted stress tensor force calculation
				// Pp = 0, and result simply doubled since energy of sphere2
				//   is identical (and it doesn't have a probe surface)
				// this calculation is Q*V/2 for sphere1, using weighted 
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
			// Print[ un, OnElementsOf VolProbeLayer, File "sphere_un.pos" ];
			// Print[ V, OnElementsOf VolAll, File "sphere_v.pos" ];

			Print[{prob, quarters, iabc}, Format "Prob=%g, Quarters=%g, IABC=%g:", File > "output.txt"]; 

			Print[ We, OnGlobal, StoreInVariable $We ];
			Print[ {$We, We[], ($We-We[])/We[]*10^6}, Format " We  = %.8g [J] (analyt %.8g, %.3g ppm)", File > "output.txt" ];

			Print[ We2, OnGlobal, StoreInVariable $We2 ];
			Print[ {$We2, We[], ($We2-We[])/We[]*10^6}, Format " We2 = %.8g [J] (analyt %.8g, %.3g ppm)", File > "output.txt" ];

			Print[ F, OnGlobal, StoreInVariable $F ]; 
			Print[ {$Fz = CompZ[$F], $FzA = F[], ($Fz-$FzA)/$FzA*10^6}, Format " Fz = %.8g [N] (analyt %.8g, %.3g ppm)", File > "output.txt" ];
		}
	} 
}
