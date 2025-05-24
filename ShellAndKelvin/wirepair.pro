Include "wirepair_dat.pro";

eps0  = 8.8541878188E-12; // ref. permittivity (F/m)
vCond = 10.0;             // wire potential (V)

Capacitance = DefineNumber[Pi * eps0 / Log[a/rc + Sqrt[(a/rc)^2-1]] * 1E12, Name "Output/Cap (Analyt.)", Units "pF/m", ReadOnly 1];
Coef = Flag_3D == 1 ? 4 : 4 * ell;


Group {
	VolInt    = #100;
	VolExt    = #101; 
	SurCond   = #120;
	SurLefts  = #121;
	SurBndInt = #122; // unused in shell version
	SurBndExt = #123;
	VolAll    = #{VolInt,VolExt};
}

Constraint {
	{ Name CstV; Case {
		If(Flag_Shell == 1)
			{ Region SurBndExt; Type Assign; Value 0.0; }
		Else
			{ Region    SurBndExt; Type Link; // Kelvin link
			  RegionRef SurBndInt;  Coefficient 1.0;
			  Function  Vector[ X[] - xe, Y[], Z[] ]; }
		EndIf
		{ Region SurCond; Type Assign; Value vCond; }
		{ Region SurLefts; Type Assign; Value 0.0; }
	}}
}

FunctionSpace {
	{ Name HgradV; Type Form0;
		BasisFunction {
			{ Name sn; NameOfCoef vn; Function BF_Node;
			  Support VolAll; Entity NodesOf[All]; }
		}
		Constraint {
			{ NameOfCoef vn; EntityType NodesOf; 
			  NameOfConstraint CstV; }
		}
	}
}

Jacobian {
	{ Name J1; Case { 
		{ Region VolInt; Jacobian Vol; }
		If (Flag_3D == 1)
			If(Flag_Shell == 1)
				{ Region VolExt; Jacobian VolCylShell{rb,re}; }
			Else // Kelvin
				{ Region VolExt; Jacobian VolCylShell{rb,0, 3, xe, 0, 0}; }
			EndIf
		Else // 2D
			If(Flag_Shell == 1)
				{ Region VolExt; Jacobian VolSphShell{rb,re}; }
			Else // Kelvin
				{ Region VolExt; Jacobian VolSphShell{rb,0, xe, 0, 0}; }
			EndIf
		EndIf
	}}
}

Integration {
	{ Name I1; Case {{ Type Gauss; Case { 
		{ GeoElement Triangle; 
		  NumberOfPoints (Flag_MinPts == 1 ? 1 : 3); }
		{ GeoElement Triangle2; 
		  NumberOfPoints (Flag_MinPts == 1 ? 3 : 6); }
		{ GeoElement Tetrahedron; 
		  NumberOfPoints (Flag_MinPts == 1 ? 1 : 4); }
		{ GeoElement Tetrahedron2; 
		  NumberOfPoints (Flag_MinPts == 1 ? 4 : 15); }
	}}}}
}

Formulation {
	{ Name FrmV; Type FemEquation;
		Quantity {
			{ Name v; Type Local; NameOfSpace HgradV; }
		}
		Equation { 
			Integral {[ eps0 * Dof{d v}, {d v} ];
			In VolAll; Jacobian J1; Integration I1; }
		}
	}
}

Resolution {
	{ Name ResMain; 
		System {      
			{ Name Sy1; NameOfFormulation FrmV; }
		}
		Operation {
			Generate[Sy1]; Solve[Sy1]; SaveSolution[Sy1];
		}
	}
} 

PostProcessing {
	{ Name PostV; NameOfFormulation FrmV;
		Quantity {
			{ Name v ; Value { Local { [ {v} ]; In VolAll ; Jacobian J1; } } }
			{ Name Cap; Value {Integral {Type Global; 
			  [ Coef * 1E12 * eps0*SquNorm[-{d v}] / (ell*(2*vCond)^2) ];
			  In VolAll; Integration I1; Jacobian J1;}}
			}
		}
	}
}

e = rb/1000;
PostOperation {
	If(Flag_Shell == 1)
		outstr = "Output/Cap (Shell)";
	Else
		outstr = "Output/Cap (Kelvin)";
	EndIf
	{ Name PostV; NameOfPostProcessing PostV;
		Format Table;	
		Operation {
			Print[ Cap, OnGlobal, SendToServer outstr, Units "pF/m", StoreInVariable $Cap];
			Print[{Flag_3D, Flag_Shell, Flag_2ndOrder, Flag_MinPts, $Cap, ($Cap - Capacitance)/Capacitance * 100}, Format "3D=%g, Shell=%g, 2ndOrd=%g, MinPts=%g, Cap=%gpF/m, Diff=%gpercent", File > "wirepair_out.txt"];
		}
	} 
	{ Name PostV2; NameOfPostProcessing PostV;
		Operation {
			Print[ v, OnElementsOf VolAll, File "wirepair_v.pos" ];
		}
	} 
}