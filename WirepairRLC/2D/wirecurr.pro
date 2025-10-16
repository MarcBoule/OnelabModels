// How to run: open this file in Gmsh, press the Run button
Include "wirecurr_common.pro";

// Simulation parameters
vIn   = 10.0;         // wire2 potential top (V)
vOut  = 9.9999;       // wire2 potential bot (V)
sigma = 1.0 / 1.7E-8; // conductivity (S/m)

Group {
	VolCond  = #{100}; // conductor volume
	SurTop   = #{120}; // top of wire
	SurBot   = #{121}; // bot of wire
}

Function {
}

Constraint {
	{ Name CstV; Case {
		{ Region SurTop;   Type Assign; Value vIn; }
		{ Region SurBot;   Type Assign; Value vOut; }
	}} // Neumann V: Left and right sides
}

FunctionSpace {
	{ Name HgradV; Type Form0; // V conductor
		BasisFunction {
			{ Name sn;  NameOfCoef vn;  Support VolCond;
				Function  BF_Node; // 1st order Basis Func.    
				Entity    NodesOf[All]; 
			}
		}  // o2 BF automat. added for o2 mesh in Form0
		Constraint {
			{ NameOfCoef vn;  EntityType NodesOf; 
				NameOfConstraint CstV; 
			}
		}
	}
}

Jacobian {
	{ Name J1; Case { 
		{ Region All; Jacobian VolAxi; }
	}}      
}

Integration {
	{ Name I1; Case {{ Type Gauss; Case { 
		{ GeoElement Triangle;  NumberOfPoints 1; }
		{ GeoElement Triangle2; NumberOfPoints 3; }
	}}}}
}

Formulation {
	{ Name FrmV; Type FemEquation; // V conductor
		Quantity {
			{ Name v; Type Local; NameOfSpace HgradV; }
		}
		Equation {// Electrokinetics (Dirichlet source)
			Integral { [ sigma * Dof{d v}, {d v} ]; 
			Integration I1; Jacobian J1; In VolCond; }
		}
	} 
}

Resolution {
	{ Name ResV;
		System {
			{ Name SV; NameOfFormulation FrmV;}
		}
		Operation {
			Generate[SV]; Solve[SV]; SaveSolution[SV];
		}
	}
} 


PostProcessing {
	{ Name PostR; NameOfFormulation FrmV;
		Quantity {
			{ Name Pc; Value {Integral {Type Global; 
				[ 2 * 2*Pi * sigma * SquNorm[-{d v}] ];
				Integration I1; Jacobian J1; In VolCond; }}
			}
			{ Name VolC; Value {Integral {Type Global; 
				[ 2 * 2*Pi * 1 ];
				Integration I1; Jacobian J1; In VolCond; }}
			}
		}
	}
}
PostOperation {
	{ Name PostR; NameOfPostProcessing PostR; 
		Format Table;
		Operation {
			Print[ Pc, OnGlobal, StoreInVariable $Pc ];
			Print[{hasDefect}, Format "hasDefect=%g", File > "output.txt"];
			Print[ {$I = $Pc / (2*(vIn-vOut))}, Format "I = %.10g [A]", File > "output.txt" ];
			Print[ {2*(vIn-vOut) / ($I*ell) *10^3}, Format "R = %.10g [mOhm/m] (two wires)", File > "output.txt" ];
			Print[ VolC, OnGlobal, StoreInVariable $Vol ];
			Print[ {$Vol *10^9}, Format "Vol = %.10g [mm^3] (two wires)", File > "output.txt" ];
		}
	} 
}

// Extra post (fields)

PostProcessing {
	{ Name PostVf; NameOfFormulation FrmV;
		Quantity {
			{ Name V; Value {Local {
				[ {v} ]; In VolCond; Jacobian J1; }}
			}
			{ Name J; Value {Local {
				[ sigma*-{d v} ]; In VolCond; Jacobian J1; }}
			}
		}
	}
}
PostOperation {
	{ Name PostVf; NameOfPostProcessing PostVf; 
		Operation {
			Print[ V, OnElementsOf VolCond, File "wirecurr_V.pos" ];
			Print[ J, OnElementsOf VolCond, File "wirecurr_J.pos" ];
		}
	} 
}



