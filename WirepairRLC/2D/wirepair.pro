// How to run: open this file in Gmsh, press the Run button
Include "wirepair_common.pro";

// Constants
eps0  = 8.8541878188E-12; // ref permittivity (F/m)
mu0   = 1.25663706127E-6; // ref permeability (H/m)

// Simulation parameters
vIn   = 10.0;         // wire2 potential at z=0 (V)
I = 0.4790189599;     // current (A), from wirecurr sim
Jz = I / (Pi*rc^2);   // current density (A/m^2)
epsRelInsul = 2.25;   // relative eps (unitless)
muRelCond   = 1.0;    // relative mu (unitless)

Group {
	VolCond   = #{100}; // conductor volume
	VolInsul  = #{101}; // insulator volume
	VolVacInt = #{102}; // interior vacuum
	VolExt    = #{103}; // exterior vacuum
	SurLeft   = #{120}; // left boundary
	SurExt    = #{121}; // exterior curved boundary
	VolNoCond = #{VolInsul, VolVacInt, VolExt};
	VolAll    = #{VolCond, VolNoCond};
}

Function {
	mu[VolCond]   = mu0 * muRelCond;
	mu[All]       = mu0; // all others
	eps[VolInsul] = eps0 * epsRelInsul; 
	eps[All]      = eps0; // all others
}

Constraint {
	{ Name CstV; Case {
		{ Region SurLeft; Value 0; }
		{ Region VolCond; Value vIn; }
	}} // Neumann V: y=0, ext
	{ Name CstA; Case {
		{ Region #{SurLeft, SurExt}; Value 0; }
	}} // Neumann A: y=0
}

FunctionSpace {
	{ Name HgradV; Type Form0; // V nonconductor
		BasisFunction {
			{ Name sn; NameOfCoef vn; Support VolNoCond;
				Function BF_Node; Entity NodesOf[All]; 
			}//o1
		}  // o2 BF automat. added for o2 mesh in Form0
		Constraint {
			{ NameOfCoef vn; EntityType NodesOf; 
				NameOfConstraint CstV; 
			}
		}
	}
	{ Name HcurlA; Type Form1P; // A all
		BasisFunction {
			{ Name sn;  NameOfCoef an;  Support VolAll;
				Function  BF_PerpendicularEdge;// 1st ord. BF         
				Entity    NodesOf[All]; 
			}
		}  // 2nd order BF automatic for 2nd order mesh for Form1P
		Constraint {
			{ NameOfCoef an;  EntityType NodesOf;
				NameOfConstraint CstA; 
			}
		}
	}
}

Jacobian {
	{ Name J1; Case { 
		{ Region VolExt; Jacobian VolCylShell{rb,re}; }
		{ Region All;    Jacobian Vol; }
	}}      
}

Integration {
	{ Name I1; Case {{ Type Gauss; Case { 
		{ GeoElement Triangle2; NumberOfPoints 4; }
	}}}}
}

Formulation {
	{ Name FrmV; Type FemEquation; // V nonconductor
		Quantity {
			{ Name v; Type Local; NameOfSpace HgradV; }
		}
		Equation {// Electrostatics (Dirichlet source)
			Integral { [ eps[] * Dof{d v}, {d v} ]; 
			Integration I1; Jacobian J1; In VolNoCond; }
		}
	} 
	{ Name FrmA; Type FemEquation; // A all
		Quantity {// v not used, for post-process only
			{ Name v; Type Local; NameOfSpace HgradV; }
			{ Name a;  Type Local; NameOfSpace HcurlA; }
		}
		Equation {// Magnetostatics (J source)
			Integral { [ 1 / mu[] * Dof{d a}, {d a} ];
			Integration I1; Jacobian J1; In VolAll; }
			Integral { [ -Vector[0,0,Jz], {a} ];
			Integration I1; Jacobian J1; In VolCond; }
		}
	}
}

Resolution {
	{ Name ResAll;
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
	{ Name PostLC; NameOfFormulation FrmA;
		Quantity {
			{ Name Wb; Value {Integral {Type Global; 
				[ 4/(2*mu[]) * SquNorm[{d a}] ];
				Integration I1; Jacobian J1; In VolAll; }}
			}
			{ Name We; Value {Integral {Type Global; 
				[ 4*(eps[]/2) * SquNorm[-{d v}] ];
				Integration I1; Jacobian J1; In VolNoCond;}}
			}
		}
	}
}

PostOperation {
	{ Name PostLC; NameOfPostProcessing PostLC; 
		Format Table;
		Operation {
			Print[ Wb, OnGlobal, StoreInVariable $Wb ];     
			Print[ We, OnGlobal, StoreInVariable $We ];
			Print[{epsRelInsul}, Format "epsRelInsul=%g :", File > "output.txt"];
			Print[ { 2*$Wb / (I^2*1) *10^9 }, Format "L = %.10g [nH/m]", File > "output.txt" ];
			Print[ { 2*$We / (4*vIn^2*1) *10^12 }, Format "C = %.10g [pF/m]", File > "output.txt" ];
		}
	} 
}

// Extra post

PostProcessing {
	{ Name Fields; NameOfFormulation FrmA;
		Quantity {
			{ Name Vn; Value {Local {
				[ {v} ]; In VolNoCond; Jacobian J1; }}
			}
			{ Name En; Value {Local {
				[ -{d v} ]; In VolNoCond; Jacobian J1; }}
			}
			{ Name B; Value {Local {
				[ {d a} ]; In VolAll; Jacobian J1; }}
			}
		}
	}
}

PostOperation { 
	{ Name Fields; NameOfPostProcessing Fields; 
		Operation {
			Print[ Vn, OnElementsOf VolNoCond, File "wirepair_Vn.pos" ];
			Print[ En, OnElementsOf VolNoCond, File "wirepair_En.pos" ];
			Print[ B, OnElementsOf VolAll, File "wirepair_B.pos" ];
		}
	}
}
