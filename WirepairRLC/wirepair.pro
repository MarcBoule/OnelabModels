// How to run: open this file in Gmsh, press the Run button
Include "wirepair_common.pro";

// Constants
eps0  = 8.8541878188E-12;// ref permittivity (F/m)
mu0   = 1.25663706127E-6;// ref permeability (H/m)

// Simulation parameters
vIn   = 10.0;        // wire2 potential at z=0 (V)
vOut  = 9.9999;      // wire2 pot. at z=ell (V)
sigma = 1.0 / 1.7E-8;// conductivity (S/m)
epsRelInsul = 2.25;  // relative eps (unitless)
muRelCond   = 1.0;   // relative mu (unitless)

Group {
	VolCond   = #{100}; // conductor volume
	VolInsul  = #{101}; // insulator volume
	VolInt    = #{102}; // interior vacuum
	VolExt    = #{103}; // exterior vacuum (shell)
	SurLeft   = #{120}; // left boundary at x=0
	SurFront  = #{121}; // front plane at z=ell
	SurBack   = #{122}; // back plane at z=0
	SurExt    = #{123}; // curved ext. bound. at r=re
	VolNoCond = #{VolInsul, VolInt, VolExt};
	VolAll    = #{VolCond, VolNoCond};
	SurDiricA = #{SurLeft, SurExt, SurFront,SurBack};
}

Function {
	mu[VolCond]   = mu0 * muRelCond;
	mu[All]       = mu0; // all others
	eps[VolInsul] = eps0 * epsRelInsul; 
	eps[All]      = eps0; // all others
}

Constraint {
	{ Name CstVn; Case {
		{ Region SurLeft; Value 0; }
		{ Region VolCond; Value vIn; }
	}} // Neumann Vn: y=0, front, back, ext
	{ Name CstVc; Case {
		{ Region SurBack; Value vIn; }
		{ Region SurFront; Value vOut; }
	}} // Neumann Vc: lateral conductor face
		{ Name CstA;  Case {
		{ Region SurDiricA; Value 0; }
	}} // Neumann A: y=0
		{ Name CstGaugeA; Case { 
		{ Region VolAll; SubRegion SurDiricA; Value 0; }
	}}
}

FunctionSpace {
	{ Name HgradVn; Type Form0; // V nonconductor
		BasisFunction {
			{ Name snn; NameOfCoef vnn; Support VolNoCond;
				Function  BF_Node; Entity NodesOf[All];
			}//o1
		}  // o2 BF automat. added for o2 mesh in Form0
		Constraint {
			{ NameOfCoef vnn;  EntityType NodesOf; 
				NameOfConstraint CstVn; 
			}
		}
	}
	{ Name HgradVc; Type Form0; // V conductor
		BasisFunction {
			{ Name snc; NameOfCoef vnc; Support VolCond; 
				Function  BF_Node; Entity NodesOf[All];
			}//o1
		}  // o2 BF automat. added for o2 mesh in Form0
		Constraint {
			{ NameOfCoef vnc;  EntityType NodesOf; 
				NameOfConstraint CstVc; 
			}
		}
	}
	{ Name HcurlA; Type Form1; // A all
		BasisFunction {
			{ Name se;  NameOfCoef ae;  Support VolAll; 
				Function  BF_Edge; Entity EdgesOf[All];
			}//o1
			{ Name s3a; NameOfCoef a3a; Support VolAll;
				Function  BF_Edge_3F_a; // o2 (2nd order) BF 
				Entity    FacetsOf[All]; 
			}
			{ Name s3b; NameOfCoef a3b; Support VolAll;
				Function  BF_Edge_3F_b; // o2 (2nd order) BF   
				Entity    FacetsOf[All]; 
			}
		}
		Constraint {
			{ NameOfCoef ae;   EntityType EdgesOf;
				NameOfConstraint CstA; 
			}
			{ NameOfCoef a3a;  EntityType FacetsOf; 
				NameOfConstraint CstA; 
			}
			{ NameOfCoef a3b;  EntityType FacetsOf; 
				NameOfConstraint CstA; 
			}
			{ NameOfCoef ae;   EntityType EdgesOfTreeIn; 
				EntitySubType    StartingOn;
				NameOfConstraint CstGaugeA; 
			}
		}
	}
}

Jacobian {
	{ Name J1; Case { 
		{ Region VolExt; Jacobian VolCylShell{rb,re}; }
		{ Region All;    Jacobian Vol; } // all others
	}}      
}

Integration {
	{ Name I1; Case {{ Type Gauss; Case { 
		{ GeoElement Tetrahedron2; NumberOfPoints 4; }
	}}}}
}

Formulation {
	{ Name FrmVn; Type FemEquation; // V nonconductor
		Quantity {
			{ Name vn; Type Local; NameOfSpace HgradVn; }
		}
		Equation {// Electrostatics (Dirichlet source)
			Integral { [ eps[] * Dof{d vn}, {d vn} ]; 
			Integration I1; Jacobian J1; In VolNoCond; }
		}
	} 
	{ Name FrmVc; Type FemEquation; // V conductor
		Quantity {
			{ Name vc; Type Local; NameOfSpace HgradVc; }
		}
		Equation {// Electrokinetics (Dirichlet source)
			Integral { [ sigma * Dof{d vc}, {d vc} ]; 
			Integration I1; Jacobian J1; In VolCond; }
		}
	} 
	{ Name FrmA; Type FemEquation; // A all
		Quantity {// vn not used, for post-process only
			{ Name vn; Type Local; NameOfSpace HgradVn; }
			{ Name vc; Type Local; NameOfSpace HgradVc; }
			{ Name a;  Type Local; NameOfSpace HcurlA; }
		}
		Equation {// Magnetostatics (J=sigma*E source)
			Integral { [ 1 / mu[] * Dof{d a}, {d a} ];
			Integration I1; Jacobian J1; In VolAll; }
			Integral { [ sigma * {d vc}, {a} ];
			Integration I1; Jacobian J1; In VolCond; }
		}
	}
}

Resolution {
	{ Name ResAll;
		System {
			{ Name SVn; NameOfFormulation FrmVn; }
			{ Name SVc; NameOfFormulation FrmVc; }
			{ Name SA;  NameOfFormulation FrmA; }
		}
		Operation {
			Generate[SVn]; Solve[SVn]; SaveSolution[SVn];
			Generate[SVc]; Solve[SVc]; SaveSolution[SVc];
			Generate[SA];  Solve[SA];  SaveSolution[SA];
		}
	}
} 

PostProcessing {
	{ Name PostRLC; NameOfFormulation FrmA;
		Quantity {
			{ Name Pc; Value {Integral {Type Global; 
				[ 4 * sigma * SquNorm[-{d vc}] ];
				Integration I1; Jacobian J1; In VolCond; }}
			}
			{ Name Wb; Value {Integral {Type Global; 
				[ 4 / (2*mu[]) * SquNorm[{d a}] ];
				Integration I1; Jacobian J1; In VolAll; }}
			}
			{ Name We; Value {Integral {Type Global; 
				[ 4 * (eps[]/2) * SquNorm[-{d vn}] ];
				Integration I1; Jacobian J1; In VolNoCond;}}
			}
			{ Name VolC; Value {Integral {Type Global; 
				[ 4 * 1 ];
				Integration I1; Jacobian J1; In VolCond; }}
			}
		}
	}
}

PostOperation {
	{ Name PostRLC; NameOfPostProcessing PostRLC; 
		Format Table;
		Operation {
			Print[ Pc, OnGlobal, StoreInVariable $Pc ];
			Print[ Wb, OnGlobal, StoreInVariable $Wb ];     
			Print[ We, OnGlobal, StoreInVariable $We ];
			Print[{hasDefect, epsRelInsul}, Format "hasDefect=%g, epsRelInsul=%g :", File > "output.txt"];    
			Print[ {$I = $Pc / (2*(vIn-vOut))}, Format "I = %.10g [A]", File > "output.txt" ];
			Print[ {2*(vIn-vOut) / ($I*ell) *10^3}, Format "R = %.10g [mOhm/m]", File > "output.txt" ];
			Print[ {2*$Wb / ($I^2*ell) *10^9}, Format "L = %.10g [nH/m]", File > "output.txt" ];
			Print[ {2*$We / (4*vIn^2*ell) *10^12}, Format "C = %.10g [pF/m]", File > "output.txt" ];
			Print[ VolC, OnGlobal, StoreInVariable $Vol ];
			Print[ {$Vol *10^9}, Format "Vol = %.10g [mm^3]", File > "output.txt" ];
		}
	} 
}

// Extra post

PostProcessing {
	{ Name Fields; NameOfFormulation FrmA;
		Quantity {
			{ Name Vn; Value {Local {
				[ {vn} ]; In VolNoCond; Jacobian J1; }}
			}
			{ Name En; Value {Local {
				[ -{d vn} ]; In VolNoCond; Jacobian J1; }}
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
			Print[ En, OnLine {{0,0,ell/2}{rb,0,ell/2}} {1000}, Format TimeTable, File "wirepair_En.txt" ];
			Print[ B, OnElementsOf VolAll, File "wirepair_B.pos" ];
			Print[ B, OnLine {{0,0,ell/2}{rb,0,ell/2}} {1000}, Format TimeTable, File "wirepair_B.txt" ];
		}
	}
}
