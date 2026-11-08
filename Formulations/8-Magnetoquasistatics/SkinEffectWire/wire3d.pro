Include "wire_common.pro";

// Parameters
sigma = 58e6;      // Conductivity [S/m)]
freq = 1e4;        // Frequency [1/s]
omega = 2*Pi*freq; // Angular frequency [1/s]
currMag = 5;       // Current phasor magnitude [A]
currPhase = Pi/3;  // Current phasor phase [rad]

// Constants
mu0  = 1.25663706127e-6; // Reference permeability [H/m]


Group {
	VolWire      = #{100}; 
	VolVac       = #{101};
	SurExt       = #{120};
	SurWireBack  = #{121};
	SurWireFront = #{122};
	SurVacFrontBack = #{123};

	VolAll       = #{VolWire, VolVac};
	SurDiriA     = #{SurExt,SurWireBack,SurWireFront,SurVacFrontBack};
}


Function {
	currCplx[] = currMag*Complex[Cos[currPhase],Sin[currPhase]];
	
	// Exact solution below is based on section 14.4 in 
	//   "Introductory Engineering Electromagnetics", by
	//   Branko D. Popovic, Addison-Wesley, 1971

	// Utility values
	K[] = Complex[1,-1]*Sqrt[omega*sigma*mu0/2];
	R[] = Norm[Vector[X[],Y[],0]]; 
	
	// Exact current density as a function of radial position R[]
	Jex[] = K[]*currCplx[] * JnComplex[0,K[]*R[]] / (2*Pi*r* JnComplex[1,K[]*r]);
}


Constraint {
	{ Name CstA; Case {
		{ Region SurDiriA; Value 0; }
	}}
	{ Name CstV; Case {
		{ Region SurWireFront; Value 0; }
	}}
	{ Name CstI; Case {
		{ Region SurWireBack; Value 1/coef * currMag*Complex[Cos[currPhase],Sin[currPhase]]; }
	}}
	{ Name CstGaugeA; Case {
		{ Region VolAll; SubRegion SurDiriA; Value 0; }
	}}
}


FunctionSpace {
	{ Name Hcurl; Type Form1;
		BasisFunction {
			{ Name se;  NameOfCoef ae;  Function BF_Edge;    
				Support VolAll; Entity EdgesOf[All]; }
			If (order > 1)
			{ Name s3a; NameOfCoef a3a; Function BF_Edge_3F_a;
				Support VolAll; Entity FacetsOf[All]; }
			{ Name s3b; NameOfCoef a3b; Function BF_Edge_3F_b;
				Support VolAll; Entity FacetsOf[All]; }
			EndIf
		}
		Constraint {
			{ NameOfCoef ae;  EntityType EdgesOf; NameOfConstraint CstA; }
			If (order > 1)
			{ NameOfCoef a3a; EntityType FacetsOf; NameOfConstraint CstA; }
			{ NameOfCoef a3b; EntityType FacetsOf; NameOfConstraint CstA; }
			EndIf
			{ NameOfCoef ae; EntityType EdgesOfTreeIn; 
				EntitySubType StartingOn; NameOfConstraint CstGaugeA; }
		}
	}
	{ Name Hgrad; Type Form0;
		BasisFunction {
			{ Name sn;  NameOfCoef vn; Function BF_Node;    
				Support VolWire; Entity NodesOf[All, Not SurWireBack]; }
			{ Name sr; NameOfCoef vr; Function BF_GroupOfNodes;
				Support VolWire; Entity GroupsOfNodesOf[SurWireBack]; }
		}
		GlobalQuantity {
			{ Name Voltage; Type AliasOf; NameOfCoef vr; }
			{ Name Current; Type AssociatedWith; NameOfCoef vr; }
		}
		Constraint {
			{ NameOfCoef vn; EntityType NodesOf; NameOfConstraint CstV; }
			{ NameOfCoef Current; EntityType GroupsOfNodesOf; NameOfConstraint CstI; }
		}
	}
}


Jacobian {
	{ Name J1; Case { 
		{ Region VolAll; Jacobian Vol; }
	}}
}


Integration {
	{ Name I1; Case {{ Type Gauss; Case { 
		{ GeoElement Tetrahedron; NumberOfPoints 4; } 
		{ GeoElement Tetrahedron2; NumberOfPoints 5; } 
	}}}}
	{ Name I2; Case {{ Type Gauss; Case { 
		// For relative L2 error norm 
		{ GeoElement Tetrahedron; NumberOfPoints 5; } 
		{ GeoElement Tetrahedron2; NumberOfPoints 15; } 
	}}}}
}


Formulation {
	{ Name MagFrm; Type FemEquation;
		Quantity {
			{ Name a; Type Local; NameOfSpace Hcurl; }
			{ Name v; Type Local; NameOfSpace Hgrad; }
			{ Name Voltage; Type Global; NameOfSpace Hgrad[Voltage]; }
			{ Name Current; Type Global; NameOfSpace Hgrad[Current]; }
		}
		Equation {
			// Electrokinetic:
			Integral {[ sigma * Dof{d v}, {d v} ];
				In VolWire ; Jacobian J1; Integration I1; }
			Integral {[ Complex[0,1] * omega * sigma * Dof{a}, {d v} ];
				In VolWire ; Jacobian J1; Integration I1; }	
			GlobalTerm { [ -Dof{Current}, {Voltage} ]; In SurWireBack; }
			// Negative sign above is because we want inward current on SurWireBack			
			
			// Magnetic:
			Integral {[ 1/mu0 * Dof{d a}, {d a} ];
				In VolAll; Jacobian J1; Integration I1; }
			Integral {[ Complex[0,1] * omega* sigma * Dof{a}, {a} ];
				In VolWire ; Jacobian J1; Integration I1; }
			Integral {[ sigma * Dof{d v}, {a} ];
				In VolWire ; Jacobian J1; Integration I1; }
		}
	} 
}


Resolution {
	{ Name EmRes;
		System {
			{ Name S1; NameOfFormulation MagFrm; Type Complex; }
		}
		Operation {
			CreateDir["pos"];
			Generate[S1]; Solve[S1]; SaveSolution[S1];
		}
	} 
}


PostProcessing {
	{ Name Post; NameOfFormulation MagFrm;
		Quantity {
			{ Name Bre; Value {Local {
				[ Re[{d a}] ]; 
				In VolAll; Jacobian J1; }}
			}
			{ Name Bim; Value {Local {
				[ Im[{d a}] ]; 
				In VolAll; Jacobian J1; }}
			}
			{ Name P; Value {Local { 
				[ 0.5* sigma*SquNorm[-Complex[0,1]*omega*{a} -{d v}] ];
				// 0.5*sigma*|E|^2, average power hence 0.5 factor
				In VolWire; Jacobian J1; }}
			}
			{ Name J; Value {Local {
				[ Norm[sigma*(-Complex[0,1]*omega*{a} -{d v})] ]; 
				// sigma*|E|
				In VolWire; Jacobian J1; }}
			}
			{ Name Vrms; Value {Local { 
				[ {Voltage} / Sqrt[2] ];
				In SurWireBack; }}
			}
			{ Name Irms; Value {Local { 
				[ coef* {Current} / Sqrt[2] ]; 
				In SurWireBack; }}
			}
			{ Name L2error; Value {Integral {Type Global; 
				[ SquNorm[Vector[0,0,1]*Jex[]-
					sigma*(-Complex[0,1]*omega*{a} -{d v})] 
				]; // square root in PostOperation
				In VolWire; Jacobian J1; Integration I2; }}
			}
			{ Name JexSqu; Value {Integral {Type Global;
				[ SquNorm[Jex[]] ];
				In VolWire; Jacobian J1; Integration I2; }}
			}
		}
	}
}


PostOperation {
	{ Name PostMain; NameOfPostProcessing Post;
		Format Table;
		Operation {			
			Print[ Bre, OnElementsOf VolAll, Format Gmsh, File "pos/wire3d_Bre.pos"];
			Print[ Bim, OnElementsOf VolAll, Format Gmsh, File "pos/wire3d_Bim.pos"];
			Print[ P, OnElementsOf VolWire, Format Gmsh, File "pos/wire3d_P.pos"];
			Print[ J, OnElementsOf VolWire, Format Gmsh, File "pos/wire3d_J.pos"];

			Print[{r, d, meshSize, order, freq}, Format "r=%g, d=%g, ms=%g, o=%g, freq=%g:", File > "output3d.txt"]; 

			Print[ Vrms, OnRegion SurWireBack, StoreInVariable $Vrms ];
			Print[ Irms, OnRegion SurWireBack, StoreInVariable $Irms ];

			Print[ {Re[$Irms], Im[$Irms]}, Format " I = %.8g + i*%.8g [A] (RMS)", File > "output3d.txt"];
			Print[ {Re[$Vrms], Im[$Vrms]}, Format " V = %.8g + i*%.8g [V] (RMS)", File > "output3d.txt" ]; 

			Print[ {$Psim=$Vrms*Conj[$Irms]}, Format " P = %.8g [W]", File > "output3d.txt" ];

			Print[ {$Rsim=Re[$Vrms/$Irms]}, Format " R = %.8g [Ohm]", File > "output3d.txt" ];

			Print[ L2error, OnGlobal, StoreInVariable $L2error ];
			Print[ JexSqu, OnGlobal, StoreInVariable $JexSqu ];
			Print[ {Sqrt[$L2error/$JexSqu]}, Format " RelL2e = %.8g [1]", File > "output3d.txt" ];
		}
	}  
}
