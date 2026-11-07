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
	VolWire    = #{100}; 
	VolVacuum  = #{101};
	SurExt     = #{120};

	VolAll     = #{VolWire, VolVacuum};
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
		{ Region SurExt; Value 0; }
	}}
	{ Name CstI; Case {
		{ Region VolWire; Value currCplx[]; }
	}}
}


FunctionSpace {
	{ Name Hcurl_A; Type Form1P;
		BasisFunction {
			{ Name sn; NameOfCoef an; Function BF_PerpendicularEdge;    
				Support VolAll; Entity NodesOf[All]; }
	}
		Constraint {
			{ NameOfCoef an; EntityType NodesOf; NameOfConstraint CstA; }
		}
	}
	{ Name Hcurl_GradV; Type Form1P; 
		BasisFunction {
			{ Name sr; NameOfCoef dvr; Function BF_RegionZ;
				Support VolWire; Entity VolWire; }
		}
		GlobalQuantity {
			{ Name GradVoltage; Type AliasOf; NameOfCoef dvr; }
			{ Name Current; Type AssociatedWith; NameOfCoef dvr; }
		}
		Constraint {
			{ NameOfCoef Current; EntityType Region; NameOfConstraint CstI; }
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
		{ GeoElement Triangle; NumberOfPoints 3;}
		{ GeoElement Triangle2; NumberOfPoints 4;}
	}}}}
	{ Name I2; Case {{ Type Gauss; Case { 
		// For relative L2 error norm 
		{ GeoElement Triangle; NumberOfPoints 4;}
		{ GeoElement Triangle2; NumberOfPoints 7;}
	}}}}
}


Formulation {
	{ Name MagFrm; Type FemEquation;
		Quantity {
			{ Name a;  Type Local; NameOfSpace Hcurl_A; }
			{ Name dv; Type Local; NameOfSpace Hcurl_GradV; }
			{ Name GradVoltage; Type Global; NameOfSpace Hcurl_GradV[GradVoltage]; }
			{ Name Current; Type Global; NameOfSpace Hcurl_GradV[Current]; }
		}
		Equation {
			// Electrokinetic:
			Integral { [ sigma * Dof{dv}, {dv} ];
				In VolWire; Jacobian J1; Integration I1; }
			Integral { [ Complex[0,1] * omega * sigma * Dof{a}, {dv} ];
				In VolWire; Jacobian J1; Integration I1; }
			GlobalTerm { [ Dof{Current} , {GradVoltage} ]; In VolWire; }

			// Magnetic:
			Integral {[ 1/mu0* Dof{d a}, {d a} ];
				In VolAll; Jacobian J1; Integration I1; }
			Integral { [ sigma * Dof{dv}, {a} ];
				In VolWire; Jacobian J1; Integration I1; }
			Integral {[ Complex[0,1] * omega* sigma * Dof{a}, {a} ];
				In VolWire; Jacobian J1; Integration I1; }
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
				[ 0.5* sigma*SquNorm[-Complex[0,1]*omega*{a} -{dv}] ];
				// 0.5*sigma*|E|^2, average power hence 0.5 factor
				In VolWire; Jacobian J1; }}
			}
			{ Name J; Value {Local {
				[ Norm[sigma*(-Complex[0,1]*omega*{a} -{dv})] ]; 
				// sigma*|E|
				In VolWire; Jacobian J1; }}
			}
			{ Name Vrms; Value {Local { 
				[ -d*{GradVoltage} / Sqrt[2] ];
				In VolWire; }}
			}
			{ Name Irms; Value {Local { 
				[ {Current} / Sqrt[2] ]; 
				In VolWire; }}
			}
			{ Name L2error; Value {Integral {Type Global; 
				[ SquNorm[Vector[0,0,1]*Jex[]-
					sigma*(-Complex[0,1]*omega*{a} -{dv})] 
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
			Print[ Bre, OnElementsOf VolAll, Format Gmsh, File "pos/plate2d_Bre.pos"];
			Print[ Bim, OnElementsOf VolAll, Format Gmsh, File "pos/plate2d_Bim.pos"];
			Print[ P, OnElementsOf VolWire, Format Gmsh, File "pos/plate2d_P.pos"];
			Print[ J, OnElementsOf VolWire, Format Gmsh, File "pos/plate2d_J.pos"];

			Print[{r, d, meshSize, order, freq}, Format "r=%g, d=%g, ms=%g, o=%g, freq=%g:", File > "output2d.txt"]; 

			Print[ Vrms, OnRegion VolWire, StoreInVariable $Vrms ];
			Print[ Irms, OnRegion VolWire, StoreInVariable $Irms ];

			Print[ {Re[$Irms], Im[$Irms]}, Format " I = %.8g + i*%.8g [A] (RMS)", File > "output2d.txt"];
			Print[ {Re[$Vrms], Im[$Vrms]}, Format " V = %.8g + i*%.8g [V] (RMS)", File > "output2d.txt" ]; 

			Print[ {$Psim=$Vrms*Conj[$Irms]}, Format " P = %.8g [W]", File > "output2d.txt" ];

			Print[ {$Rsim=Re[$Vrms/$Irms]}, Format " R = %.8g [Ohm])", File > "output2d.txt" ];

			Print[ L2error, OnGlobal, StoreInVariable $L2error ];
			Print[ JexSqu, OnGlobal, StoreInVariable $JexSqu ];
			Print[ {Sqrt[$L2error/$JexSqu]}, Format " RelL2e = %.8g [1]", File > "output2d.txt" ];
		}
	}  
}
