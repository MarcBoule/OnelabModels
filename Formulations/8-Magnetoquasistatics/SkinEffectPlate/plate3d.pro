Include "plate_common.pro";

// Parameters
order = 2;         // Interpolation order (1 and 2 are supported)
sigma = 58e6;      // Conductivity [S/m)]
freq = 1000;       // Frequency [1/s]
omega = 2*Pi*freq; // Angular frequency [1/s]
currMag = 5;       // Current phasor magnitude [A]
currPhase = Pi/3;  // Current phasor phase [rad]

// Constants
mu0  = 1.25663706127e-6; // Reference permeability [H/m]


Group {
	VolPlate      = #{100};
	VolVacuum     = #{101}; 
	SurTopBot     = #{120};
	SurVacFrontBack = #{121};
	SurLeftRight  = #{122};
	SurPlateBack  = #{123};
	SurPlateFront = #{124};

	VolAll        = #{VolPlate,VolVacuum};
	SurPlateFrontBack = #{SurPlateFront,SurPlateBack};
	SurDiriA      = #{SurTopBot,SurVacFrontBack,SurPlateFrontBack};
}


Function {
	currCplx[] = currMag*Complex[Cos[currPhase],Sin[currPhase]];
	
	// Exact solution below is based on section 14.2 in 
	//   "Introductory Engineering Electromagnetics", by
	//   Branko D. Popovic, Addison-Wesley, 1971
	
	// Utility values
	k[] = Sqrt[omega*sigma*mu0/2];
	K[] = Complex[1,1]*k[];
	
	// Exact current density as a function of position Y[]
	Jex[] = currCplx[] * K[]/(2*w) * (Cosh[K[]*Y[]]/Sinh[K[]*h/2]);

	// Exact resistance and exact active power (average real power)
	Rex[] = d*k[]/(2*w*sigma) * (Sinh[k[]*h]+Sin[k[]*h]) / (Cosh[k[]*h]-Cos[k[]*h]);
	Pex[] = 0.5*Rex[]*currMag^2;
}


Constraint {
	{ Name CstA; Case {
		{ Region SurDiriA; Value 0; }
	}}
	{ Name CstV; Case {
		{ Region SurPlateFront; Value 0; }
	}}
	{ Name CstI; Case {
		{ Region SurPlateBack; Value currCplx[]; }
	}}
	{ Name CstGaugeA; Case {
		{ Region VolAll; SubRegion SurDiriA; Value 0; }
	}}
}


FunctionSpace {
	{ Name Hcurl_A; Type Form1;
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
			{ NameOfCoef ae;  EntityType EdgesOf;  NameOfConstraint CstA; }
			If (order > 1)
			{ NameOfCoef a3a; EntityType FacetsOf; NameOfConstraint CstA; }
			{ NameOfCoef a3b; EntityType FacetsOf; NameOfConstraint CstA; }
			EndIf
			{ NameOfCoef ae; EntityType EdgesOfTreeIn; 
				EntitySubType StartingOn; NameOfConstraint CstGaugeA; }
		}
	}
	{ Name Hgrad_V; Type Form0;
		BasisFunction {
			{ Name sn;  NameOfCoef vn;  Function BF_Node;    
				Support VolAll; Entity NodesOf[All, Not SurPlateBack]; }
			If (order > 1)
			{ Name sn2; NameOfCoef vn2; Function BF_Node_2E; 
				Support VolAll; Entity EdgesOf[All, Not SurPlateFrontBack]; }
			EndIf
			{ Name sr; NameOfCoef vr; Function BF_GroupOfNodes;
				Support VolAll; Entity GroupsOfNodesOf[SurPlateBack]; }
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
		{ GeoElement Tetrahedron; NumberOfPoints (order == 1 ? 4 : 5); } 
	}}}}
	{ Name I2; Case {{ Type Gauss; Case { 
		// For relative L2 error norm 
		{ GeoElement Tetrahedron; NumberOfPoints (order == 1 ? 5 : 15); } 
	}}}}
}


Formulation {
	{ Name MagFrm; Type FemEquation;
		Quantity {
			{ Name a; Type Local; NameOfSpace Hcurl_A; }
			{ Name v; Type Local; NameOfSpace Hgrad_V; }
			{ Name Current; Type Global; NameOfSpace Hgrad_V[Current]; }
			{ Name Voltage; Type Global; NameOfSpace Hgrad_V[Voltage]; }
		}
		Equation { 
			// Electrokinetic:
			Integral {[ sigma * Dof{d v}, {d v} ];
				In VolPlate; Jacobian J1; Integration I1; }
			Integral {[ Complex[0,1] * omega * sigma * Dof{a}, {d v} ];
				In VolPlate; Jacobian J1; Integration I1; }
			GlobalTerm { [ -Dof{Current}, {Voltage} ]; In SurPlateBack; } 
			// Negative sign above is because we want inward current on SurPlateBack
			
			// Magnetic:
			Integral {[ 1/mu0 * Dof{d a}, {d a} ];
				In VolAll; Jacobian J1; Integration I1; }
			Integral {[ sigma * Dof{d v}, {a} ];
				In VolPlate; Jacobian J1; Integration I1; }
			Integral {[ Complex[0,1] * omega* sigma * Dof{a}, {a} ];
				In VolPlate; Jacobian J1; Integration I1; }
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
				In VolPlate; Jacobian J1; }}
			}
			{ Name J; Value {Local {
				[ Norm[sigma*(-Complex[0,1]*omega*{a} -{d v})] ]; 
				// sigma*|E|
				In VolPlate; Jacobian J1; }}
			}
			{ Name Vrms; Value {Local { 
				[ {Voltage} / Sqrt[2] ];
				In SurPlateBack; }}
			}
			{ Name Irms; Value {Local { 
				[ {Current} / Sqrt[2] ]; 
				In SurPlateBack; }}
			}
			{ Name L2error; Value {Integral {Type Global; 
				[ SquNorm[Vector[0,0,1]*Jex[]-
					sigma*(-Complex[0,1]*omega*{a} -{d v})] 
				]; // square root in PostOperation
				In VolPlate; Jacobian J1; Integration I2; }}
			}
			{ Name JexSqu; Value {Integral {Type Global;
				[ SquNorm[Jex[]] ];
				In VolPlate; Jacobian J1; Integration I2; }}
			}
		}
	}
}


PostOperation {
	{ Name PostMain; NameOfPostProcessing Post;
		Format Table;
		Operation {
			Print[ Bre, OnElementsOf VolAll, Format Gmsh, File "pos/plate3d_Bre.pos"];
			Print[ Bim, OnElementsOf VolAll, Format Gmsh, File "pos/plate3d_Bim.pos"];
			Print[ P, OnElementsOf VolPlate, Format Gmsh, File "pos/plate3d_P.pos"];
			Print[ J, OnElementsOf VolPlate, Format Gmsh, File "pos/plate3d_J.pos"];

			Print[{w, h, d, meshSize, order, freq}, Format "w=%g, h=%g, d=%g, ms=%g, o=%g, freq=%g:", File > "output3d.txt"];

			Print[ Vrms, OnRegion SurPlateBack, StoreInVariable $Vrms ];
			Print[ Irms, OnRegion SurPlateBack, StoreInVariable $Irms ];

			Print[ {Re[$Irms], Im[$Irms]}, Format " I = %.8g + i*%.8g [A] (RMS)", File > "output3d.txt"];
			Print[ {Re[$Vrms], Im[$Vrms]}, Format " V = %.8g + i*%.8g [V] (RMS)", File > "output3d.txt" ]; 

			Print[ {$Psim=$Vrms*Conj[$Irms], $Pex=Pex[], ($Psim-$Pex)/$Pex*10^6}, Format " P = %.8g [W] (analyt %.8g, %.3g ppm)", File > "output3d.txt" ];

			Print[ {$Rsim=Re[$Vrms/$Irms], $Rex=Rex[], ($Rsim-$Rex)/$Rex*10^6}, Format " R = %.8g [Ohm] (analyt %.8g, %.3g ppm)", File > "output3d.txt" ];

			Print[ L2error, OnGlobal, StoreInVariable $L2error ];
			Print[ JexSqu, OnGlobal, StoreInVariable $JexSqu ];
			Print[ {Sqrt[$L2error/$JexSqu]}, Format " RelL2e = %.8g [1]", File > "output3d.txt" ]; 
		}
	}  
}
