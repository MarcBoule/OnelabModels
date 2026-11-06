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
	VolPlate  = #{100}; 
	VolVacuum = #{101};
	SurTopBot = #{120};

	VolAll = #{VolPlate, VolVacuum};
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
		{ Region SurTopBot; Value 0; }
	}}
	{ Name CstI; Case {
		{ Region VolPlate; Value currCplx[]; }
	}}
}


FunctionSpace {
	{ Name Hcurl_A; Type Form1P;
		BasisFunction {
			{ Name sn; NameOfCoef an; Function BF_PerpendicularEdge;    
				Support VolAll; Entity NodesOf[All]; }
			If (order > 1)
			{ Name se; NameOfCoef ae; Function BF_PerpendicularEdge_2E; 
				Support VolAll; Entity EdgesOf[All, Not SurTopBot]; }
			EndIf
		}
		Constraint {
			{ NameOfCoef an; EntityType NodesOf; NameOfConstraint CstA; }
		}
	}
	{ Name Hcurl_GradV; Type Form1P; 
		BasisFunction {
			{ Name sr; NameOfCoef dvr; Function BF_RegionZ;
				Support VolPlate; Entity VolPlate; }
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
		{ GeoElement Triangle; NumberOfPoints (order == 1 ? 3 : 4); }
	}}}}
	{ Name I2; Case {{ Type Gauss; Case { 
		// For relative L2 error norm 
		{ GeoElement Triangle; NumberOfPoints (order == 1 ? 4 : 7); }
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
			Integral {[ sigma * Dof{dv}, {dv} ];
				In VolPlate; Jacobian J1; Integration I1; }
			Integral {[ Complex[0,1] * omega* sigma * Dof{a}, {dv} ];
				In VolPlate; Jacobian J1; Integration I1; }
			GlobalTerm { [ Dof{Current}, {GradVoltage} ]; In VolPlate; }

			// Magnetic:
			Integral {[ 1/mu0 * Dof{d a}, {d a} ];
				In VolAll; Jacobian J1; Integration I1;}
			Integral {[ sigma * Dof{dv}, {a} ];
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
				In VolPlate; Jacobian J1; }}
			}
			{ Name J; Value {Local {
				[ Norm[sigma*(-Complex[0,1]*omega*{a} -{dv})] ]; 
				// sigma*|E|
				In VolPlate; Jacobian J1; }}
			}
			{ Name Vrms; Value {Local { 
				[ -d*{GradVoltage} / Sqrt[2] ];
				In VolPlate; }}
			}
			{ Name Irms; Value {Local { 
				[ {Current} / Sqrt[2] ]; 
				In VolPlate; }}
			}
			{ Name L2error; Value {Integral {Type Global; 
				[ SquNorm[Vector[0,0,1]*Jex[]-
					sigma*(-Complex[0,1]*omega*{a} -{dv})] 
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
			Print[ Bre, OnElementsOf VolAll, Format Gmsh, File "plate_Bre.pos"];
			Print[ Bim, OnElementsOf VolAll, Format Gmsh, File "plate_Bim.pos"];
			Print[ P, OnElementsOf VolPlate, Format Gmsh, File "plate_P.pos"];
			Print[ J, OnElementsOf VolPlate, Format Gmsh, File "plate_J.pos"];

			Print[{w, h, d, meshSize, order, freq}, Format "w=%g, h=%g, d=%g, ms=%g, o=%g, freq=%g:", File > "output.txt"];

			Print[ Vrms, OnRegion VolPlate, StoreInVariable $Vrms ];
			Print[ Irms, OnRegion VolPlate, StoreInVariable $Irms ];

			Print[ {Re[$Irms], Im[$Irms]}, Format " I = %.8g + i*%.8g [A] (RMS)", File > "output.txt"];
			Print[ {Re[$Vrms], Im[$Vrms]}, Format " V = %.8g + i*%.8g [V] (RMS)", File > "output.txt" ]; 

			Print[ {$Psim=$Vrms*Conj[$Irms], $Pex=Pex[], ($Psim-$Pex)/$Pex*10^6}, Format " P = %.8g [W] (analyt %.8g, %.3g ppm)", File > "output.txt" ];

			Print[ {$Rsim=Re[$Vrms/$Irms], $Rex=Rex[], ($Rsim-$Rex)/$Rex*10^6}, Format " R = %.8g [Ohm] (analyt %.8g, %.3g ppm)", File > "output.txt" ];

			Print[ L2error, OnGlobal, StoreInVariable $L2error ];
			Print[ JexSqu, OnGlobal, StoreInVariable $JexSqu ];
			Print[ {Sqrt[$L2error/$JexSqu]}, Format " RelL2e = %.8g [1]", File > "output.txt" ]; 
		}
	}  
}
