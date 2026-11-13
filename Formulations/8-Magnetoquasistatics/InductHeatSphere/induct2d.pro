Include "induct_common.pro";

// Parameters
sigma = 36e6;      // Conductivity [S/m)]
freq = 1000;       // Frequency [1/s]
omega = 2*Pi*freq; // Angular frequency [1/s]
js = 5e6;          // Source current density [A/m]

// Constants
mu0  = 1.25663706127e-6;


Group {
 VolSphere = #{100}; 
 VolCoil   = #{101}; 
 VolVacInt = #{102}; 
 VolVacExt = #{103}; 
 
 SurLeft   = #{121};
 SurExt    = #{122};
 
 VolAll    = #{VolSphere,VolCoil,VolVacInt,VolVacExt};
 SurDiriA  = #{SurLeft,SurExt};
}


Function {
	Js[] = js * Vector[0,0,-1];
}


Constraint {
	{ Name CstA; Case {
		{ Region SurDiriA; Value 0; }
	}}
	{ Name CstGV; Case {
		{ Region VolSphere; Value 0; }
	}}
}


FunctionSpace {
	{ Name Hcurl_A; Type Form1P;
		BasisFunction {
			{ Name sn;  NameOfCoef an;  Function BF_PerpendicularEdge;    
				Support VolAll; Entity NodesOf[All]; }
		}
		Constraint {
			{ NameOfCoef an;  EntityType NodesOf; NameOfConstraint CstA; }
		}
	}
	{ Name Hcurl_GradV; Type Form1P; 
		BasisFunction {
			{ Name sr; NameOfCoef dvr; Function BF_RegionZ;
				Support VolSphere; Entity VolSphere; }
		}
		GlobalQuantity {
			{ Name GradVoltage; Type AliasOf; NameOfCoef dvr; }
			{ Name Current; Type AssociatedWith; NameOfCoef dvr; }
		}
		Constraint {
			{ NameOfCoef GradVoltage; EntityType Region; NameOfConstraint CstGV; }
		}
	}
}


Jacobian {
	{ Name J1; Case { 
		{ Region VolVacExt; Jacobian VolAxiSphShell{rb, re}; }
		{ Region All; Jacobian VolAxi; } // All that are unassigned
	}}
}


Integration {
	{ Name I1; Case {{ Type Gauss; Case { 
		{ GeoElement Triangle; NumberOfPoints 4;} 
		{ GeoElement Triangle2; NumberOfPoints 7;} 
	}}}}
}


Formulation {
	{ Name MagFrm; Type FemEquation;
		Quantity {
			{ Name a;  Type Local; NameOfSpace Hcurl_A;}
			{ Name dv; Type Local; NameOfSpace Hcurl_GradV; }
			{ Name GradVoltage; Type Global; NameOfSpace Hcurl_GradV[GradVoltage]; }
			{ Name Current; Type Global; NameOfSpace Hcurl_GradV[Current]; }
		}
		Equation { 
			Integral { [ sigma * Dof{dv}, {dv} ];
				In VolSphere; Jacobian J1; Integration I1; }
			Integral { [ Complex[0,1] * omega * sigma * Dof{a}, {dv} ];
				In VolSphere; Jacobian J1; Integration I1; }
			GlobalTerm { [ Dof{Current}, {GradVoltage} ]; In VolSphere; }
			
			Integral {[ 1/mu0 * Dof{d a}, {d a} ];
				In VolAll; Jacobian J1; Integration I1; }
			Integral { [ sigma * Dof{dv}, {a} ];
				In VolSphere; Jacobian J1; Integration I1; }
			Integral {[ Complex[0,1] * omega* sigma * Dof{a}, {a} ];
				In VolSphere; Jacobian J1; Integration I1; }
			Integral {[ -Js[], {a} ];
				In VolCoil; Jacobian J1; Integration I1; }
		}
	} 
}


Resolution {
	{ Name EmRes;
		System {
			{ Name S1; NameOfFormulation MagFrm; Type Complex;}
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
				In VolSphere; Jacobian J1; }}
			}
			{ Name J; Value {Local {
				[ Norm[sigma*(-Complex[0,1]*omega*{a} -{dv})] ]; 
				// sigma*|E|
				In VolSphere; Jacobian J1; }}
			}
			{ Name Irms; Value {Local { 
				[ 2* -{Current} / Sqrt[2] ];
				// current crossing the half plane, hence mult by 2
				In VolSphere; }}
			}
			{ Name PAvg; Value {Integral {Type Global;
				[ 2*2*Pi* 0.5* sigma*SquNorm[-Complex[0,1]*omega*{a} -{dv}] ];
				// 0.5*sigma*|E|^2, average power hence 0.5 factor
				In VolSphere; Jacobian J1; Integration I1; }}
			}
		}
	}
}


PostOperation {
	{ Name PostMain; NameOfPostProcessing Post;
		Format Table;
		Operation {			
			Print[ Bre, OnElementsOf VolAll, Format Gmsh, File "pos/sphere2d_Bre.pos"];
			Print[ Bim, OnElementsOf VolAll, Format Gmsh, File "pos/sphere2d_Bim.pos"];
			Print[ P, OnElementsOf VolSphere, Format Gmsh, File "pos/sphere2d_P.pos"];
			Print[ J, OnElementsOf VolSphere, Format Gmsh, File "pos/sphere2d_J.pos"];

			Print[{meshSize, order, freq}, Format "ms=%g, o=%g, freq=%g:", File > "output2d.txt"]; 

			Print[ Irms, OnRegion VolSphere, StoreInVariable $Irms ];
			Print[ {Re[$Irms], Im[$Irms]}, Format " I = %.8g + i*%.8g [A] (RMS)", File > "output2d.txt"];

			Print[ PAvg, OnRegion VolSphere, StoreInVariable $PAvg ];
			Print[ {$PAvg}, Format " P = %.8g [W]", File > "output2d.txt" ];
		}
	}  
}


/*
PostProcessing {
 { Name Post; NameOfFormulation MagFrm;
   Quantity {
     { Name V; Value {Local {
      [ {v} ]; 
	  In VolSphere; Jacobian J1;}}
     }
     { Name A; Value {Local {
      [ CompZ[{a}] ]; 
	  In VolAll; Jacobian J1;}}
     }
     { Name Js; Value {Local {
      [ Js[] ]; 
	  In VolCoil; Jacobian J1;}}
     }
     { Name B; Value {Local {
      [ {d a} ]; 
	  In VolAll; Jacobian J1;}}
     }
     { Name Ere; Value {Local {
      [ Re[-Complex[0,1]*omega*{a} -{d v}] ]; 
	  In VolPlate; Jacobian J1;}}
     }
     { Name Eim; Value {Local {
      [ Im[-Complex[0,1]*omega*{a} -{d v}] ]; 
	  In VolPlate; Jacobian J1;}}
     }
     { Name Hre; Value {Local {
      [ Re[{d a}/mu0] ]; 
	  In VolAll; Jacobian J1;}}
     }
     { Name Him; Value {Local {
      [ Im[{d a}/mu0] ]; 
	  In VolAll; Jacobian J1;}}
     }
     { Name J; Value {Local {
      [ Norm[sigma*CompZ[-Complex[0,1]*omega*{a} -{d v}]] ]; 
	  In VolPlate; Jacobian J1;}}
     }
      { Name PcAvg; Value {Integral {Type Global; 
       [ 2*2*Pi * 0.5*sigma* SquNorm[-Complex[0,1]*omega*{a} -{dv}] ];
	   // .5*sigma*|E|^2 
	   // average power hence 0.5 factor
	   // with real scalar sigma, Pc is purely real
       In VolSphere; Jacobian J1; Integration I1; }}
     }
		{ Name Curr; Value { Term {
			[ 2* {Current} ]; In VolSphere;}}
		}


   }
 }
}


PostOperation {
 { Name PostMain; NameOfPostProcessing Post;
   Format Table;
   Operation {
	Print[{order, freq, intPts}, 
	Format "o=%g, freq=%g, intPts=%g:", File > "output.txt"];

	Print[ PcAvg, OnGlobal, StoreInVariable $Pc ];
	Print[ {$Pc}, Format 
	 " Pc = %.8g [W] (FEMM Resistive Loss 6.0841e-007 Watts)", File > "output.txt" ];
	
	Print[ Curr, OnRegion VolSphere, StoreInVariable $Iout ];
	Print[ {Re[$Iout], Im[$Iout]}, Format " Ire = %.8g [A], Iim = %.8g [A]", File > "output.txt"];

	// Print[ B, OnElementsOf VolAll, Format Gmsh, File "plate_B.pos" ];
	//Print[ A, OnElementsOf VolAll, Format Gmsh, File "plate_A.pos", Smoothing ];
	//Print[ Js, OnElementsOf VolCoil, Format Gmsh, File "plate_Js.pos" ];
	

	Print[ VoutRms, OnRegion SurOut, StoreInVariable $Vout ];
	Print[ {Re[$Vout], Im[$Vout]}, Format " Vre = %.8g [V], Vim = %.8g [V]", File > "output.txt" ]; 
	Print[ {$Rsim=Re[$Vout/$Iout], $Rana=R[], ($Rsim-$Rana)/$Rana*10^6}, Format 
	 " R2 = %.8g [Ohm] (analyt %.8g, %.3g ppm)", File > "output.txt" ];
   }
 }  
}
*/