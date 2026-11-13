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

	SurFrontExSph = #{121};
	SurBackExSph  = #{122};
	SurExt        = #{123};
	SurFrontSph   = #{124};
	SurBackSph    = #{125};

	VolAll    = #{VolSphere,VolCoil,VolVacInt,VolVacExt};
	SurDiriA  = #{SurFrontExSph,SurBackExSph,SurFrontSph,SurBackSph,SurExt};
}


Function {
	Js[] = js * Vector[-Y[],X[],0] / Sqrt[X[]*X[]+Y[]*Y[]];
}


Constraint {
	{ Name CstA; Case {
		{ Region SurDiriA; Value 0; }
	}}
	{ Name CstVBack; Case {
		{ Region SurBackSph; Value 0; }
	}}
	{ Name CstVFront; Case {
		{ Region SurFrontSph; Value 0; }
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
			{ Name sn; NameOfCoef vn; Function BF_Node;    
				Support VolAll; Entity NodesOf[VolSphere, Not SurBackSph]; }
			{ Name sr; NameOfCoef vr; Function BF_GroupOfNodes;
				Support VolAll; Entity GroupsOfNodesOf[SurBackSph]; }
		}
		GlobalQuantity {
			{ Name Voltage; Type AliasOf; NameOfCoef vr; }
			{ Name Current; Type AssociatedWith; NameOfCoef vr; }
		}
		Constraint {
			{ NameOfCoef vn; EntityType NodesOf; NameOfConstraint CstVFront; }
			{ NameOfCoef Voltage; EntityType GroupsOfNodesOf; NameOfConstraint CstVBack; }
		}
	}
}


Jacobian {
	{ Name J1; Case { 
		{ Region VolVacExt; Jacobian VolSphShell{rb, re}; }
		{ Region All; Jacobian Vol; } // All that are unassigned
	}}
}


Integration {
	{ Name I1; Case {{ Type Gauss; Case { 
		{ GeoElement Tetrahedron; NumberOfPoints 4; } 
		{ GeoElement Tetrahedron2; NumberOfPoints 5; } 
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
			Integral {[ sigma * Dof{d v}, {d v} ];
				In VolSphere; Jacobian J1; Integration I1; }
			Integral {[ Complex[0,1] * omega* sigma * Dof{a}, {d v} ];
				In VolSphere; Jacobian J1; Integration I1; }
			GlobalTerm { [ Dof{Current}, {Voltage} ]; In SurBackSph; }

			Integral {[ 1/mu0 * Dof{d a}, {d a} ];
				In VolAll; Jacobian J1; Integration I1; }
			Integral {[ sigma * Dof{d v}, {a} ];
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
				[ 0.5* sigma*SquNorm[-Complex[0,1]*omega*{a} -{d v}] ];
				// 0.5*sigma*|E|^2, average power hence 0.5 factor
				In VolSphere; Jacobian J1; }}
			}
			{ Name J; Value {Local {
				[ Norm[sigma*(-Complex[0,1]*omega*{a} -{d v})] ]; 
				// sigma*|E|
				In VolSphere; Jacobian J1; }}
			}
			{ Name Irms; Value {Local { 
				[ 2* {Current} / Sqrt[2] ];
				// current crossing the half plane, hence mult by 2
				In SurBackSph; }}
			}
			{ Name PAvg; Value {Integral {Type Global; 
				[ 2*coef* 0.5* sigma*SquNorm[-Complex[0,1]*omega*{a} -{d v}] ];
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
			Print[ Bre, OnElementsOf VolAll, Format Gmsh, File "pos/sphere3d_Bre.pos"];
			Print[ Bim, OnElementsOf VolAll, Format Gmsh, File "pos/sphere3d_Bim.pos"];
			Print[ P, OnElementsOf VolSphere, Format Gmsh, File "pos/sphere3d_P.pos"];
			Print[ J, OnElementsOf VolSphere, Format Gmsh, File "pos/sphere3d_J.pos"];

			Print[{meshSize, order, freq}, Format "ms=%g, o=%g, freq=%g:", File > "output3d.txt"]; 

			Print[ Irms, OnRegion SurBackSph, StoreInVariable $Irms ];
			Print[ {Re[$Irms], Im[$Irms]}, Format " I = %.8g + i*%.8g [A] (RMS)", File > "output3d.txt"];

			Print[ PAvg, OnRegion VolSphere, StoreInVariable $PAvg ];
			Print[ {$PAvg}, Format " P = %.8g [W]", File > "output3d.txt" ];
		}
	}  
}
