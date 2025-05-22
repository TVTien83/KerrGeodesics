(* ::Package:: *)

(* ::Title:: *)
(*SpecialOrbits subpackage of KerrGeodesics*)


(* ::Chapter:: *)
(*Define usage for public functions*)


(* ::Section:: *)
(*Create Package*)


BeginPackage["KerrGeodesics`SpecialOrbits`",
	{"KerrGeodesics`ConstantsOfMotion`", "KerrGeodesics`OrbitalFrequencies`"}];


(* ::Subsection::Closed:: *)
(*Usage messages*)


KerrGeoPhotonSphereRadius::usage = "KerrGeoPhotonSphereRadius[a,x] returns the radius of the photon sphere."

KerrGeoISCO::usage = "KerrGeoISCO[a,x] returns the location of the innermost stable circular orbit (ISCO) for pro- and retrograde orbits."
KerrGeoISSO::usage = "KerrGeoISSO[a,x] returns the location of the innermost stable spherical orbit (ISSO)."
KerrGeoIBSO::usage = "KerrGeoISBO[a,x] returns the location of the innermost bound spherical orbit (IBSO)."

KerrGeoSeparatrix::usage = "KerrGeoSeparatrix[a,e,x] returns the value of p at the separatrix."

KerrGeoFindResonance::usage = "KerrGeoFindResonance[assoc,{\[Beta]r,\[Beta]\[Theta],\[Beta]\[Phi]}] finds the location of a resonance given {a,x} and one of {p,e} as an association."

KerrGeoOrbitType::usage = "KerrGeoOrbitType[a,p,e,x] outputs whether the parameters correspond to a bound, scatter or plunge orbit."

(*KerrGeoBoundOrbitQ::usage = "KerrGeoBoundOrbitQ[a,p,e,x] tests if the orbital parameters correspond to a bound orbit."
KerrGeoScatterOrbitQ::usage = "KerrGeoScatterOrbitQ[a,p,e,x] tests if the orbital parameters correspond to a scatter orbit."
KerrGeoPlungeOrbitQ::usage = "KerrGeoPlungeOrbitQ[a,p,e,x] tests if the orbital parameters correspond to a plunge orbit."*)


(* ::Subsection:: *)
(*Error messages*)


KerrGeoFindResonance::noresonance = "Resonant orbits only occur for semi-latus rectum in range `1` \[LessEqual] p \[LessEqual] `2`"
KerrGeoFindResonance::noTripleResonance = "No triple resonance; try `1` \[LessEqual] r-integer \[LessEqual] `2`, `3` \[LessEqual] \[Theta]-integer \[LessEqual] `4` or `5` \[LessEqual] \[Phi]-integer \[LessEqual] `6`"
KerrGeoFindResonance::invalida = "Invalid black hole spin parameter. Choose 0 < a \[LessEqual] 1."
KerrGeoFindResonance::invalide = "Invalid orbital eccentricity. Choose 0 \[LessEqual] e \[LessEqual] 1."
KerrGeoFindResonance::invalidx = "Invalid orbital inclination. Choose 0 \[LessEqual] |x| \[LessEqual] 1."
KerrGeoFindResonance::invalidRatio = "Invalid resonant integers. Choose |\[Beta]r| > |\[Beta]\[Theta]|, |\[Beta]r| > |\[Beta]\[Phi]| and \[Beta]\[Theta]*\[Beta]\[Phi] > 0."
KerrGeoFindResonance::assocErr = "Association should have 3 keys including both a and two of {p, e, x}"
KerrGeoFindResonance::assocErrTriple = "Only support association {a, e} for input"
KerrGeoFindResonance::missing = "Resonaces involving the \[Phi] frequency are not yet implemented"
KerrGeoFindResonance::invalidPrograde = "Invalid prograde resonant integers. Choose \[Beta]r/\[Beta]\[Phi] > 0, \[Beta]r/\[Beta]\[Phi] > 0 and \[Beta]\[Theta] < \[Beta]\[Phi]"
KerrGeoFindResonance::invalidRetrograde = "Invalid retrograde resonant integers. Choose \[Beta]r/\[Beta]\[Phi] < 0, \[Beta]r/\[Beta]\[Phi] < 0 and \[Beta]\[Theta] > \[Beta]\[Phi]"  
KerrGeoFindResonance::closeSeparatrix = "No resonance or the exact result is too close to the separatrix, the value is set to the separatrix+10^(-5) instead"


(* ::Subsection::Closed:: *)
(*Being Private section*)


Begin["`Private`"];


(* ::Section::Closed:: *)
(*Innermost stable circular orbit (ISCO)*)


(* ::Text:: *)
(*Schwarzschild ISCO is at r=6M*)


KerrGeoISCO[_?PossibleZeroQ,x_]:=6


(* ::Text:: *)
(*Kerr inner-most circular orbit ISCO from Bardeen, Press, Teukolsky ApJ, 178, p347 (1972), Eq. 2.21*)


KerrGeoISCO[a_,x_/;x^2==1]:=Module[{M=1,Z1,Z2},
	Z1=1+(1-a^2/M^2)^(1/3) ((1+a/M)^(1/3)+(1-a/M)^(1/3));
	Z2=(3a^2/M^2 + Z1^2)^(1/2);
	M(3+Z2-x a ((3-Z1)(3+Z1+2Z2)/(a x)^2)^(1/2))
];


(* ::Section::Closed:: *)
(*Photon Sphere*)


(* ::Text:: *)
(*The photon sphere is at 3M for all radii in Schwarzschild*)


KerrGeoPhotonSphereRadius[_?PossibleZeroQ,x_]:=3


(* ::Text:: *)
(*Radius of photon sphere  for equatorial orbits from Bardeen, Press, Teukolsky ApJ, 178, p347 (1972), Eq. 2.18*)


KerrGeoPhotonSphereRadius[a_,1]:=2(1+Cos[2/3 ArcCos[-a]])
KerrGeoPhotonSphereRadius[a_,-1]:=2(1+Cos[2/3 ArcCos[a]])


(* ::Text:: *)
(*For polar orbits the radius was given by E. Teo, General Relativity and Gravitation, v. 35, Issue 11, p. 1909-1926 (2003), Eq. (14)*)


KerrGeoPhotonSphereRadius[a_,_?PossibleZeroQ]:=1+2Sqrt[1-1/3 a^2]Cos[1/3 ArcCos[(1-a^2)/(1-1/3 a^2)^(3/2)]]


(* ::Text:: *)
(*In the extremal limit we can find the photon sphere radius exactly*)


KerrGeoPhotonSphereRadius[1,x_]:=If[x < Sqrt[3]-1, 1+Sqrt[2] Sqrt[1-x]-x, 1];
KerrGeoPhotonSphereRadius[-1,x_]:=KerrGeoPhotonSphereRadius[1,-x]


(* ::Text:: *)
(*For all other inclinations we have to numerically find the photon sphere radius*)


KerrGeoPhotonSphereRadius[a1_?NumericQ,x0_?NumericQ/;Abs[x0]<=1]/;Precision[{a1,x0}]!=\[Infinity]:=Module[{M=1,a=a1,req,rpolar,\[CapitalPhi],Q,r,u0Sq,prec},
prec=Precision[{a1,x0}];
req=KerrGeoPhotonSphereRadius[a,Sign[x0]];
rpolar=KerrGeoPhotonSphereRadius[a,0];

\[CapitalPhi]=-((r^3-3M r^2+a^2 r+a^2 M)/(a(r-M)));
Q=-((r^3 (r^3-6M r^2+9M^2 r-4a^2 M))/(a^2 (r-M)^2));

u0Sq=((a^2-Q-\[CapitalPhi]^2)+Sqrt[(a^2-Q-\[CapitalPhi]^2)^2+4a^2  Q])/(2a^2);

r/.FindRoot[1-u0Sq-x0^2,Flatten[{r,(req+rpolar)/2,Sort[{req,rpolar}]}],WorkingPrecision->Max[MachinePrecision,prec-1]]//Quiet 
(*The final Quiet[] is there to stop FindRoot complaining about the precision of the argument. 
This seems to be fine near the equatorial plane but might not be ideal for inclincation near the polar orbit*)

]


(* ::Section::Closed:: *)
(*Innermost bound spherical orbits (IBSO)*)


KerrGeoIBSO[_?PossibleZeroQ,x_]:= 4


(* ::Text:: *)
(*Equatorial IBSO results from Bardeen, Press, Teukolsky 1972*)


KerrGeoIBSO[a_,1]:= 2-a+2(1-a)^(1/2)
KerrGeoIBSO[a_,-1]:= 2+a+2(1+a)^(1/2)


(* ::Text:: *)
(*At the IBSO E=1. Solve[KerrGeo[a,p,0,0]==1,p] to get the formula for the IBSO for polar orbits*)


KerrGeoIBSO[a_,0]:=Module[{\[Delta]},
	\[Delta]=27 a^4-8 a^6+3 Sqrt[3] Sqrt[27 a^8-16 a^10];
	1+Sqrt[12-4 a^2-(6 Sqrt[6] (-2+a^2))/Sqrt[6-2 a^2+(4 a^4)/\[Delta]^(1/3)+\[Delta]^(1/3)]-(4 a^4)/\[Delta]^(1/3)-\[Delta]^(1/3)]/Sqrt[6]+Sqrt[6-2 a^2+(4 a^4)/\[Delta]^(1/3)+\[Delta]^(1/3)]/Sqrt[6]
]


KerrGeoIBSO[1,(0|0.)]:=1/3 (3+(54-6 Sqrt[33])^(1/3)+(6 (9+Sqrt[33]))^(1/3))


(* ::Text:: *)
(*The below methods come from L. Stein and N. Warburton arXiv:1912.07609*)


IBSOPoly=(-4+p)^2 p^6+a^8 (-1+x^2)^2+2 a^2 p^5 (-8+2 p+4 x^2-3 p x^2)+2 a^6 p^2 (2-5 x^2+3 x^4)+a^4 p^3 (-8 (1-3 x^2+2 x^4)+p (6-14 x^2+9 x^4));


KerrGeoIBSO[a1_?NumericQ,x1_?NumericQ]/;((Precision[{a1,x1}]!=\[Infinity])&&(1>=x1>=0)):=With[{prec=Precision[{a1,x1}]},
p/.FindRoot[IBSOPoly/.{a->a1,x->x1},{p,KerrGeoIBSO[a1,1],KerrGeoIBSO[a1,0]},WorkingPrecision->Max[MachinePrecision,prec-1]]];

KerrGeoIBSO[a1_?NumericQ,x1_?NumericQ]/;(Precision[{a1,x1}]!=\[Infinity])&&(-1<=x1<0):=With[{prec=Precision[{a1,x1}]},
p/.FindRoot[IBSOPoly/.{a->a1,x->x1},{p,KerrGeoIBSO[a1,0],KerrGeoIBSO[a1,-1]},WorkingPrecision->Max[MachinePrecision,prec-1]]];


(* ::Section::Closed:: *)
(*Separatrix*)


(* ::Text:: *)
(*Negative spin*)


KerrGeoSeparatrix[a_?Negative,e_,x_]:=KerrGeoSeparatrix[-a,e,-x]


(* ::Text:: *)
(*Schwarzschild*)


KerrGeoSeparatrix[_?PossibleZeroQ,e_,x_]:= 6+2e;


(* ::Text:: *)
(*From Glampedakis and Kennefick arXiv:gr-qc/0203086, for a=M we have Subscript[p, s]=1+e*)


KerrGeoSeparatrix[a_/;a==1,e_,x_/;x==1]:= 1+e


(* ::Text:: *)
(*Polar ISSO in extremal case found from playing around with the equations (see L. Stein and N. Warburton arXiv:1912.07609)*)


KerrGeoSeparatrix[a_/;a==1,_?PossibleZeroQ,_?PossibleZeroQ]:=1+Sqrt[3]+Sqrt[3+2 Sqrt[3]]
KerrGeoSeparatrix[a_/;a==1,e_/;e==1,_?PossibleZeroQ]:=2/3 (3+(54-6Sqrt[33])^(1/3)+(6(9+Sqrt[33]))^(1/3))


(* ::Text:: *)
(*For e=1 the Subscript[p, s] is at 2 Subscript[r, ibso]*)


KerrGeoSeparatrix[a_,e_/;e==1,x_]:=2KerrGeoIBSO[a,x]


(* ::Text:: *)
(*The below methods come from L. Stein and N. Warburton arXiv:1912.07609*)


SepPoly=-4 (3+e) p^11+p^12+a^12 (-1+e)^4 (1+e)^8 (-1+x)^4 (1+x)^4-4 a^10 (-3+e) (-1+e)^3 (1+e)^7 p (-1+x^2)^4-4 a^8 (-1+e) (1+e)^5 p^3 (-1+x)^3 (1+x)^3 (7-7 x^2-e^2 (-13+x^2)+e^3 (-5+x^2)+7 e (-1+x^2))+8 a^6 (-1+e) (1+e)^3 p^5 (-1+x^2)^2 (3+e+12 x^2+4 e x^2+e^3 (-5+2 x^2)+e^2 (1+2 x^2))-8 a^4 (1+e)^2 p^7 (-1+x) (1+x) (-3+e+15 x^2-5 e x^2+e^3 (-5+3 x^2)+e^2 (-1+3 x^2))+4 a^2 p^9 (-7-7 e+e^3 (-5+4 x^2)+e^2 (-13+12 x^2))+2 a^8 (-1+e)^2 (1+e)^6 p^2 (-1+x^2)^3 (2 (-3+e)^2 (-1+x^2)+a^2 (e^2 (-3+x^2)-3 (1+x^2)+2 e (1+x^2)))-2 p^10 (-2 (3+e)^2+a^2 (-3+6 x^2+e^2 (-3+2 x^2)+e (-2+4 x^2)))+a^6 (1+e)^4 p^4 (-1+x^2)^2 (-16 (-1+e)^2 (-3-2 e+e^2) (-1+x^2)+a^2 (15+6 x^2+9 x^4+e^2 (26+20 x^2-2 x^4)+e^4 (15-10 x^2+x^4)+4 e^3 (-5-2 x^2+x^4)-4 e (5+2 x^2+3 x^4)))-4 a^4 (1+e)^2 p^6 (-1+x) (1+x) (-2 (11-14 e^2+3 e^4) (-1+x^2)+a^2 (5-5 x^2-9 x^4+4 e^3 x^2 (-2+x^2)+e^4 (5-5 x^2+x^4)+e^2 (6-6 x^2+4 x^4)))+a^2 p^8 (-16 (1+e)^2 (-3+2 e+e^2) (-1+x^2)+a^2 (15-36 x^2+30 x^4+e^4 (15-20 x^2+6 x^4)+4 e^3 (5-12 x^2+6 x^4)+4 e (5-12 x^2+10 x^4)+e^2 (26-72 x^2+44 x^4)));
SepEquat=a^4 (-3-2 e+e^2)^2+p^2 (-6-2 e+p)^2-2 a^2 (1+e) p (14+2 e^2+3 p-e p);
SepPolar=a^6 (-1+e)^2 (1+e)^4+p^5 (-6-2 e+p)+a^2 p^3 (-4 (-1+e) (1+e)^2+(3+e (2+3 e)) p)-a^4 (1+e)^2 p (6+2 e^3+2 e (-1+p)-3 p-3 e^2 (2+p));


pEquatPro[a1_?NumericQ,e1_?NumericQ]/;(Precision[{a1,e1}]!=\[Infinity]):=With[{prec=Precision[{a1,e1}]},
p/.FindRoot[SepEquat/.{a->a1,e->e1},{p,1+e1,6+2e1},WorkingPrecision->Max[MachinePrecision,prec-2]]]
pEquatRet[a1_?NumericQ,e1_?NumericQ]/;(Precision[{a1,e1}]!=\[Infinity]):=With[{prec=Precision[{a1,e1}]},
p/.FindRoot[SepEquat/.{a->a1,e->e1},{p,6+2e1,5+e1+4Sqrt[1+e1]},WorkingPrecision->Max[MachinePrecision,prec-2]]]
pPolar[a1_?NumericQ,e1_?NumericQ]/;(Precision[{a1,e1}]!=\[Infinity]):=With[{prec=Precision[{a1,e1}]},
p/.FindRoot[SepPolar/.{a->a1,e->e1},{p,1+Sqrt[3]+Sqrt[3+2Sqrt[3]],8},WorkingPrecision->Max[MachinePrecision,prec-2]]]


KerrGeoSeparatrix[a1_?NumericQ,e1_?NumericQ,1]/;(Precision[{a1,e1}]!=\[Infinity]):=pEquatPro[a1,e1]
KerrGeoSeparatrix[a1_?NumericQ,e1_?NumericQ,0]/;(Precision[{a1,e1}]!=\[Infinity]):=pPolar[a1,e1]
KerrGeoSeparatrix[a1_?NumericQ,e1_?NumericQ,-1]/;(Precision[{a1,e1}]!=\[Infinity]):=pEquatRet[a1,e1]
KerrGeoSeparatrix[a1_?NumericQ,e1_?NumericQ,x1_?NumericQ]/;((Precision[{a1,e1,x1}]!=\[Infinity])&&(1>x1>0)):=With[{prec=Precision[{a1,e1,x1}]},
p/.FindRoot[SepPoly/.{a->a1,x->x1,e->e1},{p,pEquatPro[a1,e1],pPolar[a1,e1]},WorkingPrecision->Max[MachinePrecision,prec-2]]]
KerrGeoSeparatrix[a1_?NumericQ,e1_?NumericQ,x1_?NumericQ]/;((Precision[{a1,e1,x1}]!=\[Infinity])&&(-1<x1<0)):=With[{prec=Precision[{a1,e1,x1}]},
p/.FindRoot[SepPoly/.{a->a1,x->x1,e->e1},{p,pPolar[a1,e1],12},WorkingPrecision->Max[MachinePrecision,prec-2]]]


(* ::Section::Closed:: *)
(*Innermost stable spherical orbit (ISSO)*)


KerrGeoISSO[a_,x_/;Abs[x]==1]:=KerrGeoISCO[a,x]


KerrGeoISSO[a_,x_]:=KerrGeoSeparatrix[a,0,x]


(* ::Section::Closed:: *)
(*Bound Orbit Q*)


KerrGeoBoundOrbitQ[a_?NumericQ, p_?NumericQ, e_?NumericQ, x_?NumericQ] := Module[{ps},
	If[e > 0, ps = KerrGeoSeparatrix[a,e,x], ps = KerrGeoIBSO[a,x]];
	If[p >= ps && 0 <= e < 1, True, False]
]


(* ::Section::Closed:: *)
(*Scatter Orbit Q*)


(* ::Text:: *)
(*Test whether an orbit is a scatter orbit*)


KerrGeoScatterOrbitQ[a_?NumericQ, p_?NumericQ, e_?NumericQ, x_?NumericQ] := If[p >= KerrGeoSeparatrix[a,e,x] && e >= 1, True, False]


(* ::Section::Closed:: *)
(*Plunge Orbit Q*)


(* ::Text:: *)
(*Test whether an orbit is a plunge orbit. This test is currently not sufficient as we can have unstable orbits below the LSO for circular and spherical orbits. There are also parts of the parameter space which don't correspond to any orbit, i.e., p=0*)


KerrGeoPlungeOrbitQ[a_?NumericQ, p_?NumericQ,e_?NumericQ, x_?NumericQ]:=
	If[KerrGeoBoundOrbitQ[0,p,e,1] == KerrGeoScatterOrbitQ[0,p,e,1] == False, True, False]


(* ::Section::Closed:: *)
(*Orbit type*)


(* ::Text:: *)
(*Output the type of orbit based on the orbital parameters. *)


KerrGeoOrbitType[a_?NumericQ, p_?NumericQ, e_?NumericQ, x_?NumericQ]:=Module[{output,IBSO,ISSO,rph},

	If[PossibleZeroQ[e],
		rph = KerrGeoPhotonSphereRadius[a,x];
		IBSO = KerrGeoIBSO[a,x];
		ISSO = KerrGeoISSO[a,x];
		If[rph < p <= IBSO, output = {"Unbound", "Circular", "Unstable"}];
		If[p == IBSO, output = {"MarginallyBound", "Circular", "Unstable"}];
		If[IBSO < p < ISSO, output = {"Bound", "Circular", "Unstable"}];
		If[p == ISSO, output = {"Bound", "Circular", "MarginallyStable"}];
		If[p > ISSO, output = {"Bound", "Circular", "Stable"}];
		If[!PossibleZeroQ[Abs[x]-1] && p > rph && !PossibleZeroQ[a], AppendTo[output,"Spherical"]];
		
		(*If none of the above. At the moment we say NotClassified as the PlungeOrbitQ is not complete*) 
		If[p <= rph, output = {"NotClassified"}];

		,
		(*If not a circular orbit*)
		If[KerrGeoBoundOrbitQ[a,p,e,x] == True, 
			output = {"Bound","Eccentric"};
			,
			(*If not a bound orbit*)
			If[KerrGeoScatterOrbitQ[a,p,e,x] == True, 
				output = {"Scatter"};
				If[PossibleZeroQ[e-1], AppendTo[output, "Parabolic"]];
				If[e>1, AppendTo[output, "Hyperbolic"]];
				,
				(*If none of the above. At the moment we say NotClassified as the PlungeOrbitQ is not complete*) 
				output = {"NotClassified"}
			];
		];
	];

	If[ output[[1]] != "NotClassified",
		If[PossibleZeroQ[Abs[x]-1], AppendTo[output,"Equatorial"], AppendTo[output,"Inclined"]];
	];

	
	output
		
]


(* ::Section:: *)
(*Resonances*)


(* ::Subsection:: *)
(*Resonance solver*)


Options[ResonanceSolver]={PrecisionGoal->Automatic};


ResonanceSolver[f_, {p_, p0_?NumericQ, pSep_?NumericQ}, opts:OptionsPattern[]]:=
Module[{pp},
If[p0<pSep,
	Quiet[Check[Return[p/.FindRoot[f[p], {p, pSep+10^(-5), pSep, Infinity}]],{}]];
	Message[KerrGeoFindResonance::closeSeparatrix];
	Return[pSep+10^(-5)];
];
Quiet[Check[Return[p/.FindRoot[f[p], {p, p0, pSep, Infinity}]],{}]];
Quiet[Check[Return[p/.FindRoot[f[p], {p, pSep+10^(-5), pSep, Infinity}]],{}]];
Message[KerrGeoFindResonance::closeSeparatrix];
Return[pSep+10^(-5)];
];


(* ::Subsection:: *)
(*r\[Theta]-resonances*)


(* ::Text:: *)
(*Reference: Brink, Geyer, and Hinderer, Phys. Rev. D 91 (2015), arXiv:1501.07728*)
(*	- Root finding methods are based on those described in Sec.  VE*)


(* ::Subsubsection:: *)
(*Testing functions*)


ValidAQ[a_]:=Module[{aFlag=True},
	If[a^2>1||a^2==0, Message[KerrGeoFindResonance::invalida]; aFlag=False];
	aFlag
];


ValidEQ[e_]:=Module[{eFlag=True},
	If[e<0||e>1, Message[KerrGeoFindResonance::invalide]; eFlag=False];
	eFlag
];


ValidXQ[x_]:=Module[{xFlag=True},
	If[x^2>1, Message[KerrGeoFindResonance::invalidx]; xFlag=False];
	xFlag
];


ValidResIntQ[\[Beta]r_/;\[Beta]r!=0,\[Beta]\[Theta]_/;\[Beta]\[Theta]!=0,\[Beta]\[Phi]_/;\[Beta]\[Phi]!=0]:=Module[{resFlag},
	If[(Abs[\[Beta]r]>Abs[\[Beta]\[Theta]]), Message[KerrGeoFindResonance::invalidRatio]; resFlag=False];
	If[(Abs[\[Beta]r]>Abs[\[Beta]\[Phi]]), Message[KerrGeoFindResonance::invalidRatio]; resFlag=False];
	If[(\[Beta]\[Theta] \[Beta]\[Phi]<0), Message[KerrGeoFindResonance::invalidRatio]; resFlag=False];
	resFlag
];
ValidResIntQ[\[Beta]r_/;\[Beta]r!=0,\[Beta]\[Theta]_/;\[Beta]\[Theta]!=0,0]:=Module[{resFlag},
	If[(Abs[\[Beta]r]>Abs[\[Beta]\[Theta]]), Message[KerrGeoFindResonance::invalidRatio]; resFlag=False];
	resFlag
];
ValidResIntQ[\[Beta]r_/;\[Beta]r!=0,0,\[Beta]\[Phi]_/;\[Beta]\[Phi]!=0]:=Module[{resFlag},
	If[(Abs[\[Beta]r]>Abs[\[Beta]\[Phi]]), Message[KerrGeoFindResonance::invalidRatio]; resFlag=False];
	resFlag
];
ValidResIntQ[0,\[Beta]\[Theta]_/;\[Beta]\[Theta]!=0,\[Beta]\[Phi]_/;\[Beta]\[Phi]!=0]:=Module[{resFlag},
	If[(\[Beta]\[Theta] \[Beta]\[Phi]<0), Message[KerrGeoFindResonance::invalidRatio]; resFlag=False];
	resFlag
];


(* ::Subsubsection:: *)
(*Useful functions for root-finding*)


(* ::Text:: *)
(*See Eq 25 for definitions of these quantities*)


del1y1[a_/;a!=0,p_,e_,x_/;x^2<1]:=-e p Sqrt[\[Pi]p[a,p,e,x]^2-4 \[Pi]x[a,p,e,x]]/((1-e^2)\[Pi]x[a,p,e,x]+p^2-p \[Pi]p[a,p,e,x]);
del2y2[a_/;a!=0,p_,e_,x_/;x^2<1]:=(1-e^2)a^4 (1-x^2)^2/((1-e^2)a^4 (1-x^2)^2-2p^2 \[Pi]x[a,p,e,x]);
y1y2[a_/;a!=0,p_,e_,x_/;x^2<1]:=2a^2 (1-x^2)((1-e^2)\[Pi]x[a,p,e,x]+p^2-p \[Pi]p[a,p,e,x])/(a^4(e^2-1)(1-x^2)^2+2 p^2 \[Pi]x[a,p,e,x]);

del1y1[a_,p_,e_,x_/;x^2==1]:=-e \[Pi]p[a,p,e,x]/(p-\[Pi]p[a,p,e,x]);
del2y2[a_,p_,e_,x_/;x^2==1]:=0;
y1y2[a_,p_,e_,x_/;x^2==1]:=(p-\[Pi]p[a,p,e,x])/(p+2\[Pi]p[a,p,e,x]);

del1y1[0,p_,e_,x_/;x^2<1]:=del1y1[0,p,e,1];
del2y2[0,p_,e_,x_/;x^2<1]:=del2y2[0,p,e,1];
y1y2[0,p_,e_,x_/;x^2<1]:=y1y2[0,p,e,1];


(* ::Text:: *)
(*\[Pi]p = r3 + r4, \[Pi]x = r3*r4, where (dr/d\[Lambda])^2 = -(\[Mu]^2+ \[ScriptCapitalE]^2)(r-r1)(r-r2)(r-r3)(r-r4) with r1 >= r2 >= r3 > r4. *)


\[Pi]p[a_/;a!=0,p_,e_,x_/;x^2<1]:=((-p^2 + a^2*(-1 + e^2)*(-1 + x^2))*(a^2*(-1 + x^2) + \[Pi]x[a,p,e,x]))/(2*a^2*p*(-1 + x^2));

\[Pi]p[0,p_,e_,x_/;x^2<1]:=2p/(p-4);

\[Pi]p[a_,p_,e_/;e<1,1]:=(2*(a^4*(-1+e^2)*p+(-4+p)*p^3+a^2*p^2*(3+e^2+p)-2*Sqrt[a^2*p^3*(a^4*(-1+e^2)^2+(-4*e^2+(-2+p)^2)*p^2+2*a^2*p*(-2+p+e^2*(2+p)))]))/(a^4*(-1+e^2)^2+(-4+p)^2*p^2+2*a^2*(-1+e^2)*p*(4+p));
\[Pi]p[a_,p_,e_/;e==1,1]:=(-4*Sqrt[p^3*(a^3+a*(-2+p)*p)^2]+2*p*(-a^4+(-4+p)*p^2+a^2*p*(3+p)))/(((a-p)^2-4*p)*(-4*p+(a+p)^2));

\[Pi]p[a_,p_,e_/;e<1,-1]:=(2*(a^4*(-1+e^2)*p+(-4+p)*p^3+a^2*p^2*(3+e^2+p)+2*Sqrt[a^2*p^3*(a^4*(-1+e^2)^2+(-4*e^2+(-2+p)^2)*p^2+2*a^2*p*(-2+p+e^2*(2+p)))]))/(a^4*(-1+e^2)^2+(-4+p)^2*p^2+2*a^2*(-1+e^2)*p*(4+p));
\[Pi]p[a_,p_,e_/;e==1,-1]:=(2*(-(a^4*p)+(-4+p)*p^3+a^2*p^2*(3+p)+2*Sqrt[p^3*(a^3+a*(-2+p)*p)^2]))/(a^4+(-4+p)^2*p^2-2*a^2*p*(4+p));


\[Pi]x[a_,p_,e_,x_/;x>=0]:=(a^2*(1 - x^2)*((-4 + p)*p^7 + a^8*(-1 + e^2)^4*(-1 + x^2)^2 + 2*a^6*(-1 + e^2)^2*p*(1 - x^2)*((1 + e^2)*p + (-2 + p + e^2*(2 + p))*(1 - x^2)) + 
   a^4*p^3*((-3 + 2*e^2 + e^4)*p + 4*(-4 + 3*p + e^4*(4 + p))*(1 - x^2) + (-1 + e^2)*(-4 + e^2*(-12 + p) + 3*p)*(-1 + x^2)^2) + 
   2*a^2*p^2*((-1 - e^2)*p^4*(-2 + x^2) + 8*(-1 + e^2)*p^2*(-1 + x^2) + 2*p^3*(-3 - e^2 + 4*(1 + e^2)*x^2) - 
     4*Sqrt[-1/(a^2*p)(((a^4*(-1 + e^2)^2 + (-4*e^2 + (-2 + p)^2)*p^2 + 2*a^2*p*(-2 + p + e^2*(2 + p)))*x^2*(-p^2 + a^2*(-1 + e)^2*(-1 + x^2))*(-p^2 + a^2*(1 + e)^2*(-1 + x^2))*
          (-p^2 + a^2*(-1 + e^2)*(-1 + x^2))))])))/((-4 + p)^2*p^6 + a^8*(-1 + e^2)^4*(-1 + x^2)^2 + 2*a^2*p^5*((-1 + e^2)*(4 + p) + (-4 + e^2*(-12 + p) + 3*p)*(1 - x^2)) + 
  2*a^6*(-1 + e^2)^2*p^2*(-1 + x^2)*(-2 + 3*x^2 + e^2*(-2 + x^2)) + a^4*p^3*(-8*(-1 + e^2)^2*(1 - 3*x^2 + 2*x^4) + p*((-1 + e^2)^2 - 4*(-1 + e^4)*(-1 + x^2) + (3 + e^2)^2*(-1 + x^2)^2)));

\[Pi]x[a_,p_,e_,x_/;x<0]:=(a^2*(1 - x^2)*((-4 + p)*p^7 + a^8*(-1 + e^2)^4*(-1 + x^2)^2 + 2*a^6*(-1 + e^2)^2*p*(1 - x^2)*((1 + e^2)*p + (-2 + p + e^2*(2 + p))*(1 - x^2)) + 
   a^4*p^3*((-3 + 2*e^2 + e^4)*p + 4*(-4 + 3*p + e^4*(4 + p))*(1 - x^2) + (-1 + e^2)*(-4 + e^2*(-12 + p) + 3*p)*(-1 + x^2)^2) + 
   2*a^2*p^2*(-((1 + e^2)*p^4*(-2 + x^2)) + 8*(-1 + e^2)*p^2*(-1 + x^2) + 2*p^3*(-3 - e^2 + 4*(1 + e^2)*x^2) + 
     4*Sqrt[-1/(a^2*p)(((a^4*(-1 + e^2)^2 + (-4*e^2 + (-2 + p)^2)*p^2 + 2*a^2*p*(-2 + p + e^2*(2 + p)))*x^2*(-p^2 + a^2*(-1 + e)^2*(-1 + x^2))*(-p^2 + a^2*(1 + e)^2*(-1 + x^2))*
          (-p^2 + a^2*(-1 + e^2)*(-1 + x^2))))])))/((-4 + p)^2*p^6 + a^8*(-1 + e^2)^4*(-1 + x^2)^2 + 2*a^2*p^5*((-1 + e^2)*(4 + p) + (-4 + e^2*(-12 + p) + 3*p)*(1 - x^2)) + 
  2*a^6*(-1 + e^2)^2*p^2*(-1 + x^2)*(-2 + 3*x^2 + e^2*(-2 + x^2)) + a^4*p^3*(-8*(-1 + e^2)^2*(1 - 3*x^2 + 2*x^4) + p*((-1 + e^2)^2 - 4*(-1 + e^4)*(-1 + x^2) + (3 + e^2)^2*(-1 + x^2)^2)));


(* ::Text:: *)
(*Carlson's integral of the first kind, which is defined in Eq. (B12) of reference listed above*)


Rf[0,alpha_,beta_]:=beta^(-1/2)EllipticK[1-alpha/beta];


(* ::Text:: *)
(*r\[Theta] ratio*)


r\[Theta]Ratio[a_, p_, e_, x_]:=Sqrt[y1y2[a,p,e,x]]Rf[0,1+del2y2[a,p,e,x],1-del2y2[a,p,e,x]]/Rf[0,1+del1y1[a,p,e,x],1-del1y1[a,p,e,x]];


(* ::Subsubsection::Closed:: *)
(*Given (a,e,x) and Subscript[\[CapitalOmega], r]/Subscript[\[CapitalOmega], \[Theta]]= Subscript[\[Beta], r]/Subscript[\[Beta], \[Theta]] find p*)


(* ::Text:: *)
(*Reference: Brink, Geyer, and Hinderer, Phys. Rev. D 91 (2015), arXiv:1501.07728*)
(*	Root finding methods are based on those described in Sec.  VE*)


Options[KerrGeoOrbitRThetaResonantP]={PrecisionGoal->Automatic};


KerrGeoOrbitRThetaResonantP[a_?NumericQ, e_?NumericQ, x_?NumericQ, {\[Beta]r_Integer, \[Beta]\[Theta]_Integer}, opts:OptionsPattern[]]:=
Module[{pg,ratio,argpg,resonantEqn,pStar,pp,pgTest,pSep},
	(*See if the user specified some precision goal *)
	If[Not@ValidAQ[a]||Not@ValidEQ[e]||Not@ValidXQ[x]||Not@ValidResIntQ[\[Beta]r,\[Beta]\[Theta],0],Abort[]];
	pg=OptionValue[PrecisionGoal];
	ratio=\[Beta]r/\[Beta]\[Theta];
	If[x>0&&ratio<0, Message[KerrGeoFindResonance::invalidPrograde]; Abort[]];
	If[x<0&&ratio>0, Message[KerrGeoFindResonance::invalidRetrograde]; Abort[]];
	(* pStar provides an initial guess for p. Note that p = pStar for e = 0 and a = 0*)
	pStar=6/(1-ratio^2); (* See Eq 30 *)
	pSep=KerrGeoSeparatrix[a,e,x];
	
	(* Resonant condition defined by the equation below *)
	resonantEqn[p_?NumericQ]:=r\[Theta]Ratio[a, p, e, x]-Abs[ratio];
	
	(* Working precision of the root-finding method is based on the precision specified
	 by the PrecisionGoal option, or the precision of the arguments. MachinePrecision
	 is the default precision if not other precision specifications are made. *)
	argpg=Precision[{resonantEqn[pStar],a,e,x,ratio}];
	If[argpg==$MachinePrecision,pg=argpg];
	If[(Not@NumericQ[pg]||pg>argpg),pg=argpg];
	If[pg==Infinity,pg=$MachinePrecision];
	
	If[pg==$MachinePrecision,
		ResonanceSolver[resonantEqn,{pp, pStar, pSep}],
		ResonanceSolver[resonantEqn,{pp, pStar, pSep, PrecisionGoal->pg, WorkingPrecision->pg}],
		ResonanceSolver[resonantEqn,{pp, pStar, pSep, PrecisionGoal->pg, WorkingPrecision->pg}]
	]
];


(* ::Subsubsection::Closed:: *)
(*Given (a,p,x) and Subscript[\[CapitalOmega], r]/Subscript[\[CapitalOmega], \[Theta]]= Subscript[\[Beta], r]/Subscript[\[Beta], \[Theta]] find e*)


Options[KerrGeoOrbitRThetaResonantE]={PrecisionGoal->Automatic};


KerrGeoOrbitRThetaResonantE[a_?NumericQ, p_?NumericQ, x_?NumericQ, {\[Beta]r_Integer, \[Beta]\[Theta]_Integer}, opts:OptionsPattern[]]:=
Module[{pg,argpg,resonantEqn,e0Test,e1Test,eGuess,ee,ratio},
	(*See if the user specified some precision goal *)
	If[Not@ValidAQ[a]||Not@ValidXQ[x]||Not@ValidResIntQ[\[Beta]r,\[Beta]\[Theta],0],Abort[]];
	pg=OptionValue[PrecisionGoal];
	ratio=\[Beta]r/\[Beta]\[Theta];
	If[x>0&&ratio<0, Message[KerrGeoFindResonance::invalidPrograde]; Abort[]];
	If[x<0&&ratio>0, Message[KerrGeoFindResonance::invalidRetrograde]; Abort[]];
	
	(* Test to see if there is an eccentricity that will lead to a bound orbital resonance 
	 based on the provided values of a, p, x *)
	e0Test=KerrGeoOrbitRThetaResonantP[a,0,x,{\[Beta]r,\[Beta]\[Theta]}, opts];
	e1Test=KerrGeoOrbitRThetaResonantP[a,1,x,{\[Beta]r,\[Beta]\[Theta]}, opts];
	If[p<e0Test||e1Test<p, Message[KerrGeoFindResonance::noresonance, e0Test, e1Test]; Abort[];];
	If[p==e0Test,Return[0]];
	If[p==e1Test,Return[1]];
	eGuess=(p-e0Test)/(e1Test-e0Test);

	(* Resonant condition defined by the equation below *)
	resonantEqn[e_?NumericQ]:=r\[Theta]Ratio[a, p, e, x]-Abs[ratio];
	
	(* Working precision of the root-finding method is based on the precision specified
	 by the PrecisionGoal option, or the precision of the arguments. MachinePrecision
	 is the default precision if not other precision specifications are made. *)
	argpg=Precision[{resonantEqn[eGuess],a,p,x,ratio}];
	If[argpg==$MachinePrecision,pg=argpg];
	If[(Not@NumericQ[pg]||pg>argpg),pg=argpg];
	If[pg==Infinity,pg=$MachinePrecision];

	If[pg==$MachinePrecision,
		Re[ee/.FindRoot[resonantEqn[ee],{ee,eGuess}]],
		Re[ee/.FindRoot[resonantEqn[ee],{ee,eGuess},PrecisionGoal->pg,WorkingPrecision->pg]],
		Re[ee/.FindRoot[resonantEqn[ee],{ee,eGuess},PrecisionGoal->pg,WorkingPrecision->pg]]
	]
];


(* ::Subsubsection::Closed:: *)
(*Given (a,p,e) and Subscript[\[CapitalOmega], r]/Subscript[\[CapitalOmega], \[Theta]]= Subscript[\[Beta], r]/Subscript[\[Beta], \[Theta]] find x*)


Options[KerrGeoOrbitRThetaResonantX]={PrecisionGoal->Automatic};


KerrGeoOrbitRThetaResonantX[a_?NumericQ, p_?NumericQ, e_?NumericQ, {\[Beta]r_Integer, \[Beta]\[Theta]_Integer}, opts:OptionsPattern[]]:=
Module[{pg,argpg,resonantEqn,x0Test,x1Test,xGuess,xx,xxx,ratio},
	(*See if the user specified some precision goal *)
	If[Not@ValidAQ[a]||Not@ValidEQ[e]||Not@ValidResIntQ[\[Beta]r,\[Beta]\[Theta],0],Abort[]];

	pg=OptionValue[PrecisionGoal];
	ratio=\[Beta]r/\[Beta]\[Theta];
	(* Test to see if there is an inclination angle that will lead to a bound orbital resonance 
	 based on the provided values of a, p, e *)
	x0Test=KerrGeoOrbitRThetaResonantP[a,e,0,{\[Beta]r,\[Beta]\[Theta]},opts];
	x1Test=KerrGeoOrbitRThetaResonantP[a,e,Sign[ratio],{\[Beta]r,\[Beta]\[Theta]},opts];
	
	If[ratio>0 && (p<x1Test||x0Test<p), Message[KerrGeoFindResonance::noresonance, x1Test, x0Test]; Abort[];];
	If[ratio<0 && (p<x0Test||x1Test<p), Message[KerrGeoFindResonance::noresonance, x0Test, x1Test]; Abort[];];
	If[p==x0Test,Return[0]];
	If[p==x1Test,Return[Sign[ratio]]];
	xGuess=Sign[ratio](p-x0Test)/(x1Test-x0Test);

	(* Resonant condition defined by the equation below *)
	resonantEqn[x_?NumericQ]:=r\[Theta]Ratio[a, p, e, x]-Abs[ratio];
	
	(* Working precision of the root-finding method is based on the precision specified
	 by the PrecisionGoal option, or the precision of the arguments. MachinePrecision
	 is the default precision if not other precision specifications are made. *)
	argpg=Precision[{a,p,e,ratio}];
	If[argpg==$MachinePrecision,pg=0.9argpg];
	If[(Not@NumericQ[pg]||pg>argpg),pg=0.9argpg];
	If[pg==Infinity,pg=$MachinePrecision];

	If[pg==$MachinePrecision,
		xxx=Re[xx/.FindRoot[resonantEqn[xx],{xx,xGuess}]],
		xxx=Re[xx/.FindRoot[resonantEqn[xx],{xx,xGuess},PrecisionGoal->pg,WorkingPrecision->0.95argpg]],
		xxx=Re[xx/.FindRoot[resonantEqn[xx],{xx,xGuess},PrecisionGoal->pg,WorkingPrecision->0.95argpg]]
	]
];


(* ::Subsection:: *)
(*r\[Phi]-resonances*)


(* ::Subsubsection::Closed:: *)
(*Testing functions*)


(* ::Text:: *)
(*Same as those of r\[Theta]-resonance*)


(* ::Subsubsection::Closed:: *)
(*Useful functions for root-finding*)


En2[a_, p_, e_/;e<1, x_]:=1-(2(1-e^2))/(2p+(1-e^2)\[Pi]p[a,p,e,x]);
En2[a_,p_,1, x_]:=1;

L2[a_, p_, e_, x_]:=(2 (a^2 (-1+e^2)+p^2) (a^2-\[Pi]x[a,p,e,x])+4 a^2 p \[Pi]p[a,p,e,x])/(a^2 (2 p+(1-e^2)\[Pi]p[a,p,e,x]));
L2[a_, p_, e_, 0]:=0;

Q[a_, p_, e_, x_/;x^2<1]:=(2 p^2 \[Pi]x[a,p,e,x])/(a^2 (2 p+(1-e^2)\[Pi]p[a,p,e,x]));
Q[a_, p_, e_, x_/;x^2==1]:=0;

r3[a_, p_, e_, x_/;x^2<1]:=(\[Pi]p[a,p,e,x]+Sqrt[\[Pi]p[a,p,e,x]^2-4 \[Pi]x[a,p,e,x]])/2;
r3[a_, p_, e_, x_/;x^2==1]:=\[Pi]p[a,p,e,x];

r4[a_, p_, e_, x_/;x^2<1]:=(\[Pi]p[a,p,e,x]-Sqrt[\[Pi]p[a,p,e,x]^2-4 \[Pi]x[a,p,e,x]])/2;
r4[a_, p_, e_, x_/;x^2==1]:=0;


\[Pi]xReduced[a_,p_,e_,x_/;x>=0]:=(a^2*((-4 + p)*p^7 + a^8*(-1 + e^2)^4*(-1 + x^2)^2 + 2*a^6*(-1 + e^2)^2*p*(1 - x^2)*((1 + e^2)*p + (-2 + p + e^2*(2 + p))*(1 - x^2)) + 
   a^4*p^3*((-3 + 2*e^2 + e^4)*p + 4*(-4 + 3*p + e^4*(4 + p))*(1 - x^2) + (-1 + e^2)*(-4 + e^2*(-12 + p) + 3*p)*(-1 + x^2)^2) + 
   2*a^2*p^2*((-1 - e^2)*p^4*(-2 + x^2) + 8*(-1 + e^2)*p^2*(-1 + x^2) + 2*p^3*(-3 - e^2 + 4*(1 + e^2)*x^2) - 
     4*Sqrt[-1/(a^2*p)(((a^4*(-1 + e^2)^2 + (-4*e^2 + (-2 + p)^2)*p^2 + 2*a^2*p*(-2 + p + e^2*(2 + p)))*x^2*(-p^2 + a^2*(-1 + e)^2*(-1 + x^2))*(-p^2 + a^2*(1 + e)^2*(-1 + x^2))*
          (-p^2 + a^2*(-1 + e^2)*(-1 + x^2))))])))/((-4 + p)^2*p^6 + a^8*(-1 + e^2)^4*(-1 + x^2)^2 + 2*a^2*p^5*((-1 + e^2)*(4 + p) + (-4 + e^2*(-12 + p) + 3*p)*(1 - x^2)) + 
  2*a^6*(-1 + e^2)^2*p^2*(-1 + x^2)*(-2 + 3*x^2 + e^2*(-2 + x^2)) + a^4*p^3*(-8*(-1 + e^2)^2*(1 - 3*x^2 + 2*x^4) + p*((-1 + e^2)^2 - 4*(-1 + e^4)*(-1 + x^2) + (3 + e^2)^2*(-1 + x^2)^2)));
\[Pi]xReduced[a_,p_,e_,x_/;x<0]:=(a^2*((-4 + p)*p^7 + a^8*(-1 + e^2)^4*(-1 + x^2)^2 + 2*a^6*(-1 + e^2)^2*p*(1 - x^2)*((1 + e^2)*p + (-2 + p + e^2*(2 + p))*(1 - x^2)) + 
   a^4*p^3*((-3 + 2*e^2 + e^4)*p + 4*(-4 + 3*p + e^4*(4 + p))*(1 - x^2) + (-1 + e^2)*(-4 + e^2*(-12 + p) + 3*p)*(-1 + x^2)^2) + 
   2*a^2*p^2*(-((1 + e^2)*p^4*(-2 + x^2)) + 8*(-1 + e^2)*p^2*(-1 + x^2) + 2*p^3*(-3 - e^2 + 4*(1 + e^2)*x^2) + 
     4*Sqrt[-1/(a^2*p)(((a^4*(-1 + e^2)^2 + (-4*e^2 + (-2 + p)^2)*p^2 + 2*a^2*p*(-2 + p + e^2*(2 + p)))*x^2*(-p^2 + a^2*(-1 + e)^2*(-1 + x^2))*(-p^2 + a^2*(1 + e)^2*(-1 + x^2))*
          (-p^2 + a^2*(-1 + e^2)*(-1 + x^2))))])))/((-4 + p)^2*p^6 + a^8*(-1 + e^2)^4*(-1 + x^2)^2 + 2*a^2*p^5*((-1 + e^2)*(4 + p) + (-4 + e^2*(-12 + p) + 3*p)*(1 - x^2)) + 
  2*a^6*(-1 + e^2)^2*p^2*(-1 + x^2)*(-2 + 3*x^2 + e^2*(-2 + x^2)) + a^4*p^3*(-8*(-1 + e^2)^2*(1 - 3*x^2 + 2*x^4) + p*((-1 + e^2)^2 - 4*(-1 + e^4)*(-1 + x^2) + (3 + e^2)^2*(-1 + x^2)^2)));


kr[a_,p_,e_,x_]:=(p-p (1-e)/(1+e))/(p-r3[a,p,e,x](1-e))(r3[a,p,e,x]-r4[a,p,e,x])/(p/(1+e)-r4[a,p,e,x])
k\[Theta][a_,p_,e_,x_]:=(1-x^2)(1-e^2)a^4/(p^2 \[Pi]xReduced[a,p,e,x]);
rPlus[a_]:=1+Sqrt[1-a^2];
rMinus[a_]:=1-Sqrt[1-a^2];
hPlus[a_,p_,e_,x_]:=(p-p (1-e)/(1+e))/(p-r3[a,p,e,x](1-e))(r3[a,p,e,x]-rPlus[a])/(p/(1+e)-rPlus[a]);
hMinus[a_,p_,e_,x_]:=(p-p (1-e)/(1+e))/(p-r3[a,p,e,x](1-e))(r3[a,p,e,x]-rMinus[a])/(p/(1+e)-rMinus[a]);
prefactor\[Theta][a_,p_,e_,x_]:=2/Pi Sqrt[((a^2(-1+e^2)+p^2)(a^2-\[Pi]x[a,p,e,x])+2 a^2 p \[Pi]p[a,p,e,x])/(p^2 \[Pi]xReduced[a,p,e,x])];
prefactorR[a_,p_,e_,x_]:=2a/(Pi(rPlus[a]-rMinus[a])) Sqrt[(2p+(1-e^2)\[Pi]p[a,p,e,x])/(2(p-r3[a,p,e,x](1-e))(p-r4[a,p,e,x](1+e)))];
termRPlus[a_,p_,e_,x_]:=(2Sqrt[En2[a,p,e,x]] rPlus[a]-a Sign[x] Sqrt[L2[a,p,e,x]])/(r3[a,p,e,x]-rPlus[a])(EllipticK[kr[a,p,e,x]]-(p/(1+e)-r3[a,p,e,x])/(p/(1+e)-rPlus[a])EllipticPi[hPlus[a,p,e,x], kr[a,p,e,x]]);
termRMinus[a_,p_,e_,x_]:=(2Sqrt[En2[a,p,e,x]] rMinus[a]-a Sign[x] Sqrt[L2[a,p,e,x]])/(r3[a,p,e,x]-rMinus[a])(EllipticK[kr[a,p,e,x]]-(p/(1+e)-r3[a,p,e,x])/(p/(1+e)-rMinus[a])EllipticPi[hMinus[a,p,e,x], kr[a,p,e,x]]);


term\[Theta][a_,p_,e_,x_]:=prefactor\[Theta][a,p,e,x] EllipticPi[1-x^2, k\[Theta][a,p,e,x]];
termR[a_,p_,e_,x_]:=prefactorR[a,p,e,x](termRPlus[a,p,e,x]-termRMinus[a,p,e,x]);


r\[Phi]Ratio[a_, p_, e_, x_/;x!=0]:=Sign[x]/(Sign[x] term\[Theta][a,p,e,x]/r\[Theta]Ratio[a,p,e,x]+termR[a,p,e,x]);
r\[Phi]Ratio[a_, p_, e_, x_/;x==0]:={1/(1/r\[Theta]Ratio[a,p,e,x]+termR[a,p,e,x]),-1/(-1/r\[Theta]Ratio[a,p,e,x]+termR[a,p,e,x])};


(* ::Subsubsection::Closed:: *)
(*Given (a,e,x) and Subscript[\[CapitalOmega], r]/Subscript[\[CapitalOmega], \[Phi]]= Subscript[\[Beta], r]/Subscript[\[Beta], \[Phi]] find p*)


Options[KerrGeoOrbitRPhiResonantP]={PrecisionGoal->Automatic};


KerrGeoOrbitRPhiResonantP[a_?NumericQ, e_?NumericQ, x_?NumericQ, {\[Beta]r_Integer, \[Beta]\[Phi]_Integer}, opts:OptionsPattern[]]:=
Module[{pg,ratio,argpg,resonantEqn,pStar,pp,pgTest, pSep},
	(*See if the user specified some precision goal *)
	If[Not@ValidAQ[a]||Not@ValidEQ[e]||Not@ValidXQ[x]||Not@ValidResIntQ[\[Beta]r,0,\[Beta]\[Phi]],Abort[]];
	
	pg=OptionValue[PrecisionGoal];
	ratio=\[Beta]r/\[Beta]\[Phi];
	If[x>0&&ratio<0, Message[KerrGeoFindResonance::invalidPrograde]; Abort[]];
	If[x<0&&ratio>0, Message[KerrGeoFindResonance::invalidRetrograde]; Abort[]];
	(* pStar provides an initial guess for p. Note that p = pStar for e = 0 and a = 0*)
	pStar=6/(1-ratio^2);
	pSep=KerrGeoSeparatrix[a,e,x];
	
	(* Resonant condition defined by the equation below *)
	If[x==0,
		If[ratio>0,
			resonantEqn[p_?NumericQ]:= r\[Phi]Ratio[a, p, e, x][[1]]-Abs[ratio],
			resonantEqn[p_?NumericQ]:= r\[Phi]Ratio[a, p, e, x][[2]]-Abs[ratio]
		],
		resonantEqn[p_?NumericQ]:= r\[Phi]Ratio[a, p, e, x]-Abs[ratio]
	];
	
	(* Working precision of the root-finding method is based on the precision specified
	 by the PrecisionGoal option, or the precision of the arguments. MachinePrecision
	 is the default precision if not other precision specifications are made. *)
	argpg=Precision[{resonantEqn[pStar],a,e,x,ratio}];
	If[argpg==$MachinePrecision,pg=argpg];
	If[(Not@NumericQ[pg]||pg>argpg),pg=argpg];
	If[pg==Infinity,pg=$MachinePrecision];
	
	If[pg==$MachinePrecision,
		ResonanceSolver[resonantEqn,{pp, pStar, pSep}],
		ResonanceSolver[resonantEqn,{pp, pStar, pSep, PrecisionGoal->pg, WorkingPrecision->pg}],
		ResonanceSolver[resonantEqn,{pp, pStar, pSep, PrecisionGoal->pg, WorkingPrecision->pg}]
	]
];


(* ::Subsubsection::Closed:: *)
(*Given (a,p,x) and Subscript[\[CapitalOmega], r]/Subscript[\[CapitalOmega], \[Phi]]= Subscript[\[Beta], r]/Subscript[\[Beta], \[Phi]] find e*)


Options[KerrGeoOrbitRPhiResonantE]={PrecisionGoal->Automatic};


KerrGeoOrbitRPhiResonantE[a_?NumericQ, p_?NumericQ, x_?NumericQ, {\[Beta]r_Integer, \[Beta]\[Phi]_Integer}, opts:OptionsPattern[]]:=
Module[{pg,argpg,resonantEqn,e0Test,e1Test,eGuess,ee,ratio},
	(*See if the user specified some precision goal *)
	If[Not@ValidAQ[a]||Not@ValidXQ[x]||Not@ValidResIntQ[\[Beta]r,0,\[Beta]\[Phi]],Abort[]];
	
	pg=OptionValue[PrecisionGoal];
	ratio=\[Beta]r/\[Beta]\[Phi];
	If[x>0&&ratio<0, Message[KerrGeoFindResonance::invalidPrograde]; Abort[]];
	If[x<0&&ratio>0, Message[KerrGeoFindResonance::invalidRetrograde]; Abort[]];
	
	(* Test to see if there is an eccentricity that will lead to a bound orbital resonance 
	 based on the provided values of a, p, x *)
	e0Test=KerrGeoOrbitRPhiResonantP[a,0,x,{\[Beta]r,\[Beta]\[Phi]},opts];
	e1Test=KerrGeoOrbitRPhiResonantP[a,1,x,{\[Beta]r,\[Beta]\[Phi]},opts];
	If[p<e0Test||e1Test<p, Message[KerrGeoFindResonance::noresonance, e0Test, e1Test]; Abort[];];
	If[p==e0Test,Return[0]];
	If[p==e1Test,Return[1]];
	eGuess=(p-e0Test)/(e1Test-e0Test);

	(* Resonant condition defined by the equation below *)
	If[x==0,
		If[ratio>0,
			resonantEqn[e_?NumericQ]:= r\[Phi]Ratio[a, p, e, x][[1]]-Abs[ratio],
			resonantEqn[e_?NumericQ]:= r\[Phi]Ratio[a, p, e, x][[2]]-Abs[ratio]
		],
		resonantEqn[e_?NumericQ]:= r\[Phi]Ratio[a, p, e, x]-Abs[ratio]
	];
	
	(* Working precision of the root-finding method is based on the precision specified
	 by the PrecisionGoal option, or the precision of the arguments. MachinePrecision
	 is the default precision if not other precision specifications are made. *)
	argpg=Precision[{resonantEqn[eGuess],a,p,x,ratio}];
	If[argpg==$MachinePrecision,pg=argpg];
	If[(Not@NumericQ[pg]||pg>argpg),pg=argpg];
	If[pg==Infinity,pg=$MachinePrecision];

	If[pg==$MachinePrecision,
		Re[ee/.FindRoot[resonantEqn[ee],{ee,eGuess}]],
		Re[ee/.FindRoot[resonantEqn[ee],{ee,eGuess},PrecisionGoal->pg,WorkingPrecision->pg]],
		Re[ee/.FindRoot[resonantEqn[ee],{ee,eGuess},PrecisionGoal->pg,WorkingPrecision->pg]]
	]
];


(* ::Subsubsection::Closed:: *)
(*Given (a,p,e) and Subscript[\[CapitalOmega], r]/Subscript[\[CapitalOmega], \[Phi]]= Subscript[\[Beta], r]/Subscript[\[Beta], \[Phi]] find x*)


Options[KerrGeoOrbitRPhiResonantX]={PrecisionGoal->Automatic};


KerrGeoOrbitRPhiResonantX[a_?NumericQ, p_?NumericQ, e_?NumericQ, {\[Beta]r_Integer, \[Beta]\[Phi]_Integer}, opts:OptionsPattern[]]:=
Module[{pg,rtt,argpg,resonantEqn,x0Test, x1Test, x0Test1,x0Test2,x1Test1,x1Test2,xGuess,xx1,xx2,xx,ratio},
	(*See if the user specified some precision goal *)
	If[Not@ValidAQ[a]||Not@ValidEQ[e]||Not@ValidResIntQ[\[Beta]r,0,\[Beta]\[Phi]],Abort[]];
	pg=OptionValue[PrecisionGoal];
	ratio=\[Beta]r/\[Beta]\[Phi];
	
	(* Test to see if there is an eccentricity that will lead to a bound orbital resonance 
	based on the provided values of a, p, x *)
	(*See if the user specified some precision goal *)
	If[Not@ValidAQ[a]||Not@ValidEQ[e]||Not@ValidResIntQ[\[Beta]r,\[Beta]\[Theta]],Abort[]];

	pg=OptionValue[PrecisionGoal];
	ratio=\[Beta]r/\[Beta]\[Phi];
	(* Test to see if there is an inclination angle that will lead to a bound orbital resonance 
	 based on the provided values of a, p, e *)
	x0Test=KerrGeoOrbitRPhiResonantP[a,e,0,{\[Beta]r,\[Beta]\[Phi]},opts];
	x1Test=KerrGeoOrbitRPhiResonantP[a,e,Sign[ratio],{\[Beta]r,\[Beta]\[Phi]},opts];
	
	If[ratio>0 && (p<x1Test||x0Test<p), Message[KerrGeoFindResonance::noresonance, x1Test, x0Test]; Abort[];];
	If[ratio<0 && (p<x0Test||x1Test<p), Message[KerrGeoFindResonance::noresonance, x0Test, x1Test]; Abort[];];
	If[p==x0Test,Return[0]];
	If[p==x1Test,Return[Sign[ratio]]];
	xGuess=Sign[ratio](p-x0Test)/(x1Test-x0Test);

	(* Resonant condition defined by the equation below *)
	resonantEqn[x_?NumericQ]:=r\[Phi]Ratio[a, p, e, x]-Abs[ratio];
	
	(* Working precision of the root-finding method is based on the precision specified
	 by the PrecisionGoal option, or the precision of the arguments. MachinePrecision
	 is the default precision if not other precision specifications are made. *)
	argpg=Precision[{a,p,e,ratio}];
	If[argpg==$MachinePrecision,pg=0.9argpg];
	If[(Not@NumericQ[pg]||pg>argpg),pg=0.9argpg];
	If[pg==Infinity,pg=$MachinePrecision];

	If[pg==$MachinePrecision,
		xxx=Re[xx/.FindRoot[resonantEqn[xx],{xx,xGuess}]],
		xxx=Re[xx/.FindRoot[resonantEqn[xx],{xx,xGuess},PrecisionGoal->pg,WorkingPrecision->0.95argpg]],
		xxx=Re[xx/.FindRoot[resonantEqn[xx],{xx,xGuess},PrecisionGoal->pg,WorkingPrecision->0.95argpg]]
	]
];


(* ::Subsection:: *)
(*\[Phi]\[Theta]-resonance*)


(* ::Subsubsection::Closed:: *)
(*Testing functions*)


(* ::Text:: *)
(*Same as those of r\[Theta]-resonance*)


(* ::Subsubsection::Closed:: *)
(*Useful functions for root-finding*)


\[Phi]\[Theta]Ratio[a_, p_, e_, x_/;x!=0]:=term\[Theta][a,p,e,x]+Sign[x]r\[Theta]Ratio[a,p,e,x]termR[a,p,e,x];
\[Phi]\[Theta]Ratio[a_, p_, e_, x_/;x==0]:={1+r\[Theta]Ratio[a,p,e,x]termR[a,p,e,x],1-r\[Theta]Ratio[a,p,e,x]termR[a,p,e,x]};


(* ::Subsubsection::Closed:: *)
(*Given (a,e,x) and Subscript[\[CapitalOmega], \[Phi]]/Subscript[\[CapitalOmega], \[Theta]]= Subscript[\[Beta], \[Phi]]/Subscript[\[Beta], \[Theta]] find p*)


Options[KerrGeoOrbitPhiThetaResonantP]={PrecisionGoal->Automatic};


KerrGeoOrbitPhiThetaResonantP[a_?NumericQ, e_?NumericQ, x_?NumericQ, {\[Beta]\[Theta]_Integer, \[Beta]\[Phi]_Integer}, opts:OptionsPattern[]]:=
Module[{pg,ratio,argpg,resonantEqn,pStar,pp,pgTest,pSep},
	(*See if the user specified some precision goal *)
	If[Not@ValidAQ[a]||Not@ValidEQ[e]||Not@ValidXQ[x]||Not@ValidResIntQ[0,\[Beta]\[Theta],\[Beta]\[Phi]],Abort[]];
	
	pg=OptionValue[PrecisionGoal];
	ratio=\[Beta]\[Phi]/\[Beta]\[Theta];
	If[x>0&&ratio<1, Message[KerrGeoFindResonance::invalidPrograde]; Abort[]];
	If[x<0&&ratio>1, Message[KerrGeoFindResonance::invalidRetrograde]; Abort[]];
	(* pStar provides an initial guess for p. Note that p = pStar for e = 0*)
	pStar=(2 a/Abs[ratio-1])^(2/3);
	pSep=KerrGeoSeparatrix[a,e,x];
	(* Resonant condition defined by the equation below *)
	If[x==0,
		If[ratio>1,
			resonantEqn[p_?NumericQ]:= \[Phi]\[Theta]Ratio[a, p, e, x][[1]]-Abs[ratio],
			resonantEqn[p_?NumericQ]:= \[Phi]\[Theta]Ratio[a, p, e, x][[2]]-Abs[ratio]
		],
		resonantEqn[p_?NumericQ]:= \[Phi]\[Theta]Ratio[a, p, e, x]-Abs[ratio]
	];
	
	(* Working precision of the root-finding method is based on the precision specified
	 by the PrecisionGoal option, or the precision of the arguments. MachinePrecision
	 is the default precision if not other precision specifications are made. *)
	argpg=Precision[{resonantEqn[pStar],a,e,x,ratio}];
	If[argpg==$MachinePrecision,pg=argpg];
	If[(Not@NumericQ[pg]||pg>argpg),pg=argpg];
	If[pg==Infinity,pg=$MachinePrecision];
	
	If[pg==$MachinePrecision,
		ResonanceSolver[resonantEqn,{pp, pStar, pSep}],
		ResonanceSolver[resonantEqn,{pp, pStar, pSep, PrecisionGoal->pg, WorkingPrecision->pg}],
		ResonanceSolver[resonantEqn,{pp, pStar, pSep, PrecisionGoal->pg, WorkingPrecision->pg}]
	]
];


(* ::Subsubsection::Closed:: *)
(*Given (a,p,x) and Subscript[\[CapitalOmega], \[Phi]]/Subscript[\[CapitalOmega], \[Theta]]= Subscript[\[Beta], \[Phi]]/Subscript[\[Beta], \[Theta]] find e*)


Options[KerrGeoOrbitPhiThetaResonantE]={PrecisionGoal->Automatic};


KerrGeoOrbitPhiThetaResonantE[a_?NumericQ, p_?NumericQ, x_?NumericQ, {\[Beta]\[Theta]_Integer, \[Beta]\[Phi]_Integer}, opts:OptionsPattern[]]:=
Module[{pg,argpg,resonantEqn,e0Test,e1Test,eGuess,ee,ratio},
	(*See if the user specified some precision goal *)
	If[Not@ValidAQ[a]||Not@ValidXQ[x]||Not@ValidResIntQ[0,\[Beta]\[Theta],\[Beta]\[Phi]],Abort[]];
	
	pg=OptionValue[PrecisionGoal];
	ratio=\[Beta]\[Phi]/\[Beta]\[Theta];
	If[x>0&&ratio<1, Message[KerrGeoFindResonance::invalidPrograde]; Abort[]];
	If[x<0&&ratio>1, Message[KerrGeoFindResonance::invalidRetrograde]; Abort[]];
	
	(* Test to see if there is an eccentricity that will lead to a bound orbital resonance 
	 based on the provided values of a, p, x *)
	e0Test=KerrGeoOrbitPhiThetaResonantP[a,0,x,{\[Beta]\[Theta],\[Beta]\[Phi]},opts];
	e1Test=KerrGeoOrbitPhiThetaResonantP[a,1,x,{\[Beta]\[Theta],\[Beta]\[Phi]},opts];
	If[p<e0Test||e1Test<p, Message[KerrGeoFindResonance::noresonance, e0Test, e1Test]; Abort[];];
	If[p==e0Test,Return[0]];
	If[p==e1Test,Return[1]];
	eGuess=(p-e0Test)/(e1Test-e0Test);

	(* Resonant condition defined by the equation below *)
	If[x==0,
		If[ratio>1,
			resonantEqn[e_?NumericQ]:= \[Phi]\[Theta]Ratio[a, p, e, x][[1]]-Abs[ratio],
			resonantEqn[e_?NumericQ]:= \[Phi]\[Theta]Ratio[a, p, e, x][[2]]-Abs[ratio]
		],
		resonantEqn[e_?NumericQ]:= \[Phi]\[Theta]Ratio[a, p, e, x]-Abs[ratio]
	];
	
	(* Working precision of the root-finding method is based on the precision specified
	 by the PrecisionGoal option, or the precision of the arguments. MachinePrecision
	 is the default precision if not other precision specifications are made. *)
	argpg=Precision[{resonantEqn[eGuess],a,p,x,ratio}];
	If[argpg==$MachinePrecision,pg=argpg];
	If[(Not@NumericQ[pg]||pg>argpg),pg=argpg];
	If[pg==Infinity,pg=$MachinePrecision];

	If[pg==$MachinePrecision,
		Re[ee/.FindRoot[resonantEqn[ee],{ee,eGuess}]],
		Re[ee/.FindRoot[resonantEqn[ee],{ee,eGuess},PrecisionGoal->pg,WorkingPrecision->pg]],
		Re[ee/.FindRoot[resonantEqn[ee],{ee,eGuess},PrecisionGoal->pg,WorkingPrecision->pg]]
	]
];


(* ::Subsubsection::Closed:: *)
(*Given (a,p,e) and Subscript[\[CapitalOmega], \[Theta]]/Subscript[\[CapitalOmega], \[Phi]]= Subscript[\[Beta], \[Theta]]/Subscript[\[Beta], \[Phi]] find x*)


Options[KerrGeoOrbitPhiThetaResonantX]={PrecisionGoal->Automatic};


KerrGeoOrbitPhiThetaResonantX[a_?NumericQ, p_?NumericQ, e_?NumericQ, {\[Beta]\[Theta]_Integer, \[Beta]\[Phi]_Integer}, opts:OptionsPattern[]]:=
Module[{pg,argpg,resonantEqn,x0Test, x1Test,xGuess,xx,ratio},
	(* Test to see if there is an eccentricity that will lead to a bound orbital resonance 
	based on the provided values of a, p, x *)
	(*See if the user specified some precision goal *)
	If[Not@ValidAQ[a]||Not@ValidEQ[e]||Not@ValidResIntQ[0,\[Beta]\[Theta],\[Beta]\[Phi]],Abort[]];

	pg=OptionValue[PrecisionGoal];
	ratio=\[Beta]\[Phi]/\[Beta]\[Theta];
	(* Test to see if there is an inclination angle that will lead to a bound orbital resonance 
	 based on the provided values of a, p, e *)
	x0Test=KerrGeoOrbitPhiThetaResonantP[a,e,0,{\[Beta]\[Theta],\[Beta]\[Phi]},opts];
	x1Test=KerrGeoOrbitPhiThetaResonantP[a,e,Sign[ratio-1],{\[Beta]\[Theta],\[Beta]\[Phi]},opts];
	
	If[ratio>1 && (p<x1Test||x0Test<p), Message[KerrGeoFindResonance::noresonance, x1Test, x0Test]; Abort[];];
	If[ratio<1 && (p<x0Test||x1Test<p), Message[KerrGeoFindResonance::noresonance, x0Test, x1Test]; Abort[];];
	If[p==x0Test,Return[0]];
	If[p==x1Test,Return[Sign[ratio-1]]];
	xGuess=Sign[ratio-1](p-x0Test)/(x1Test-x0Test);

	(* Resonant condition defined by the equation below *)
	resonantEqn[x_?NumericQ]:=\[Phi]\[Theta]Ratio[a, p, e, x]-Abs[ratio];
	
	(* Working precision of the root-finding method is based on the precision specified
	 by the PrecisionGoal option, or the precision of the arguments. MachinePrecision
	 is the default precision if not other precision specifications are made. *)
	argpg=Precision[{a,p,e,ratio}];
	If[argpg==$MachinePrecision,pg=0.9argpg];
	If[(Not@NumericQ[pg]||pg>argpg),pg=0.9argpg];
	If[pg==Infinity,pg=$MachinePrecision];

	If[pg==$MachinePrecision,
		Re[xx/.FindRoot[resonantEqn[xx],{xx,xGuess}]],
		Re[xx/.FindRoot[resonantEqn[xx],{xx,xGuess},PrecisionGoal->pg,WorkingPrecision->0.95argpg]],
		Re[xx/.FindRoot[resonantEqn[xx],{xx,xGuess},PrecisionGoal->pg,WorkingPrecision->0.95argpg]]
	]
];


(* ::Subsection:: *)
(*r\[Theta]\[Phi]-resonance*)


(* ::Subsubsection::Closed:: *)
(*Testing functions*)


(* ::Text:: *)
(*Same as those of r\[Theta]-resonance*)


(* ::Subsubsection::Closed:: *)
(*Useful functions for root-finding*)


(* ::Text:: *)
(*Use those of double resonances*)


(* ::Subsubsection:: *)
(*Given (a,e,\[Beta]r,\[Beta]\[Theta]) find (\[Beta]rmin, \[Beta]rmax)*)


Options[KerrGeoOrbitTripleResonantPX\[Beta]r]={PrecisionGoal->Automatic};


KerrGeoOrbitTripleResonantPX\[Beta]r[a_?NumericQ, e_?NumericQ, {\[Beta]\[Theta]_Integer, \[Beta]\[Phi]_Integer}, opts:OptionsPattern[]]:=
Module[{pg,argpg,resonantEqnr\[Theta],resonantEqnr\[Phi],p0Test,p1Test,r\[Theta]ratio0,r\[Theta]ratio1,\[Phi]\[Theta]ratio},
	If[Not@ValidAQ[a]||Not@ValidEQ[e]||Not@ValidXQ[x]||Not@ValidResIntQ[0,\[Beta]\[Theta],\[Beta]\[Phi]],Abort[]];
	\[Phi]\[Theta]ratio=\[Beta]\[Phi]/\[Beta]\[Theta];
	p0Test=KerrGeoOrbitPhiThetaResonantP[a, e, 0, {\[Beta]\[Theta], \[Beta]\[Phi]}, opts];
	p1Test=KerrGeoOrbitPhiThetaResonantP[a, e, Sign[\[Phi]\[Theta]ratio-1], {\[Beta]\[Theta], \[Beta]\[Phi]}, opts];
	r\[Theta]ratio0=Sign[\[Phi]\[Theta]ratio-1]r\[Theta]Ratio[a, p0Test, e, 0];
	r\[Theta]ratio1=Sign[\[Phi]\[Theta]ratio-1]r\[Theta]Ratio[a, p1Test, e, Sign[\[Phi]\[Theta]ratio-1]];
	Sort[{\[Beta]\[Theta] r\[Theta]ratio0, \[Beta]\[Theta] r\[Theta]ratio1}, Less]
];


(* ::Subsubsection::Closed:: *)
(*Given (a,e,\[Beta]r,\[Beta]\[Phi]) find (\[Beta]\[Theta]min, \[Beta]\[Theta]max)*)


Options[KerrGeoOrbitTripleResonantPX\[Beta]\[Theta]]={PrecisionGoal->Automatic};


KerrGeoOrbitTripleResonantPX\[Beta]\[Theta][a_?NumericQ, e_?NumericQ, {\[Beta]r_Integer, \[Beta]\[Phi]_Integer}, opts:OptionsPattern[]]:=
Module[{pg,argpg,resonantEqnr\[Theta],resonantEqnr\[Phi],p0Test,p1Test,r\[Theta]ratio0,r\[Theta]ratio1,r\[Phi]ratio},
	If[Not@ValidAQ[a]||Not@ValidEQ[e]||Not@ValidXQ[x]||Not@ValidResIntQ[\[Beta]r,0,\[Beta]\[Phi]],Abort[]];
	r\[Phi]ratio=\[Beta]r/\[Beta]\[Phi];
	p0Test=KerrGeoOrbitRPhiResonantP[a, e, 0, {\[Beta]r, \[Beta]\[Phi]}, opts];
	p1Test=KerrGeoOrbitRPhiResonantP[a, e, Sign[r\[Phi]ratio], {\[Beta]r, \[Beta]\[Phi]}, opts];
	r\[Theta]ratio0=Sign[r\[Phi]ratio]r\[Theta]Ratio[a, p0Test, e, 0];
	r\[Theta]ratio1=Sign[r\[Phi]ratio]r\[Theta]Ratio[a, p1Test, e, Sign[r\[Phi]ratio]];
	Sort[{\[Beta]r/r\[Theta]ratio0, \[Beta]r/r\[Theta]ratio1}, Less]
];


(* ::Subsubsection::Closed:: *)
(*Given (a,e,\[Beta]r,\[Beta]\[Theta]) find (\[Beta]\[Phi]min, \[Beta]\[Phi]max)*)


Options[KerrGeoOrbitTripleResonantPX\[Beta]\[Phi]]={PrecisionGoal->Automatic};


KerrGeoOrbitTripleResonantPX\[Beta]\[Phi][a_?NumericQ, e_?NumericQ, {\[Beta]r_Integer, \[Beta]\[Theta]_Integer}, opts:OptionsPattern[]]:=
Module[{pg,argpg,resonantEqnr\[Theta],resonantEqnr\[Phi],p0Test,p1Test,r\[Phi]ratio0,r\[Phi]ratio1,r\[Theta]ratio},
	If[Not@ValidAQ[a]||Not@ValidEQ[e]||Not@ValidXQ[x]||Not@ValidResIntQ[\[Beta]r,\[Beta]\[Theta],0],Abort[]];
	r\[Theta]ratio=\[Beta]r/\[Beta]\[Theta];
	p0Test=KerrGeoOrbitRThetaResonantP[a, e, 0, {\[Beta]r, \[Beta]\[Theta]}, opts];
	p1Test=KerrGeoOrbitRThetaResonantP[a, e, Sign[r\[Theta]ratio], {\[Beta]r, \[Beta]\[Theta]}, opts];
	r\[Phi]ratio0=Sign[r\[Theta]ratio]If[r\[Theta]ratio>0, r\[Phi]Ratio[a, p0Test, e, 0][[1]], r\[Phi]Ratio[a, p0Test, e, 0][[2]]];
	r\[Phi]ratio1=Sign[r\[Theta]ratio]r\[Phi]Ratio[a, p1Test, e, Sign[r\[Theta]ratio]];
	Sort[{\[Beta]r/r\[Phi]ratio0, \[Beta]r/r\[Phi]ratio1}, Less]
];


(* ::Subsubsection::Closed:: *)
(*Given (a,e,\[Beta]r,\[Beta]\[Theta],\[Beta]\[Phi]) find (p, x)*)


Options[KerrGeoOrbitTripleResonantPX]={PrecisionGoal->Automatic};


KerrGeoOrbitTripleResonantPX[a_?NumericQ, e_?NumericQ, {\[Beta]r_Integer, \[Beta]\[Theta]_Integer, \[Beta]\[Phi]_Integer}, opts:OptionsPattern[]]:=
Module[{pg,argpg,resonantEqnr\[Theta],resonantEqnr\[Phi],pGuess,xGuess,pp,xx,r\[Theta]ratio,r\[Phi]ratio,p0r\[Theta]Test,p1r\[Theta]Test,p0r\[Phi]Test,p1r\[Phi]Test,p0\[Phi]\[Theta]Test,p1\[Phi]\[Theta]Test,checkSolution, rInts, \[Theta]Ints, \[Phi]Ints},
	(*See if the user specified some precision goal *)
	If[Not@ValidAQ[a]||Not@ValidEQ[e]||Not@ValidXQ[x]||Not@ValidResIntQ[0,\[Beta]\[Theta],\[Beta]\[Phi]],Abort[]];
	
	pg=OptionValue[PrecisionGoal];
	r\[Theta]ratio=\[Beta]r/\[Beta]\[Theta];
	r\[Phi]ratio=\[Beta]r/\[Beta]\[Phi];
	p0r\[Theta]Test=KerrGeoOrbitRThetaResonantP[a, e, 0, {\[Beta]r, \[Beta]\[Theta]},opts];
	p1r\[Theta]Test=KerrGeoOrbitRThetaResonantP[a, e, Sign[r\[Theta]ratio], {\[Beta]r, \[Beta]\[Theta]},opts];
	p0r\[Phi]Test=KerrGeoOrbitRPhiResonantP[a, e, 0, {\[Beta]r, \[Beta]\[Phi]},opts];
	p1r\[Phi]Test=KerrGeoOrbitRPhiResonantP[a, e, Sign[r\[Phi]ratio], {\[Beta]r, \[Beta]\[Phi]},opts];
	checkSolution=(p0r\[Theta]Test-p0r\[Phi]Test)(p1r\[Theta]Test-p1r\[Phi]Test);
	If[checkSolution>0,
		rInts=KerrGeoOrbitTripleResonantPX\[Beta]r[a, e, {\[Beta]\[Theta], \[Beta]\[Phi]}, opts];
		\[Theta]Ints=KerrGeoOrbitTripleResonantPX\[Beta]\[Theta][a, e, {\[Beta]r, \[Beta]\[Phi]}, opts];
		\[Phi]Ints=KerrGeoOrbitTripleResonantPX\[Beta]\[Phi][a, e, {\[Beta]r, \[Beta]\[Theta]}, opts];
		Message[KerrGeoFindResonance::noTripleResonance,rInts[[1]],rInts[[2]],\[Theta]Ints[[1]],\[Theta]Ints[[2]],\[Phi]Ints[[1]],\[Phi]Ints[[2]]];
		Abort[]];
	
	xGuess=Sign[r\[Theta]ratio]/((p1r\[Phi]Test-p1r\[Theta]Test)/(p0r\[Theta]Test-p0r\[Phi]Test)+1);
	pGuess=p0r\[Theta]Test-Sign[r\[Theta]ratio]xGuess(p0r\[Theta]Test-p1r\[Theta]Test);
	resonantEqnr\[Theta][p_?NumericQ, x_?NumericQ]:=r\[Theta]Ratio[a, p, e, x]-Abs[r\[Theta]ratio];
	resonantEqnr\[Phi][p_?NumericQ, x_?NumericQ]:=r\[Phi]Ratio[a, p, e, x]-Abs[r\[Phi]ratio];
	argpg=Precision[{a,e,r\[Theta]ratio,r\[Phi]ratio,resonantEqnr\[Theta][pGuess, xGuess],resonantEqnr\[Phi][pGuess, xGuess]}];
	If[argpg==$MachinePrecision,pg=0.9argpg];
	If[(Not@NumericQ[pg]||pg>argpg),pg=0.9argpg];
	If[pg==Infinity,pg=$MachinePrecision];
	If[pg==$MachinePrecision,
		Re[{pp,xx}/.FindRoot[{resonantEqnr\[Theta][pp,xx],resonantEqnr\[Phi][pp,xx]},{{pp,pGuess},{xx,xGuess}}]],
		Re[{pp,xx}/.FindRoot[{resonantEqnr\[Theta][pp,xx],resonantEqnr\[Phi][pp,xx]},{{pp,pGuess},{xx,xGuess}},PrecisionGoal->pg,WorkingPrecision->0.95argpg]],
		Re[{pp,xx}/.FindRoot[{resonantEqnr\[Theta][pp,xx],resonantEqnr\[Phi][pp,xx]},{{pp,pGuess},{xx,xGuess}},PrecisionGoal->pg,WorkingPrecision->0.95argpg]]
	]
];


(* ::Subsection::Closed:: *)
(*Generic resonance interface*)


Options[KerrGeoFindResonance]={PrecisionGoal->Automatic};


KerrGeoFindResonance[assoc_Association,{\[Beta]r_Integer, \[Beta]\[Theta]_Integer, \[Beta]\[Phi]_Integer},opts:OptionsPattern[]]:=Module[{pp, xx},
	If[ContainsExactly[Keys[assoc],{"a","e"}],
		{pp, xx}=KerrGeoOrbitTripleResonantPX["a"/.assoc, "e"/.assoc, {\[Beta]r, \[Beta]\[Theta], \[Beta]\[Phi]}, opts];
		{"p"->pp,"x"->xx},
		Message[KerrGeoFindResonance::assocErrTriple]
	]
]

KerrGeoFindResonance[assoc_Association,{\[Beta]r_Integer, \[Beta]\[Theta]_Integer, 0},opts:OptionsPattern[]]:=Block[{},
	If[ContainsExactly[Keys[assoc],{"a","p","x"}],
		"e"->KerrGeoOrbitRThetaResonantE["a"/.assoc, "p"/.assoc, "x"/.assoc, {\[Beta]r, \[Beta]\[Theta]}, opts],
		If[ContainsExactly[Keys[assoc],{"a","e","x"}],
			"p"->KerrGeoOrbitRThetaResonantP["a"/.assoc, "e"/.assoc, "x"/.assoc, {\[Beta]r, \[Beta]\[Theta]}, opts],
			If[ContainsExactly[Keys[assoc],{"a","p","e"}],
				"x"->KerrGeoOrbitRThetaResonantX["a"/.assoc, "p"/.assoc, "e"/.assoc, {\[Beta]r, \[Beta]\[Theta]}, opts],
				Message[KerrGeoFindResonance::assocErr]
			]
		]
	]
]
KerrGeoFindResonance[assoc_Association,{\[Beta]r_Integer, 0, \[Beta]\[Phi]_Integer},opts:OptionsPattern[]]:=Block[{},
	If[ContainsExactly[Keys[assoc],{"a","p","x"}],
		"e"->KerrGeoOrbitRPhiResonantE["a"/.assoc, "p"/.assoc, "x"/.assoc, {\[Beta]r, \[Beta]\[Phi]}, opts],
		If[ContainsExactly[Keys[assoc],{"a","e","x"}],
			"p"->KerrGeoOrbitRPhiResonantP["a"/.assoc, "e"/.assoc, "x"/.assoc, {\[Beta]r, \[Beta]\[Phi]}, opts],
			If[ContainsExactly[Keys[assoc],{"a","p","e"}],
				"x"->KerrGeoOrbitRPhiResonantX["a"/.assoc, "p"/.assoc, "e"/.assoc, {\[Beta]r, \[Beta]\[Phi]}, opts],
				Message[KerrGeoFindResonance::assocErr]
			]
		]
	]
]

KerrGeoFindResonance[assoc_Association,{0, \[Beta]\[Theta]_Integer, \[Beta]\[Phi]_Integer},opts:OptionsPattern[]]:=Block[{},
	If[ContainsExactly[Keys[assoc],{"a","p","x"}],
		"e"->KerrGeoOrbitPhiThetaResonantE["a"/.assoc, "p"/.assoc, "x"/.assoc, {\[Beta]\[Theta], \[Beta]\[Phi]}, opts],
		If[ContainsExactly[Keys[assoc],{"a","e","x"}],
			"p"->KerrGeoOrbitPhiThetaResonantP["a"/.assoc, "e"/.assoc, "x"/.assoc, {\[Beta]\[Theta], \[Beta]\[Phi]}, opts],
			If[ContainsExactly[Keys[assoc],{"a","p","e"}],
				"x"->KerrGeoOrbitPhiThetaResonantX["a"/.assoc, "p"/.assoc, "e"/.assoc, {\[Beta]\[Theta], \[Beta]\[Phi]}, opts],
				Message[KerrGeoFindResonance::assocErr]
			]
		]
	]
]


(* ::Section::Closed:: *)
(*Close the package*)


End[];

EndPackage[];
