(* ::Package:: *)

(* Author: Benjamin Suzzoni
Organization: University of Southampton *)

(* This package contains definitions for various curvature invariants of General Relativity and custom implementations of exterior calculus *)
BeginPackage["GRBenterre`"]

Christoffel::usage = "Christoffel[g, xx] Calculates the Christoffel symbols for the associated metric g and coordinates xx.
Conventions are such that its components are given by spacetime indices \!\(\*SubscriptBox[SuperscriptBox[\[CapitalGamma],\[Mu]], \[Nu]\[Rho]]\) in that same order.
These are computed from the metric using

\!\(\*SubscriptBox[SuperscriptBox[\[CapitalGamma],\[Mu]], \[Nu]\[Rho]]=\*FractionBox[\(1\), \(2\)]\\ \*SuperscriptBox[\(g\),\[Mu]\[Sigma]]\\ ( \*SubscriptBox[\(\[PartialD]\), \(\[Nu]\)]\*SubscriptBox[\(g\),\[Rho]\[Sigma]]+\*SubscriptBox[\(\[PartialD]\), \(\[Rho]\)]\*SubscriptBox[\(g\),\[Nu]\[Sigma]]-\*SubscriptBox[\(\[PartialD]\), \(\[Sigma]\)]\*SubscriptBox[\(g\),\[Nu]\[Rho]]) \)"

RiemannTensor::usage = "RiemannTensor[g, xx] Calculates the Riemann tensor for the associated metric g and coordinates xx.
Conventions are such that all spacetime indices are down \!\(\*SubscriptBox[\(R\), \[Mu]\[Nu]\[Rho]\[Sigma]]\).
The components are computed from the metric using

\!\(\*SubscriptBox[\(R\), \[Mu]\[Nu]\[Rho]\[Sigma]]=\*SubscriptBox[\(g\), \[Mu]\[Alpha]]\\ (\*SubscriptBox[\(\[PartialD]\), \(\[Rho]\)]\*SubscriptBox[SuperscriptBox[\[CapitalGamma],\[Alpha]], \[Sigma]\[Nu]] - \*SubscriptBox[\(\[PartialD]\), \(\[Sigma]\)]\*SubscriptBox[SuperscriptBox[\[CapitalGamma],\[Alpha]], \[Rho]\[Nu]] + \*SubscriptBox[SuperscriptBox[\[CapitalGamma],\[Alpha]], \[Rho]\[Lambda]]\*SubscriptBox[SuperscriptBox[\[CapitalGamma],\[Lambda]], \[Sigma]\[Nu]] - \*SubscriptBox[SuperscriptBox[\[CapitalGamma],\[Alpha]], \[Sigma]\[Lambda]]\*SubscriptBox[SuperscriptBox[\[CapitalGamma],\[Lambda]], \[Rho]\[Nu]] ) \)"

RicciTensor::usage =  "Calculates the Ricci tensor for the associated metric g and coordinates xx.
Conventions are such that all spacetime indices are down \!\(\*SubscriptBox[\(R\), \[Mu]\[Nu]]\).
The components are computed from the metric using

\!\(\*SubscriptBox[\(R\), \[Mu]\[Nu]]=\*SuperscriptBox[\(g\), \[Rho]\[Sigma]]\\ \*SubscriptBox[\(R\), \[Rho]\[Mu]\[Sigma]\[Nu]] \)"

RicciScalar::usage =  "Calculates the Ricci scalar for the associated metric g and coordinates xx.
The components are computed from the metric using

\!\( \(R\)=\*SuperscriptBox[\(g\), \[Mu]\[Nu]]\\ \*SubscriptBox[\(R\), \[Mu]\[Nu]] \)"

KretschmannScalar::usage =  "Calculates the Kretschmann scalar for the associated metric g and coordinates xx.
Its expression is identical to the third Lovelock invariant at level 2, I32.
The components are computed from the metric using

\!\( \*SubscriptBox[\(K\),\(1\)] = \*SubscriptBox[\(R\),\[Mu]\[Nu]\[Rho]\[Sigma]]\*SuperscriptBox[\(R\),\[Mu]\[Nu]\[Rho]\[Sigma]] \)"

WeylTensor::usage =   "Calculates the Weyl tensor for the associated metric g and coordinates xx.
Conventions are such that all spacetime indices are down \!\(\*SubscriptBox[\(C\), \[Mu]\[Nu]\[Rho]\[Sigma]]\).
The components are computed from the metric using

\!\( \*SubscriptBox[\(C\), \[Mu]\[Nu]\[Rho]\[Sigma]] = \*SubscriptBox[\(R\), \[Mu]\[Nu]\[Rho]\[Sigma]] + \*FractionBox[\(1\), \(d-2\)]\\ (\*SubscriptBox[\(R\), \[Mu]\[Sigma]]\*SubscriptBox[\(g\), \[Nu]\[Rho]] - \*SubscriptBox[\(R\), \[Mu]\[Rho]]\*SubscriptBox[\(g\), \[Nu]\[Sigma]] + \*SubscriptBox[\(R\), \[Nu]\[Rho]]\*SubscriptBox[\(g\), \[Mu]\[Sigma]] - \*SubscriptBox[\(R\), \[Nu]\[Sigma]]\*SubscriptBox[\(g\), \[Mu]\[Rho]]) + \*FractionBox[\(1\), \((d-1)(d-2)\)]\\ R(\*SubscriptBox[\(g\), \[Mu]\[Rho]]\*SubscriptBox[\(g\), \[Nu]\[Sigma]] - \*SubscriptBox[\(g\), \[Mu]\[Sigma]]\*SubscriptBox[\(g\), \[Nu]\[Rho]]) \)"

EulerScalar::usage =  "Calculates the Euler scalar for the associated metric g and coordinates xx.
The components are computed from the metric through the Weyl tensor as

\!\( \*SubscriptBox[\(K\),\(3\)] = -\*SubscriptBox[\(C\),\[Mu]\[Nu]\[Rho]\[Sigma]]\*SuperscriptBox[\(C\),\[Mu]\[Nu]\[Rho]\[Sigma]] + 2\*SubscriptBox[\(R\),\[Mu]\[Nu]]\*SuperscriptBox[\(R\),\[Mu]\[Nu]] - \*FractionBox[\(2\), \(3\)]\*SuperscriptBox[\(R\),\(2\)]\)"

CottonTensor::usage =  "Calculates the Cotton tensor for the associated metric g and coordinates xx.
The components are computed from the metric using

\!\( \*SubscriptBox[\(C\),\[Mu]\[Nu]\[Rho]] = \*SubscriptBox[\(\[Del]\),\[Rho]]\*SubscriptBox[\(R\),\[Mu]\[Nu]] - \*SubscriptBox[\(\[Del]\),\[Nu]]\*SubscriptBox[\(R\),\[Mu]\[Rho]]  + \*FractionBox[\(1\),\(2(d-1)\)] ( \*SubscriptBox[\(\[Del]\),\[Nu]] R \*SubscriptBox[\(g\),\[Mu]\[Rho]] - \*SubscriptBox[\(\[Del]\),\[Rho]] R \*SubscriptBox[\(g\),\[Mu]\[Nu]])\)"

SchoutenTensor::usage =  "Calculates the Schouten tensor for the associated metric g and coordinates xx.
The components are computed from the metric using

\!\( \*SubscriptBox[\(P\),\[Mu]\[Nu]] = \*FractionBox[\(1\),\(d-2\)] ( \*SubscriptBox[\(R\),\[Mu]\[Nu]] - \*FractionBox[\(R\),\(2(d-1)\)] \*SubscriptBox[\(g\),\[Mu]\[Nu]] )\)"


Lovelock::usage = "Lovelock[k,n][g,xx] returns the nth Lovelock invariant at order k.
For n greater than the number of Lovelock invariant at a given order k, Lovelock might return additional invariants built from the Weyl tensor.

Type Lovelock1, Lovelock2,... for more information about each specific order."

Lovelock1::usage = "The unique Lovelock invariant at order 1 is equal to the Ricci scalar

\!\(\*SubscriptBox[SuperscriptBox[\(Lov\),\(1\)], \(1\)]=R\)"

Lovelock2::usage = "The Lovelock invariant at order 2 are

\!\(\*SubscriptBox[SuperscriptBox[\(Lov\),\(2\)], \(1\)]=\*SuperscriptBox[\(R\),\(2\)]\)
\!\(\*SubscriptBox[SuperscriptBox[\(Lov\),\(2\)], \(2\)]=\*SubscriptBox[\(R\),\[Mu]\[Nu]]\*SuperscriptBox[\(R\),\[Mu]\[Nu]]\)
\!\(\*SubscriptBox[SuperscriptBox[\(Lov\),\(2\)], \(3\)]=\*SubscriptBox[\(R\),\[Mu]\[Nu]\[Rho]\[Sigma]]\*SuperscriptBox[\(R\),\[Mu]\[Nu]\[Rho]\[Sigma]]=\*SubscriptBox[\(K\),\(1\)]\)

Plus an extra one constructed out of the Weyl tensor (which we abusingly also label Lovelock)

\!\(\*SubscriptBox[SuperscriptBox[\(Lov\),\(2\)], \(4\)]=\*SubscriptBox[\(C\),\[Mu]\[Nu]\[Rho]\[Sigma]]\*SuperscriptBox[\(C\),\[Mu]\[Nu]\[Rho]\[Sigma]]\)"

Lovelock3::usage =  "The Lovelock invariants at order 3 are

\!\(\*SubscriptBox[SuperscriptBox[\(Lov\),\(3\)], \(1\)]=\*SuperscriptBox[\(R\),\(3\)]\)
\!\(\*SubscriptBox[SuperscriptBox[\(Lov\),\(3\)], \(2\)]=\(R\)\*SubscriptBox[\(R\),\[Mu]\[Nu]]\*SuperscriptBox[\(R\),\[Mu]\[Nu]]=\(R\) \*SubscriptBox[SuperscriptBox[\(Lov\),\(2\)], \(2\)]\)
\!\(\*SubscriptBox[SuperscriptBox[\(Lov\),\(3\)], \(3\)]=\*SubscriptBox[\(R\),\[Nu]\[Alpha]]\*SubscriptBox[SuperscriptBox[\(R\),\[Nu]],\[Mu]]\*SuperscriptBox[\(R\),\[Alpha]\[Mu]]\)
\!\(\*SubscriptBox[SuperscriptBox[\(Lov\),\(3\)], \(4\)]=\*SubscriptBox[\(R\),\[Nu]\[Alpha]]\*SubscriptBox[\(R\),\[Mu]\[Beta]]\*SuperscriptBox[\(R\),\[Nu]\[Mu]\[Alpha]\[Beta]]\)
\!\(\*SubscriptBox[SuperscriptBox[\(Lov\),\(3\)], \(5\)]=\(R\)\*SubscriptBox[\(R\),\[Mu]\[Nu]\[Alpha]\[Beta]]\*SuperscriptBox[\(R\),\[Mu]\[Nu]\[Alpha]\[Beta]]=\(R\)\*SubscriptBox[\(K\),\(1\)]\)
\!\(\*SubscriptBox[SuperscriptBox[\(Lov\),\(3\)], \(6\)]=\*SubscriptBox[\(R\),\[Nu]\[Alpha]]\*SuperscriptBox[SubscriptBox[\(R\),\[Beta]\[Gamma]\[Epsilon]],\[Nu]]\*SuperscriptBox[\(R\),\[Beta]\[Gamma]\[Epsilon]\[Alpha]]\)
\!\(\*SubscriptBox[SuperscriptBox[\(Lov\),\(3\)], \(7\)]=\*SubscriptBox[\(R\),\[Mu]\[Nu]\[Alpha]\[Beta]]\*SubscriptBox[SuperscriptBox[\(R\),\[Mu]\[Nu]],\[Gamma]\[Epsilon]]\*SuperscriptBox[\(R\),\[Alpha]\[Beta]\[Gamma]\[Epsilon]]\)
\!\(\*SubscriptBox[SuperscriptBox[\(Lov\),\(3\)], \(8\)]=\*SubscriptBox[\(R\),\[Mu]\[Nu]\[Alpha]\[Beta]]\*SubscriptBox[SuperscriptBox[SubscriptBox[SuperscriptBox[\(R\),\[Mu]],\[Gamma]],\[Alpha]],\[Epsilon]]\*SuperscriptBox[\(R\),\[Nu]\[Gamma]\[Beta]\[Epsilon]]\)

Plus an extra three constructed out of the Weyl tensor (which we abusingly also label Lovelock)

\!\(\*SubscriptBox[SuperscriptBox[\(Lov\),\(3\)], \(9\)]=\*SubscriptBox[\(C\),\[Mu]\[Nu]\[Alpha]\[Beta]]\*SuperscriptBox[\(R\),\[Mu]\[Alpha]]\*SuperscriptBox[\(R\),\[Nu]\[Beta]]\)
\!\(\*SubscriptBox[SuperscriptBox[\(Lov\),\(3\)], \(10\)]=\*SubscriptBox[\(C\),\[Mu]\[Nu]\[Alpha]\[Beta]]\*SuperscriptBox[\(C\),\[Alpha]\[Beta]\[Rho]\[Nu]]\*SubscriptBox[SuperscriptBox[\(R\),\[Mu]],\[Rho]]\)
\!\(\*SubscriptBox[SuperscriptBox[\(Lov\),\(3\)], \(11\)]=\*SubscriptBox[\(C\),\[Mu]\[Nu]\[Alpha]\[Beta]]\*SuperscriptBox[\(C\),\[Alpha]\[Beta]\[Rho]\[Sigma]]\*SuperscriptBox[SubscriptBox[\(C\),\[Rho]\[Sigma]],\[Mu]\[Nu]]\)"

LovelockAll::usage = "LovelockAll[g,xx] Returns a table containing all Lovelock Invariant available in this package.

Optional Arguments:
 - \"Extend\"->Bool. If set to True, includes invariants built from the Weyl tensor (defaults to True)."

GRInvariants::usage = "GRInvariants[g,xx] Returns a table containing all the scalar curvature invariant available in this package.

Optional Arguments:
 - \"Minimal\"->Bool. If set to True, only display a basis of multiplicatively independent invariants (defaults to False)."

GRlist::usgae = "The current list of available functions is

----- Standard quantities -----
Christoffel[g, xx]: Christoffel symbols
RiemannTensor[g, xx]: Riemann tensor
RicciTensor[g, xx]: Ricci tensor
RicciScalar[g, xx]: Ricci scalar
KretschmannScalar[g, xx]: Kretschmann scalar
WeylTensor[g, xx]: Weyl tensor
EulerScalar[g, xx]: Euler scalar
CottonTensor[g, xx]: Cotton tensor
SchoutenTensor[g, xx]: Schouten tensor

----- Lovelock invariants -----
Lovelock[1,1][g,xx]

Lovelock[2,1][g,xx]
Lovelock[2,2][g,xx]
Lovelock[2,3][g,xx]

Lovelock[3,1][g, xx]
Lovelock[3,2][g, xx]
Lovelock[3,3][g, xx]
Lovelock[3,4][g, xx]
Lovelock[3,5][g, xx]
Lovelock[3,6][g, xx]
Lovelock[3,7][g, xx]
Lovelock[3,8][g, xx]

----- Mixed Invariants -----
Lovelock[2,4][g,xx]: Weyl tensor squared
Lovelock[3,9][g, xx]
Lovelock[3,10][g, xx]
Lovelock[3,11][g, xx]

LovelockAll[g,xx]: All the Lovelock and mixed invariants listed above
GRInvariants[g,xx]: All the scalar invariants in this package
"

Invariants::usage = " This package includes definitions for various curvature invariants!

All quantities require a metric g, in the form of a dxd array, and an array of d coordinates, xx.

The various conventions are such that the standard spheres have positive Ricci scalar.
The Lovelock invariants were taking from Fulling, S. A., et al. 'Normal forms for tensor polynomials. I. The Riemann tensor.' Classical and Quantum Gravity 9.5 (1992): 1151.

Type ?GRlist for a list of functions or ?'Function' for a function's description."

ExteriorAlgebra::usage = "Based on the EDC360 package. The exterior algebra part of this package defines the exterior derivative as d[]. The Wedge product can be used between forms d[x]\[Wedge]d[y].

The full list of functions:
Form2Metric
Metric
HodgeStar
SelfDualQ
AntiSelfDualQ
Form2Tensor
Tensor2Form
Lie
d
KillingQ
InteriorProduct"

Form2Metric::usage = "Form2Metric[g,xx] Converts a metric (abusively) written using squared one-forms,
eg. g=d[r]^2+r^2d[\[Theta]]^2 and xx={r,\[Theta]},
into its matrix representation."
Form2Metric::missingInput = "Cannot locate a globally defined array of coordinates, xx. Either provide one as an input or define it globally."

Metric::usage = "Metric[str,dim] Returns a metric, in tensor form, with shape dim x dim.

str='Sphere' or 'S' gives the standard metric on the unit sphere
str='Fubini-Study' or 'FS' gives the Fubini-Study metric on \!\(\*SuperscriptBox[\(CP\),\(dim\)]\)"

HodgeStar::usage = "HodgeStar[tensor, metric]: returns the Hodge dual of the tensor, with respect to the metric.
HodgeStar[form, metric, xx]: converts the form to its tensor representation then calls HodgeStar[tensor, metric]

If metric isn't specified, HodgeStar will try to use a globally defined variable called 'g' as the metric

Optional Arguments:
 - Assumptions->List. Used in all Simplify[] calls."
(*Star::usage = "Star[form, g, xx] is a custom operator for HodgeStar.";*)
HodgeStar::missingInput = "Missing required coordinates array xx or metric g when converting from form to tensor.
Common usage HodgeStar[form, g, xx]."

SelfDualQ::usage = "SelfDualQ[form, metric, coords]: return true if form is self-dual wrt the metric and false otherwise.
Accepts rank-r tensor and r-form representation for the input form. The former doesn't require coords to be specified.

SelfDualQ[form] will try to use globally defined variables called 'xx' as the set of coordinates and 'g' as the metric and then calls SelfDualQ[form, g, xx].

Optional Arguments:
  - Assumptions->None: Array of assumptions used in Simplify"
SelfDualQ::missingInput = "Cannot locate a globally defined array of coordinates, xx and/or metric g. Either provide them as an input or define them globally."

AntiSelfDualQ::usage = "AntiSelfDualQ[form, g, xx]: return true if form is anti-self-dual wrt the metric 'g' and false otherwise.
Accepts rank-r tensor and r-form representation for the input form. The former doesn't require coords to be specified.

AntiSelfDualQ[form] will try to use globally defined variables called 'xx' as the set of coordinates and 'g' as the metric and then calls AntiSelfDualQ[form, g, xx].

Optional Arguments:
  - Assumptions->None: Array of assumptions used in Simplify"
AntiSelfDualQ::missingInput = "Cannot locate a globally defined array of coordinates, xx and/or metric g. Either provide them as an input or define them globally."

Form2Tensor::usage = "Form2Tensor[form, xx]: Converts a d-form, represented using d[x]\[Wedge]d[y]... into its tensor components. 'xx' is a list of coordinates.
Form2Tensor[form] will try to use a globally defined variable called 'xx' as the set of coordinates and then calls Form2Tensor[form, xx].

Optional Arguments:
  - Assumptions->None: Array of assumptions used in Simplify
  - \"Type\"->String. When equal to 'SymmetrizedArray' returns the tensor as a fully anti-symmetric SymmetrizedArray. When set to 'Normal' returns a standard array. (defaults to 'Normal')"
Form2Tensor::missingInput = "Cannot locate a globally defined array of coordinates, xx. Either provide it as an input or define one globally."

Tensor2Form::usage = "Tensor2Form[tensor, xx]: Converts a tensor of any dimension into its form representation, using the d[] function. xx is a list of coordinates.
Tensor2Form[tensor] will try to use a globally defined variable called 'xx' as the set of coordinates and then calls Tensor2Form[tensor, xx].
If the original tensor has a symmetric part, it is automatically set to zero.

Optional Arguments:
  - Assumptions->None: Array of assumptions used in Simplify"
Tensor2Form::inputForm = "The input variable cannot be a differential form. Please provide a suitable tensor object, as an Array."
Tensor2Form::missingInput = "Cannot locate a globally defined array of coordinates, xx. Either provide it as an input or define one globally."

(* ScalarEDCode functions *)
d::usage = "d[form] applies the exterior derivative to the differential form 'form'.";
Wedge::usage = "Wedge[a, b] returns the wedge product of two differential forms 'a' and 'b'.";

Lie::usage = "Lie[vec, tensor, xx] Returns the Lie derivative of 'tensor' wrt to the vector 'vec' using the set of coordinates 'xx'.
A unit vector in the coordinates can be used as vector input by specifying an integer (in the range 1 to Length[xx]) instead of 'vec'.
A form of degree r can be used in lieu of 'tensor', after which it is converted to a tensor of rank r.

Lie[vec, tensor] will try to use a globally defined variable called 'xx' as the set of coordinates and then calls Lie[vec, tensor, xx].

Optional Arguments:
  - Assumptions->None: Array of assumptions used in Simplify
  - \"Up\"->Range[rank]: Array containing the set of indices in 'tensor' which are contravariant (or up)
  - \"Down\"->{}: Array containing the set of indices in 'tensor' which are covariant (or down)
  - Type->\"Other\": Used to specify the output format. Type->\"Form\" will output the tensor as a rank-r form (defaults to that when using a form as input)
"
Lie::vectorDim = "Dimensions of vector don't match that of the coordinates. Please provide a correct vector."
Lie::missingInput = "Cannot locate a globally defined array of coordinates, xx. Either provide it as an input or define one globally."

KillingQ::usage = "KillingQ[vector, g, xx]: returns true if the vector is Killing wrt the metric 'g' given the set of coordinates 'xx'.
KillingQ[vector, g] will try to use a globally defined variable called 'xx' as the set of coordinates and then calls KillingQ[vector, g, xx].
KillingQ[vector] will try to use a globally defined variable called 'g' as the metric and then calls KillingQ[vector, g].

Optional Arguments:
  - Assumptions->None: Array of assumptions used in Simplify"
KillingQ::missingInput = "Cannot locate a globally defined metric g or coordinates xx. Either provide them as input or define them globally."

Cov::usage = "Cov[vec, tensor, g, xx] Returns the covariant derivative of 'tensor' wrt to the vector 'vec' using the metric 'g' and set of coordinates 'xx'.
A unit vector in the coordinates can be used as vector input by specifying an integer (in the range 1 to Length[xx]) instead of 'vec'.
A form of degree r can be used in lieu of 'tensor', after which it is converted to a tensor of rank r.

Cov[vec, tensor, g] will try to use a globally defined variable called 'xx' as the set of coordinates and then calls Cov[vec, tensor, g, xx].
Cov[vec, tensor] will try to use a globally defined variable called 'g' as the metric and then calls Cov[vec, tensor, xx].

Optional Arguments:
  - Assumptions->None: Array of assumptions used in Simplify
  - \"Up\"->Range[rank]: Array containing the set of indices in 'tensor' which are contravariant (or up)
  - \"Down\"->{}: Array containing the set of indices in 'tensor' which are covariant (or down)
  - Type->\"Other\": Used to specify the output format. Type->\"Form\" will output the tensor as a rank-r form (defaults to that when using a form as input)
"
Cov::vectorDim = "Dimensions of vector don't match that of the coordinates. Please provide a correct vector."
Cov::missingInput = "Cannot locate a globally defined matrix g or array of coordinates, xx. Either provide them as a inputs or define them globally."

ConfKillingQ::usage = "ConfKillingQ[vector, g, xx]: returns true if the vector is conformal Killing wrt the metric 'g' given the set of coordinates 'xx'.
ConfKillingQ[vector, g] will try to use a globally defined variable called 'xx' as the set of coordinates and then calls ConfKillingQ[vector, g, xx].
ConfKillingQ[vector] will try to use a globally defined variable called 'g' as the metric and then calls ConfKillingQ[vector, g].

Optional Arguments:
  - Assumptions->None: Array of assumptions used in Simplify"
ConfKillingQ::missingInput = "Cannot locate a globally defined metric g or coordinates xx. Either provide them as a input or define them globally."

InteriorProduct::usage = "InteriorProduct[vec, form, xx] returns the interior product of the form with respect to the vector 'vec' given the set of coordinates 'xx'.
InteriorProduct[vec, tensor, xx] returns the interior product of the tensor with respect to the vector 'vec' given the set of coordinates 'xx'.

Optional Arguments:
  - Assumptions->None: Array of assumptions used in Simplify
  - Type->\"Other\": Used to specify the ouput format. Type->\"Form\" will output a rank r-1 form (defaults to that when using rank-r form as input)"
InteriorProduct::dimMismatch = "Mismatch between the length of the input vector and the length of the coordinates."

(*Needs["Notation`"];*)

Begin["`Private`"]

(* ---------- ScalarEDCode ---------- *)
Off[General::"spell", General::"spell1", Remove::"rmnsm", 
  UpSet::"write"];

Unprotect["Global`*"]; ClearAll["Global`*"]; Remove["Global`*"]; \
Unprotect[In, Out]; Clear[In, Out]; Protect[In, Out]; $Line = 0; \
$RecursionLimit = 256; $IterationLimit = 4096;

Forms[i_] := {}; AllScalForms = {}; AllDifForms = {}; AllSymbols = \
{}; FormVars = {_Wedge, _d}; HeadList = {}; ScalHeadList = {}; \
nodHeads = {Bar, Pattern, Condition, RuleDelayed, 
  SeriesData}; $EDCversion = 360;

zeroQ[0] = True; zeroQ[x_List] := And @@ (zeroQ /@ Union[Flatten[x]]);
 zeroQ[x_SeriesData] := If[x[[3]] === {}, True, False]; 
zeroQ[x_] := False;

SetAttributes[Bar, {Listable}]; Bar[Bar[x_]] = x; 
Bar[Complex[u_, v_]] := Complex[u, -v]; 
Bar[x_Plus | x_Times | x_Wedge | x_Power | x_Rule | x_Equal] := 
 Bar /@ x; 
Bar[x_SeriesData] := 
 x /. {First[x] -> Bar[First[x]], x[[2]] -> Bar[x[[2]]], 
   x[[3]] -> Bar /@ x[[3]]}; 
Bar[DirectedInfinity[x_]] := DirectedInfinity[Bar[x]]; 
Bar[d[x_]] := d[Bar[x]]; 
Bar[Derivative[x__][y_][z__]] := 
 If[Union[Bar[{z}]] === Union[{z}], 
  Derivative[
     Sequence @@ 
      Permutations[{x}][[
       Position[Permutations[{z}], Bar[{z}]][[1, 1]]]]][Bar[y]][z], 
  Derivative[x][Bar[y]][Sequence @@ Bar[{z}]]]; 
Bar[x_] := x /; NumericQ[x];

FormDegree[x_Plus] := FormDegree[First[x]]; 
FormDegree[x_Times] := Plus @@ FormDegree /@ List @@ x;
FormDegree[x_Wedge] := Plus @@ FormDegree /@ List @@ x;
FormDegree[d[x_]] := FormDegree[d[x]] = 1 + FormDegree[x];
FormDegree[x_List] := 
  FormDegree[First[Select[Union[Flatten[x]], ! zeroQ[#] &]]];
FormDegree[Bar[x_]] := FormDegree[Bar[x]] = FormDegree[x]; 
FormDegree[x_SeriesData] := If[x[[3]] === {}, 0, FormDegree[x[[3]]]];
FormDegree[x_] := 0;

DeclareForms[z__] := 
  Block[{h, x = {{z}}, xi, rxi, k, oldHeads, newHeads}, 
   While[Head[x[[1, 1]]] === List, x = First[x]]; 
   Do[xi = x[[i]]; rxi = Rest[xi]; 
    Forms[First[xi]] = Union[Forms[First[xi]], rxi];
    AllScalForms = Union[AllScalForms, rxi];
    Do[h = Head[rxi[[j]]]; 
     If[h === Symbol, FormDegree[rxi[[j]]] = First[xi]; 
      BasicScalFormQ[rxi[[j]]] = True;, FormDegree[_h] = First[xi]; 
      BasicScalFormQ[_h] = True;], {j, Length[rxi]}], {i, 
     Length[x]}];
   AllDifForms = Union[AllScalForms]; oldHeads = ScalHeadList; 
   ScalHeadList = 
    Complement[Union[Head /@ AllScalForms, ScalHeadList], {Symbol}]; 
   newHeads = Complement[ScalHeadList, oldHeads]; 
   AllSymbols = Union[AllDifForms];
   HeadList = Union[Head /@ AllSymbols]; 
   DifFormSymbols = Drop[AllDifForms, -Length[ScalHeadList]];
   k = Thread[Blank[ScalHeadList]]; 
   FormVars = 
    Flatten[{_Wedge, _d, Union[k, Bar[k]], 
      u_ | Bar[u_] /; MemberQ[DifFormSymbols, u], _tr}];];

NoDif[z__] := (nodHeads = Union[nodHeads, Flatten[{z}]];)

BasicScalFormQ[Bar[x_]] := BasicScalFormQ[x]; 
BasicScalFormQ[_d] = True; BasicScalFormQ[x_] := False;

ScalFormQ[x_Times] := Or @@ ScalFormQ /@ List @@ x; 
ScalFormQ[x_Wedge | x_Plus] := And @@ ScalFormQ /@ List @@ x; 
ScalFormQ[x_] := BasicScalFormQ[x];

SetAttributes[d, {Listable}]; 
d[x_Times | x_Wedge] := 
 d[First[x]]\[Wedge]Rest[x] + (-1)^FormDegree[First[x]]*
   First[x]\[Wedge]d[Rest[x]]; d[x_?NumericQ | x_d] = 0; 
d[Power[y_, n_]] := n y^(n - 1) d[y] + y^n Log[y] d[n]; 
d[x_Plus] := d /@ x; 
HoldPattern[d[Bar[x_]]] := 
 With[{evalHx = d[x]}, Bar[evalHx] /; Head[evalHx] =!= d]; 
d[x_Rule | x_Equal] := reWrite[d /@ x]; 
d[x_SeriesData] := (x /. x[[3]] -> d[x[[3]]]) + 
  Wedge[d[First[x]], D[x, First[x]]]; 
d[h_[y__] /; ! MemberQ[nodHeads, h]] := 
 Sum[(Derivative[
        Sequence @@ 
         RotateRight[Append[Table[0, {Length[{y}] - 1}], 1], i]][h][
      y]) d[{y}[[i]]], {i, Length[{y}]}] /; 
  FormDegree[h[y]] === 0 && ! 
    MemberQ[{Integer, Blank, Pattern, Condition}, Head[First[{y}]]];

newSer$[x_SeriesData, k_] := 
 SeriesData[First[x], x[[2]], 
  Flatten[Transpose[
    Prepend[Table[Table[0, {Length[x[[3]]]}], {k - 1}], x[[3]]]]], 
  k x[[4]], k x[[5]], k Last[x]]

Wedge[x_] := x /; Length[{x}] < 2 && Head[{x}[[1]]] =!= Pattern

Default[Wedge] := 1; SetAttributes[Wedge, {Flat, OneIdentity}]; 
Wedge[0, y__] = 0; Wedge[x__, 0] = 0; 
Wedge[x_SeriesData, y_SeriesData] := 
 Block[{x$, y$, r1, r2, res, x3, y3}, 
   If[Last[x] === Last[y], x$ = x; y$ = y, 
    x$ = newSer$[x, LCM[Last[x], Last[y]]/Last[x]]; 
    y$ = newSer$[y, LCM[Last[x], Last[y]]/Last[y]]]; 
   r1 = x$[[-3]] + y$[[-3]]; 
   r2 = Min[x$[[-2]] + y$[[-3]], x$[[-3]] + y$[[-2]]]; 
   If[Length[x$[[3]]] < x$[[-2]] - x$[[-3]], 
    x3 = Join[x$[[3]], 
      Table[0, {x$[[-2]] - x$[[-3]] - Length[x$[[3]]]}]], 
    x3 = x$[[3]]]; 
   If[Length[y$[[3]]] < y$[[-2]] - y$[[-3]], 
    y3 = Join[y$[[3]], 
      Table[0, {y$[[-2]] - y$[[-3]] - Length[y$[[3]]]}]], 
    y3 = y$[[3]]]; 
   Which[r2 === r1, res = {}, r2 - r1 === x$[[-2]] - x$[[-3]], 
    res = Sum[
      Map[Wedge[x3[[i]], #] &, 
       Join[Table[0, {i - 1}], Drop[y3, -i + 1]]], {i, Length[x3]}], 
    True, res = 
     Sum[Map[Wedge[#, y3[[i]]] &, 
       Join[Table[0, {i - 1}], Drop[x3, -i + 1]]], {i, Length[y3]}]]; 
   SeriesData[First[x$], x$[[2]], res, r1, r2, Last[x$]]] /; 
  First[x] === First[y] && x[[2]] === y[[2]]; 
Wedge[y_, x_SeriesData] := x /. x[[3]] -> Map[Wedge[y, #] &, x[[3]]]; 
Wedge[x_SeriesData, y_] := x /. x[[3]] -> Map[Wedge[#, y] &, x[[3]]];
Wedge[x__, y_Plus] := Plus @@ Map[Wedge[x, #] &, List @@ y];
Wedge[x_Plus, y__] := Plus @@ Map[Wedge[#, y] &, List @@ x]; 
Wedge[u__, Times[x_, y_]] := 
 Times[x, Wedge[u, y]] /; NumericQ[x] || FormDegree[x] === 0;
Wedge[Times[x_, y_], z__] := 
 Times[x, Wedge[y, z]] /; NumericQ[x] || FormDegree[x] === 0; 
x_^n_.\[Wedge]y_ := x^n*y /; FormDegree[x] === 0;
y_\[Wedge]x_^n_. := x^n*y /; FormDegree[x] === 0;
Wedge[x_, y___, x_] := 0 /; OddQ[FormDegree[x]] && ScalFormQ[x]; 
Wedge[x__] := 
 Block[{doubL = Transpose[{FormDegree /@ {x}, {x}}]}, 
   Signature[Select[doubL, OddQ[First[#]] &]] Wedge @@ 
     Map[Last[#1] &, Sort[doubL]]] /; 
  Union[BasicScalFormQ /@ {x}] === {True} && 
   Map[Last[#1] &, Sort[Transpose[{FormDegree /@ {x}, {x}}]]] =!= {x};

simpRules = {Cos[z_]^2*(x_.) + (x_.)*Sin[z_]^2 -> x}; varList = {};

coll[x_] := Collect[x, {_Wedge, _tr}, Factor];

SetAttributes[reWrite, {Listable}]; reWrite[0] = 0; 
reWrite[x_Equal] := Equal[reWrite[First[x] - Last[x]], 0]; 
reWrite[x_Rule] := Rule[First[x], reWrite[Last[x]]]; 
reWrite[x_SeriesData] := (x /. x[[3]] -> reWrite[x[[3]]]); 
reWrite[x_] := (Collect[coll[x] /. simpRules, FormVars, 
     bestFacRul] /. simpRules) /; FormDegree[x] > 0; 
reWrite[x_] := bestFacRul[x /. simpRules] /. simpRules;

FreeALLQ[x_, y_List] := And @@ Map[FreeQ[x, #] &, y];

bestFacRul[x_] := 
 If[varList === {} || FreeALLQ[x, varList], minFacRul[x], 
  varFacRul[x]]; 
varFacRul[x_] := 
 Collect[Expand[x] /. simpRules, varList, Factor] /. simpRules; 
minFacRul[x_] := Factor[Expand[x] /. simpRules] /. simpRules;

Unprotect[SeriesData]; If[$VersionNumber < 5, 
 Times[x_SeriesData, 0] ^= 0]; 
SeriesData[x1_, x2_, x3_, x4_, x5_, x6_] := 
 Plus @@ If[x2 === Infinity, 
    Map[First[#]/x1^((x4 + Last[#] - 1)/x6) &, 
      Transpose[{x3, Range[Length[x3]]}]] + 
     SeriesData[x1, x2, {}, x4, x5, x6],   
    Map[First[#]*(x1 - x2)^((x4 + Last[#] - 1)/x6) &, 
      Transpose[{x3, Range[Length[x3]]}]] + 
     SeriesData[x1, x2, {}, x4, x5, x6]] /; ! 
   FreeQ[x3, x1]; Protect[SeriesData];

On[General::"spell", General::"spell1", UpSet::"write"];


(* ---------- Curvature Invariants ---------- *)

Christoffel[g_, xx_] := 
  Table[1/2 Sum[
      Inverse[g][[\[Nu], \[Sigma]]] (D[g[[\[Rho], \[Sigma]]], 
          xx[[\[Mu]]]] + D[g[[\[Mu], \[Sigma]]], xx[[\[Rho]]]] - 
         D[g[[\[Mu], \[Rho]]], xx[[\[Sigma]]]]), {\[Sigma], 1, 
       Length[xx]}], {\[Nu], 1, Length[xx]}, {\[Mu], 1, 
     Length[xx]}, {\[Rho], 1, Length[xx]}] // Simplify;

RiemannTensor[g_, xx_] := 
  TensorContract[Inactive[TensorProduct][g,Transpose[Grad[Christoffel[g,xx],xx],{1,4,2,3}]-Transpose[Grad[Christoffel[g,xx],xx],{1,3,2,4}]+Transpose[TensorContract[TensorProduct[Christoffel[g,xx],Christoffel[g,xx]],{{3,4}}],{1,3,4,2}]-Transpose[TensorContract[TensorProduct[Christoffel[g,xx],Christoffel[g,xx]],{{3,4}}],{1,4,3,2}]],{{2,3}}] // Activate // Simplify;

RicciTensor[g_, xx_] := 
  TensorContract[
    Inactive[TensorProduct][Inverse[g], 
     RiemannTensor[g, xx]], {{1, 3}, {2, 5}}] // Activate // Simplify;

RicciScalar[g_, xx_] := 
  TensorContract[
    Inactive[TensorProduct][Inverse[g], RicciTensor[g, xx]], {{1, 3}, {2, 4}}] // Activate //
    Simplify;

KretschmannScalar[g_, xx_] := 
  TensorContract[Inactive[TensorProduct][Inverse[g],Inverse[g],Inverse[g],Inverse[g],RiemannTensor[g,xx],RiemannTensor[g,xx]],{{1,9},{2,13},{3,10},{4,14},{5,11},{6,15},{7,12},{8,16}}]//Activate // Simplify;

WeylTensor[g_, xx_] := 
    If[Length[xx]<4,Table[0,{\[Mu],1,Length[xx]},{\[Nu],1,Length[xx]},{\[Rho],1,Length[xx]},{\[Sigma],1,Length[xx]}],
  RiemannTensor[g, xx] + 
    1/(Length[xx] - 
      2) (TensorTranspose[
        TensorProduct[RicciTensor[g, xx], g], {1, 4, 2, 3}] - 
       TensorTranspose[
        TensorProduct[RicciTensor[g, xx], g], {1, 3, 2, 4}] + 
       TensorTranspose[
        TensorProduct[RicciTensor[g, xx], g], {2, 3, 1, 4}] - 
       TensorTranspose[
        TensorProduct[RicciTensor[g, xx], g], {2, 4, 1, 3}]) + 
    1/((Length[xx] - 1) (Length[xx] - 2))
      RicciScalar[g, 
      xx] (TensorTranspose[TensorProduct[g, g], {1, 3, 2, 4}] - 
       TensorTranspose[TensorProduct[g, g], {1, 4, 2, 3}]) // Simplify];

CottonTensor[g_,xx_] :=
    If[Length[xx]<3,
        Table[0,{\[Mu],1,Length[xx]},{\[Nu],1,Length[xx]},{\[Rho],1,Length[xx]}],
        Grad[RicciTensor[g,xx],xx]-Transpose[Activate[TensorContract[Inactive[TensorProduct][Christoffel[g,xx],RicciTensor[g,xx]],{{1,4}}]],{3,1,2}]-Transpose[Activate[TensorContract[Inactive[TensorProduct][Christoffel[g,xx],RicciTensor[g,xx]],{{1,5}}]],{3,2,1}]-Transpose[Grad[RicciTensor[g,xx],xx],{1,3,2}]+Transpose[Activate[TensorContract[Inactive[TensorProduct][Christoffel[g,xx],RicciTensor[g,xx]],{{1,4}}]],{2,1,3}]+Transpose[Activate[TensorContract[Inactive[TensorProduct][Christoffel[g,xx],RicciTensor[g,xx]],{{1,5}}]],{2,3,1}] // Simplify
    ];

SchoutenTensor[g_,xx_] :=
    If[Length[xx]<3,Table[0,{\[Mu],1,Length[xx]},{\[Nu],1,Length[xx]},{\[Rho],1,Length[xx]},{\[Sigma],1,Length[xx]}],
    1/(Length[xx]-2) (RicciTensor[g,xx]-RicciScalar[g,xx]/(2(Length[xx]-1)) g)
    ]

EulerScalar[g_, 
   xx_] := -Lovelock[2,4][g, xx] + 2 Lovelock[2,2][g, xx] - 
    2/3 RicciScalar[g, xx]^2 // Simplify;

(*Lovelock Invariants*)
Lovelock[0,1][g_,xx_] := 1;

Lovelock[1,1][g_, xx_] := RicciScalar[g, xx];

Lovelock[2,1][g_, xx_] := RicciScalar[g, xx]^2;

Lovelock[2,2][g_, xx_] := 
  TensorContract[
     Inactive[TensorProduct][Inverse[g], Inverse[g], RicciTensor[g, xx], 
      RicciTensor[g, xx]], {{1, 5}, {2, 7}, {3, 6}, {4, 8}}] - 
    2/3 RicciScalar[g, xx]^2 // Activate // Simplify;

Lovelock[2,3][g_, xx_] := KretschmannScalar[g, xx];

Lovelock[2,4][g_, xx_] := 
    TensorContract[Inactive[TensorProduct][Inverse[g],Inverse[g],Inverse[g],Inverse[g],WeylTensor[g,xx],WeylTensor[g,xx]],{{1,9},{2,13},{3,10},{4,14},{5,11},{6,15},{7,12},{8,16}}]//Activate // Simplify;

Lovelock[3,1][g_, xx_] := RicciScalar[g, xx]^3;

Lovelock[3,2][g_, xx_] := RicciScalar[g, xx] Lovelock[2,2][g, xx];

Lovelock[3,3][g_, xx_] := 
  TensorContract[
    Inactive[TensorProduct][Inverse[g], Inverse[g], Inverse[g], 
     RicciTensor[g, xx], RicciTensor[g, xx], 
     RicciTensor[g, xx]], {{1, 7}, {2, 9}, {3, 8}, {4, 11}, {5, 
      10}, {6, 12}}] // Activate // Simplify;

Lovelock[3,4][g_, xx_] := 
  TensorContract[
    Inactive[TensorProduct][Inverse[g], Inverse[g], Inverse[g], Inverse[g], 
     RicciTensor[g, xx], RicciTensor[g, xx], 
     RiemannTensor[g, xx]], {{1, 9}, {2, 13}, {3, 10}, {4, 15}, {5, 
      11}, {6, 14}, {7, 12}, {8, 16}}] // Activate // Simplify;

Lovelock[3,5][g_, xx_] := 
  RicciScalar[g, xx] KretschmannScalar[g, xx] // Simplify;

Lovelock[3,6][g_, xx_] := 
  TensorContract[
    Inactive[TensorProduct][Inverse[g], Inverse[g], Inverse[g], Inverse[g], 
     Inverse[g], RicciTensor[g, xx], RiemannTensor[g, xx], 
     RiemannTensor[g, xx]], {{1, 11}, {2, 16}, {3, 12}, {4, 20}, {5, 
      13}, {6, 17}, {7, 14}, {8, 18}, {9, 15}, {10, 19}}] // Activate // 
   Simplify;

Lovelock[3,7][g_, xx_] := 
  TensorContract[
    Inactive[TensorProduct][Inverse[g], Inverse[g], Inverse[g], Inverse[g], 
     Inverse[g], Inverse[g], RiemannTensor[g, xx], 
     RiemannTensor[g, xx], 
     RiemannTensor[g, xx]], {{1, 13}, {2, 17}, {3, 14}, {4, 18}, {5, 
      15}, {6, 21}, {7, 16}, {8, 22}, {9, 19}, {10, 23}, {11, 
      20}, {12, 24}}] // Activate // Simplify;

Lovelock[3,8][g_, xx_] := 
  TensorContract[
    Inactive[TensorProduct][Inverse[g], Inverse[g], Inverse[g], Inverse[g], 
     Inverse[g], Inverse[g], RiemannTensor[g, xx], 
     RiemannTensor[g, xx], 
     RiemannTensor[g, xx]], {{1, 13}, {2, 17}, {3, 14}, {4, 21}, {5, 
      15}, {6, 19}, {7, 16}, {8, 23}, {9, 18}, {10, 22}, {11, 
      20}, {12, 24}}] // Activate // Simplify;

Lovelock[3,9][g_,xx_] :=
    TensorContract[
  Inactive[TensorProduct][Inverse[g], Inverse[g], Inverse[g], 
   Inverse[g], WeylTensor[g, xx], RicciTensor[g, xx], 
   RicciTensor[g, xx]], {{1, 9}, {2, 13}, {3, 10}, {4, 15}, {5, 
    11}, {6, 14}, {7, 12}, {8, 16}}] // Activate // Simplify;

Lovelock[3,10][g_,xx_] :=
    TensorContract[
   Inactive[TensorProduct][Inverse[g], Inverse[g], Inverse[g], 
    Inverse[g], Inverse[g], WeylTensor[g, xx], WeylTensor[g, xx], 
    RicciTensor[g, xx]], {{1, 11}, {2, 19}, {3, 12}, {4, 18}, {5, 
     13}, {6, 15}, {7, 14}, {8, 16}, {9, 17}, {10, 20}}] // Activate // Simplify;

Lovelock[3,11][g_,xx_] :=
    TensorContract[
   Inactive[TensorProduct][Inverse[g], Inverse[g], Inverse[g], 
    Inverse[g], Inverse[g], Inverse[g], WeylTensor[g, xx], 
    WeylTensor[g, xx], 
    WeylTensor[g, xx]], {{1, 13}, {2, 23}, {3, 14}, {4, 24}, {5, 
     15}, {6, 17}, {7, 16}, {8, 18}, {9, 19}, {10, 21}, {11, 20}, {12,
      22}}] // Activate // Simplify;

Options[LovelockAll] = {"Extend" -> True};
LovelockAll[g_,xx_,OptionsPattern[]] := 
    If[OptionValue["Extend"] == True,
    {
        Lovelock[1,1][g,xx],
        Lovelock[2,1][g,xx],
        Lovelock[2,2][g,xx],
        Lovelock[2,3][g,xx],
        Lovelock[2,4][g,xx],
        Lovelock[3,1][g,xx],
        Lovelock[3,2][g,xx],
        Lovelock[3,3][g,xx],
        Lovelock[3,4][g,xx],
        Lovelock[3,5][g,xx],
        Lovelock[3,6][g,xx],
        Lovelock[3,7][g,xx],
        Lovelock[3,8][g,xx],
        Lovelock[3,9][g,xx],
        Lovelock[3,10][g,xx],
        Lovelock[3,11][g,xx]
    },
    {
        Lovelock[1,1][g,xx],
        Lovelock[2,1][g,xx],
        Lovelock[2,2][g,xx],
        Lovelock[2,3][g,xx],
        Lovelock[3,1][g,xx],
        Lovelock[3,2][g,xx],
        Lovelock[3,3][g,xx],
        Lovelock[3,4][g,xx],
        Lovelock[3,5][g,xx],
        Lovelock[3,6][g,xx],
        Lovelock[3,7][g,xx],
        Lovelock[3,8][g,xx]
    }]

Options[GRInvariants] = {"Minimal"->False};
GRInvariants[g_,xx_,OptionsPattern[]] := 
    If[OptionValue["Minimal"] == True,
    {
        Lovelock[1,1][g,xx],
        Lovelock[2,2][g,xx],
        Lovelock[2,3][g,xx],
        Lovelock[2,4][g,xx],
        Lovelock[3,3][g,xx],
        Lovelock[3,4][g,xx],
        Lovelock[3,6][g,xx],
        Lovelock[3,7][g,xx],
        Lovelock[3,8][g,xx],
        Lovelock[3,9][g,xx],
        Lovelock[3,10][g,xx],
        Lovelock[3,11][g,xx]
    },
    {
        RicciScalar[g, xx],
        KretschmannScalar[g, xx],
        EulerScalar[g, xx],
        Lovelock[1,1][g,xx],
        Lovelock[2,1][g,xx],
        Lovelock[2,2][g,xx],
        Lovelock[2,3][g,xx],
        Lovelock[2,4][g,xx],
        Lovelock[3,1][g,xx],
        Lovelock[3,2][g,xx],
        Lovelock[3,3][g,xx],
        Lovelock[3,4][g,xx],
        Lovelock[3,5][g,xx],
        Lovelock[3,6][g,xx],
        Lovelock[3,7][g,xx],
        Lovelock[3,8][g,xx],
        Lovelock[3,9][g,xx],
        Lovelock[3,10][g,xx],
        Lovelock[3,11][g,xx]
    }]


(* ---------- Diff Geo ---------- *)
Options[Form2Metric] = {Assumptions->None};
Form2Metric[g_, xx_, OptionsPattern[]] := Simplify[Table[If[i === j, 1, 1/2] Coefficient[g, d[xx[[i]]] d[ xx[[j]]]], {i, 
    Length[xx]}, {j, Length[xx]}], OptionValue[Assumptions]];
Form2Metric[g_, OptionsPattern[]] := If[If[ValueQ[Symbol["Global`" <> SymbolName[xx]]], TensorRank[Symbol["Global`" <> SymbolName[xx]]] == 1, False],
  Form2Metric[g, Symbol["Global`" <> SymbolName[xx]]],

  Message[Form2Metric::missingInput];Return[$Failed]
]

Metric[str_, dim_] := Module[{gform, xx},
  If[str == "S" || str == "s" || str == "Sphere" || str == "sphere" || str == "Spherical" || str == "spherical",
    Print["Standard spherical metric in ", ToString[dim], " dimensions"];
    xx = Global`\[Theta][#] & /@ Range[dim];
    Print[Join[StringForm["0 < `` < \[Pi] ", #] & /@ xx[[;; -2]], {StringForm["0 < `` < 2\[Pi] ", xx[[-1]]]}] // MatrixForm];
    gform = Total[d[#]^2 Times @@ (Sin[##]^2 & /@ xx[[1 ;; Position[xx, #][[1]][[1]] - 1]]) & /@ xx];
  ];
  If[str == "FS" || str == "fs" || str == "Fs" || str == "fS" || str == "FubiniStudy" || str == "fubinistudy" || str == "Fubinistudy" || str == "fubiniStudy" || str == "Fubini-Study" || str == "Fubini-study" || str == "fubini-Study" || str == "fubini-study",
    Print[StringForm["Fubini-Study metric on \!\(\*SuperscriptBox[\(CP\),\(``\)]\)",dim], " where z[i] \[Element] C"];
    xx = Flatten[{Global`z[#],Global`zb[#]} &/@ Range[dim]];
    gform = Sum[d[xx[[2*ii - 1]]] d[xx[[2*ii]]], {ii, 1, dim}]/(1 + Sum[xx[[2*ii - 1]] xx[[2*ii]], {ii, 1, dim}])-(Sum[xx[[2*ii]]d[xx[[2*ii-1]]],{ii,1,dim}]Sum[xx[[2*ii-1]]d[xx[[2*ii]]],{ii,1,dim}])/(1 + Sum[xx[[2*ii - 1]] xx[[2*ii]], {ii, 1, dim}])//Expand
  ];
  Return[{xx, Form2Metric[gform, xx]}]
]

Options[HodgeStar] = {Assumptions->None, Type->None};
HodgeStar[formT_, gT_, xx_, OptionsPattern[]] := 
 Module[{tensor, g, m, dim, ginv, dualRules = {}, dual, sqrtDetg, i, pos, val, ind, compl, sign, type = OptionValue[Type]},

  If[FormDegree[formT] != 0,
    tensor = Form2Tensor[formT,xx,Type->"SymmetrizedArray"];,
    tensor = SymmetrizedArray[formT,Automatic, Antisymmetric[All]];
  ];
  If[FormDegree[gT] != 0,
    g = SparseArray[Form2Metric[gT,xx]];,
    g = SparseArray[gT];
  ];

  If[type == None, type = "Form"];
  
  m = TensorRank[tensor];
  If[StringContainsQ[ToString[m], "TensorRank"], m=0]; (* TensorRank gets confused if you provide a rank 0 tensor with a variable *)
  dim = Dimensions[g][[1]];
  ginv = SymmetrizedArray[Simplify[Inverse[g],OptionValue[Assumptions]], Automatic, Symmetric[All]];
  sqrtDetg = Simplify[Sqrt[Abs[Det[g]]],OptionValue[Assumptions]];

  If[m == 0,
    dual = Simplify[tensor * sqrtDetg, OptionValue[Assumptions]] * Wedge @@ d[xx[[#]] & /@ Range[dim]];
    If[type == "Form", 
      Return[dual];,
      Return[Form2Tensor[dual,xx]];
    ];,

    (* Bring up all indices for the input form, using the inverse metric *)
    tensorUp = 
    Simplify[
    SymmetrizedArray[
      TensorContract[
        Inactive[TensorProduct][
        Sequence @@ Join[{tensor}, Table[ginv, m]]], 
        Table[{i, m + 2 i - 1}, {i, 1, m}]] // Activate, Automatic, 
      Antisymmetric[All]], OptionValue[Assumptions]];

    (* Add unique elements of tensorUp to the dual tensor *)
    tensorUpRules = Most[ArrayRules[tensorUp]];
    tensorUpIndices = DeleteDuplicates[Sort[#] & /@ tensorUpRules[[All, 1]]];
    For[i = 1, i <= Length[tensorUpIndices], i++,
      pos = Position[tensorUpRules[[All, 1]], tensorUpIndices[[i]]][[1, 1]];
      val = tensorUpRules[[pos, 2]];
      ind = tensorUpRules[[pos, 1]];
      compl = Complement[Range[dim], ind];
      sign = Signature[Join[ind, compl]];
      AppendTo[dualRules, compl -> sign val];
    ];
    dual = Simplify[sqrtDetg * SparseArray[dualRules, Table[dim, dim - m]], OptionValue[Assumptions]];
    dual = SymmetrizedArray[(dim - m)! * dual, Automatic, Antisymmetric[All]];
    If[type == "Form", 
      Return[Tensor2Form[dual, xx]];,
      Return[dual];
    ];
  ];
]
HodgeStar[formT_, gT_, OptionsPattern[]] := If[FormDegree[formT] != 0 || FormDegree[gT] != 0 || OptionValue[Type] == "Form",
  Message[HodgeStar::missingInput];Return[$Failed],

  HodgeStar[formT, gT, {}, OptionsPattern[]]];
HodgeStar[formT_, OptionsPattern[]] := If[If[ValueQ[Symbol["Global`" <> SymbolName[g]]], TensorRank[Symbol["Global`" <> SymbolName[g]]] == 2, False],
  HodgeStar[form, Symbol["Global`" <> SymbolName[g]], Symbol["Global`" <> SymbolName[xx]]],

  Message[HodgeStar::missingInput];Return[$Failed]
]

(* Define the custom Star function that will wrap HodgeStar
Star[form_, g_, xx_] := HodgeStar[form, g, xx];

(* Define the custom notation for the \[Star] operator *)
InfixNotation[ParsedBoxWrapper[
  SubsuperscriptBox["\[Star]", "xx_", "g_"]], 
  Star];*)

Options[SelfDualQ] = {Assumptions->None};
SelfDualQ[formT_, gT_, xx_, OptionsPattern[]] := Simplify[formT-HodgeStar[formT,gT,xx],OptionValue[Assumptions]]===0;
SelfDualQ[formT_, OptionsPattern[]] := If[If[ValueQ[Symbol["Global`" <> SymbolName[xx]]], TensorRank[Symbol["Global`" <> SymbolName[xx]]] == 1, False] && If[ValueQ[Symbol["Global`" <> SymbolName[g]]], TensorRank[Symbol["Global`" <> SymbolName[g]]] == 2, False],
  SelfDualQ[form, Symbol["Global`" <> SymbolName[g]], Symbol["Global`" <> SymbolName[xx]]],

  Message[SelfDualQ::missingInput];Return[$Failed]
]

Options[AntiSelfDualQ] = {Assumptions->None};
AntiSelfDualQ[formT_, gT_, xx_, OptionsPattern[]] := Simplify[formT+HodgeStar[formT,gT,xx],OptionValue[Assumptions]]===0;
AntiSelfDualQ[formT_, OptionsPattern[]] := If[If[ValueQ[Symbol["Global`" <> SymbolName[xx]]], TensorRank[Symbol["Global`" <> SymbolName[xx]]] == 1, False] && If[ValueQ[Symbol["Global`" <> SymbolName[g]]], TensorRank[Symbol["Global`" <> SymbolName[g]]] == 2, False],
  AntiSelfDualQ[form, Symbol["Global`" <> SymbolName[g]], Symbol["Global`" <> SymbolName[xx]]],

  Message[AntiSelfDualQ::missingInput];Return[$Failed]
]

Options[Form2Tensor] = {Type->"Normal", Assumptions->None};
Form2Tensor[form_, xx_, OptionsPattern[]] := Module[
 {n = FormDegree[form], dim = Length[xx], terms, tensor, wedgeTerm, 
  coefficient, permutedIndices, antisymtensor},
  terms = If[Head[Expand@form] === Plus, List @@ Expand@form, {Expand@form}];
  tensor = ConstantArray[0, Array[dim &, n]];

  If[n==0, Return[form],
  Do[
  If[n==1,
    wedgeTerm = Select[terms[[i]], MatchQ[#, _d] &];,
    wedgeTerm = Select[terms[[i]], MatchQ[#, _Wedge] &];];
  
  If[ToString[wedgeTerm] == "Wedge[]" || ToString[wedgeTerm] == "d[]" || ToString[wedgeTerm] == "1",
   wedgeTerm = terms[[i]];
   coefficient = 1;
   ,
   coefficient = Coefficient[terms[[i]], wedgeTerm];
   ];
  
  permutedIndices = 
  Position[xx, #] & /@ (List @@ wedgeTerm /. 
     HoldPattern[d[var_]] :> var);
  If[MemberQ[permutedIndices, {}], coefficient = 0;];
  permutedIndices = Flatten[permutedIndices];
  
  tensor[[Sequence @@ permutedIndices]] += coefficient;
  , {i, Length[terms]}
  ];
 antisymtensor = Simplify[n! SymmetrizedArray[tensor, Automatic, 
    Antisymmetric[All]],OptionValue[Assumptions]];
  If[OptionValue[Type] == "SymmetrizedArray", Return[antisymtensor];];
  If[OptionValue[Type] == "SparseArray", Return[SparseArray[antisymtensor]];,
    Return[Normal[antisymtensor]];];
  ]
];
Form2Tensor[form_, OptionsPattern[]] := If[If[ValueQ[Symbol["Global`" <> SymbolName[xx]]], TensorRank[Symbol["Global`" <> SymbolName[xx]]] == 1, False],
  Form2Tensor[form, Symbol["Global`" <> SymbolName[xx]]],

  Message[Form2Tensor::missingInput];Return[$Failed]
]

Options[Tensor2Form] = {Assumptions->None};
Tensor2Form[tensor_, xx_, OptionsPattern[]] := Module[
  {n = TensorRank[tensor], dim = Length[xx], sparse, form, rules, indices, pos, val, ind},
  If[FormDegree[tensor] != 0,
    Message[Tensor2Form::inputForm];
    Return[$Failed]];

  If[n == 0, Return[Simplify[tensor,OptionValue[Assumptions]]],
    If[n == 1,
      Return[Simplify[Sum[Normal[tensor][[ii]] d[xx[[ii]]], {ii, 1, dim}],OptionValue[Assumptions]]],
      sparse = SymmetrizedArray[tensor ,Automatic , Antisymmetric[All]];
      form=0;
      rules = Most[ArrayRules[sparse]];
      indices = DeleteDuplicates[Sort[#] & /@ rules[[All, 1]]];
      For[i = 1, i <= Length[indices], i++, 
        pos = Position[rules[[All, 1]], indices[[i]]][[1, 1]];
        val = rules[[pos, 2]];
        ind = rules[[pos, 1]];
        form += val*Wedge @@ d[xx[[#]] & /@ ind];
        ];
      Return[Simplify[form,OptionValue[Assumptions]]]
    ];
  ];
]
Tensor2Form[tensor_, OptionsPattern[]] := If[If[ValueQ[Symbol["Global`" <> SymbolName[xx]]], TensorRank[Symbol["Global`" <> SymbolName[xx]]] == 1, False],
  Tensor2Form[tensor, Symbol["Global`" <> SymbolName[xx]]],

  Message[Tensor2Form::missingInput];Return[$Failed]
]

Options[Lie] = {Assumptions->None, "Up"->{}, "Down"->{}, Type->"Other"};
Lie[vecT_, tensorT_, xx_, OptionsPattern[]] := Module[
  {r, dim, up, down, lie, tensor, vec, type},

  up = OptionValue["Up"];
  down = OptionValue["Down"];
  type = OptionValue[Type];

  If[FormDegree[tensorT] != 0,
    r = FormDegree[tensorT];
    up = {};
    down = Range[r];
    tensor = Form2Tensor[tensorT, xx];
    If[type == "Other", type = "Form";];,

    tensor = tensorT;
    r = TensorRank[tensor];
    If[StringContainsQ[ToString[r], "TensorRank"], r=0];
  ];

  If[Length[up]==0 && Length[down]==0,
    up=Range[r];
    down={};
  ];
  If[Length[up]==0 && Length[down]!=0,
    up=Complement[Range[r],down];
  ];
  If[Length[up]!=0 && Length[down]==0,
    down=Complement[Range[r],up];
  ];

  If[Length[vecT] == 0,
    (* check if vec in xx, if yes vec = {0,0,...,1,...,0} given by pos vec in xx. Otherwise error*)
    If[vecT > Length[xx] || vecT <= 0,
      Message[Lie::vectorDim],

      vec = Table[If[i == vecT, 1, 0], {i, 1, Length[xx]}];
    ];,

    If[Length[vecT] != Length[xx] || TensorRank[vecT] > 1,
      Message[Lie::vectorDim],
      
      vec = vecT];
  ];

  lie = TensorContract[TensorProduct[Grad[tensor, xx], vec], {r + 1, r + 2}];
  
  If[r > 0,
    lie += -Total[TensorTranspose[TensorContract[TensorProduct[Grad[vec, xx], tensor], {2, 2 + #}], Join[{#}, Complement[Range[r], {#}]]] & /@ up];
    lie += Total[TensorTranspose[TensorContract[TensorProduct[Grad[vec, xx], tensor], {1, 2 + #}], Join[{#}, Complement[Range[r], {#}]]] & /@ down];
  ];
  lie = Simplify[lie, OptionValue[Assumptions]];
  If[type == "Form",
    Return[Tensor2Form[lie,xx]];,
    Return[lie];
  ];
]
Lie[vecT_, tensorT_, OptionsPattern[]] := If[If[ValueQ[Symbol["Global`" <> SymbolName[xx]]], TensorRank[Symbol["Global`" <> SymbolName[xx]]] == 1, False],
  Lie[vecT, tensorT, Symbol["Global`" <> SymbolName[xx]], Assumptions->OptionValue[Assumptions], "Up"->OptionValue["Up"], "Down"->OptionValue["Down"], Type->OptionValue[Type]],

  Message[Lie::missingInput];Return[$Failed]
]

Options[KillingQ] = {Assumptions->None};
KillingQ[vec_, g_, xx_, OptionsPattern[]] := Return[Simplify[Lie[vec, g, xx, "Down"->{1,2}] === Table[0, {i, 1, Length[xx]}, {j, 1, Length[xx]}],OptionValue[Assumptions]]]
KillingQ[vec_, g_, OptionsPattern[]] := If[If[ValueQ[Symbol["Global`" <> SymbolName[xx]]], TensorRank[Symbol["Global`" <> SymbolName[xx]]] == 1, False],
  KillingQ[vec, g, Symbol["Global`" <> SymbolName[xx]]],

  Message[KillingQ::missingInput];Return[$Failed]
]
KillingQ[vec_, OptionsPattern[]] := If[If[ValueQ[Symbol["Global`" <> SymbolName[xx]]], TensorRank[Symbol["Global`" <> SymbolName[xx]]] == 1, False] && If[ValueQ[Symbol["Global`" <> SymbolName[g]]], TensorRank[Symbol["Global`" <> SymbolName[g]]] == 2, False],
  KillingQ[vec, Symbol["Global`" <> SymbolName[g]], Symbol["Global`" <> SymbolName[xx]]],

  Message[KillingQ::missingInput];Return[$Failed]
]

Options[Cov] = {Assumptions->None, "Up"->{}, "Down"->{}, Type->"Other"};
Cov[vecT_, tensorT_, g_, xx_, OptionsPattern[]] := Module[
  {r, dim, up, down, cov, tensor, vec, type, csymbols},

  up = OptionValue["Up"];
  down = OptionValue["Down"];
  type = OptionValue[Type];

  If[FormDegree[tensorT] != 0,
    r = FormDegree[tensorT];
    up = {};
    down = Range[r];
    tensor = Form2Tensor[tensorT, xx];
    If[type == "Other", type = "Form";];,

    tensor = tensorT;
    r = TensorRank[tensor];
    If[StringContainsQ[ToString[r], "TensorRank"], r=0];
  ];

  If[Length[up]==0 && Length[down]==0,
    up=Range[r];
    down={};
  ];
  If[Length[up]==0 && Length[down]!=0,
    up=Complement[Range[r],down];
  ];
  If[Length[up]!=0 && Length[down]==0,
    down=Complement[Range[r],up];
  ];

  If[Length[vecT] == 0,
    (* check if vec in xx, if yes vec = {0,0,...,1,...,0} given by pos vec in xx. Otherwise error*)
    If[vecT > Length[xx] || vecT <= 0,
      Message[Del::vectorDim],

      vec = Table[If[i == vecT, 1, 0], {i, 1, Length[xx]}];
    ];,

    If[Length[vecT] != Length[xx] || TensorRank[vecT] > 1,
      Message[Del::vectorDim],
      
      vec = vecT];
  ];

  cov = TensorContract[TensorProduct[Grad[tensor, xx], vec], {r + 1, r + 2}];
  csymbols = TensorContract[TensorProduct[Christoffel[g, xx], vec], {2, 4}];
  
  If[r > 0,
    cov += TensorTranspose[TensorContract[TensorProduct[csymbols, tensor], {2, 2 + #}], Join[{#}, Complement[Range[r], {#}]]] & /@ up;
    cov += -TensorTranspose[TensorContract[TensorProduct[csymbols, tensor], {1, 2 + #}], Join[{#}, Complement[Range[r], {#}]]] & /@ down;
  ];
  cov = Simplify[cov, OptionValue[Assumptions]];
  If[type == "Form",
    Return[Tensor2Form[cov,xx]];,
    Return[cov];
  ];
]
Cov[vecT_, tensorT_, g_, OptionsPattern[]] := If[If[ValueQ[Symbol["Global`" <> SymbolName[xx]]], TensorRank[Symbol["Global`" <> SymbolName[xx]]] == 1, False],
  Cov[vecT, tensorT, g, Symbol["Global`" <> SymbolName[xx]]],

  Message[Cov::missingInput];Return[$Failed]
]
Cov[vecT_, tensorT_, OptionsPattern[]] := If[If[ValueQ[Symbol["Global`" <> SymbolName[g]]], TensorRank[Symbol["Global`" <> SymbolName[g]]] == 2, False],
  Cov[vecT, tensorT, Symbol["Global`" <> SymbolName[g]], Symbol["Global`" <> SymbolName[xx]]],

  Message[Cov::missingInput];Return[$Failed]
]

Options[ConfKillingQ] = {Assumptions->None};
ConfKillingQ[vec_, g_, xx_, OptionsPattern[]] := Return[Simplify[Lie[vec, g, xx, "Down"->{1,2}] - 2/Length[xx] Sum[Cov[i, vec, g, xx][[i]], {i, 1, Length[xx]}] g === Table[0, {i, 1, Length[xx]}, {j, 1, Length[xx]}],OptionValue[Assumptions]]]
ConfKillingQ[vec_, g_, OptionsPattern[]] := If[If[ValueQ[Symbol["Global`" <> SymbolName[xx]]], TensorRank[Symbol["Global`" <> SymbolName[xx]]] == 1, False],
  ConfKillingQ[vec, g, Symbol["Global`" <> SymbolName[xx]]],

  Message[ConfKillingQ::missingInput];Return[$Failed]
]
ConfKillingQ[vec_, OptionsPattern[]] := If[If[ValueQ[Symbol["Global`" <> SymbolName[xx]]], TensorRank[Symbol["Global`" <> SymbolName[xx]]] == 1, False] && If[ValueQ[Symbol["Global`" <> SymbolName[g]]], TensorRank[Symbol["Global`" <> SymbolName[g]]] == 2, False],
  ConfKillingQ[vec, Symbol["Global`" <> SymbolName[g]], Symbol["Global`" <> SymbolName[xx]]],

  Message[ConfKillingQ::missingInput];Return[$Failed]
]

Options[InteriorProduct] = {Assumptions->None, Type->None};
InteriorProduct[vec_, formT_, xx_, OptionsPattern[]] := 
 Module[{tensor, m, dim, intRules = {}, inttensor, i, pos, val, ind, type = OptionValue[Type]},

  If[FormDegree[formT] != 0,
    tensor = Form2Tensor[formT,xx,Type->"SymmetrizedArray"];,

    tensor = SymmetrizedArray[formT, Automatic, Antisymmetric[All]];
  ];

  If[Length[vec] != Length[xx],
    Message[InteriorProduct::dimMismatch];Return[$Failed]
  ];

  If[type == None, type = "Form"];
  
  m = TensorRank[tensor];
  If[StringContainsQ[ToString[m], "TensorRank"], m=0]; (* TensorRank gets confused if you provide a rank 0 tensor with a variable *)
  dim = Length[xx];

  If[m == 0,
    inttensor = SymmetrizedArray[{}, Table[dim,r-1], Antisymmetric[All]];
    If[type == "Form", 
      Return[0];,
      Return[inttensor];
    ];,

    tensorRules = Most[ArrayRules[tensor]];
    tensorIndices = DeleteDuplicates[Sort[#] & /@ tensorRules[[All, 1]]];

    vecIndices = Flatten[Position[vec, #] & /@ DeleteCases[vec, 0]];

    For[i = 1, i <= Length[tensorIndices], i++, 
      pos = Position[tensorRules[[All, 1]], tensorIndices[[i]]][[1, 1]];
      val = tensorRules[[pos, 2]];
      ind = tensorRules[[pos, 1]];
      posvec = Position[ind, #];
      If[Length[posvec] > 0,
        posvec = Flatten[posvec][[1]];
        AppendTo[intRules, 
          DeleteCases[ind, ind[[posvec]]] -> Simplify[(-1)^(1 + Length[ind[[1 ;; posvec]]]) vec[[#]] val, OptionValue[Assumptions]]];
      ];
    ] & /@ vecIndices;
    inttensor = SymmetrizedArray[intRules, Table[dim, {i, 1, m - 1}], Antisymmetric[All]];

    If[type == "Form", 
      Return[Tensor2Form[inttensor, xx]];,
      Return[inttensor];
    ];
  ];
]

(* TODO:
-coord change
-induced metric
-Lie bracket
-spin co
-Conf Killing vec eq

-gamma matrices basis
-KSE

-FG solver
*)

End[]                                                                           
                                                                                
EndPackage[]                                                                    
    
                                                                            
Print["benterre's GR package includes exterior algebra of forms and curvature invariants!

For more information about each, type ?ExteriorAlgebra or ?Invariants for more information."']