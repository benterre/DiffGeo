(* ::Package:: *)

(* Author: Benjamin Suzzoni
Organization: University of Southampton *)

(* This package contains definitions for various curvature invariants of General Relativity *)
BeginPackage["GRBenterre`"]

Christoffel::usage = "Calculates the Christoffel symbols for the associated metric g and coordinates xx.
Conventions are such that its components are given by spacetime indices \!\(\*SubscriptBox[SuperscriptBox[\[CapitalGamma],\[Mu]], \[Nu]\[Rho]]\) in that same order.
These are computed from the metric using

\!\(\*SubscriptBox[SuperscriptBox[\[CapitalGamma],\[Mu]], \[Nu]\[Rho]]=\*FractionBox[\(1\), \(2\)]\\ \*SuperscriptBox[\(g\),\[Mu]\[Sigma]]\\ ( \*SubscriptBox[\(\[PartialD]\), \(\[Nu]\)]\*SubscriptBox[\(g\),\[Rho]\[Sigma]]+\*SubscriptBox[\(\[PartialD]\), \(\[Rho]\)]\*SubscriptBox[\(g\),\[Nu]\[Sigma]]-\*SubscriptBox[\(\[PartialD]\), \(\[Sigma]\)]\*SubscriptBox[\(g\),\[Nu]\[Rho]]) \)"

RiemannTensor::usage = "Calculates the Riemann tensor for the associated metric g and coordinates xx.
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

ExteriorAlgebra::usage = "Based on the EDC360 package"

Form2Metric::usage = "Form2Metric[g,xx] Converts a metric (abusively) written using squared one-forms,
eg. g=d[r]^2+r^2d[\[Theta]]^2 and xx={r,\[Theta]},
into its matrix representation."

HodgeStar::usage = "HodgeStar[tensor, metric]: returns the Hodge dual of the tensor, with respect to the {d,d} matrix array.
HodgeStar[form, metric, xx]: converts the form to its tensor representation then calls HodgeStar[tensor, metric]

Optional Arguments:
 - Assumptions->List. Used in all Simplify[] calls."
(*Star::usage = "Star[form, g, xx] is a custom operator for HodgeStar.";*)
HodgeStar::missingInput = "Missing required coordinates array xx when converting to form to tensor.
Common usage HodgeStar[form, g, xx]."

SelfDualQ::usage = "SelfDualQ[form, metric, coords]: return true if form is self-dual wrt the metric and false otherwise.
Accepts rank-r tensor and r-form representation for the input form. The fromer doesn't require coords to be specified."
AntiSelfDualQ::usage = "AntiSelfDualQ[form, metric, coords]: return true if form is anti-self-dual wrt the metric and false otherwise.
Accepts rank-r tensor and r-form representation for the input form. The fromer doesn't require coords to be specified."

Form2Tensor::usage = "Form2Tensor[form, xx]: Converts a d-form, represented using d[x]\[Wedge]d[y]... into its tensor components. xx is a list of coordinates.

Optional Arguments:
 - \"Type\"->String. When equal to 'SymmetrizedArray' returns the tensor as a fully anti-symmetric SymmetrizedArray. When set to 'Normal' returns a standard array. (defaults to 'Normal')"

Tensor2Form::usage = "Tensor2Form[tensor, xx]: Converts a tensor of any dimension into its form representation, using the d[] function. xx is a list of coordinates.

If the original tensor has a symmetric part, it is automatically set to zero."
Tensor2Form::inputForm = "The input variable cannot be a differential form. Please provide a suitable tensor object, as an Array."

(* ScalarEDCode functions *)
d::usage = "d[form] applies the exterior derivative to the differential form 'form'.";
Wedge::usage = "Wedge[a, b] returns the wedge product of two differential forms 'a' and 'b'.";

Lie::usage = ""
Lie::vectorDim = "Dimensions of vector don't match that of the coordinates. Please provide a correct vector."

(*Needs["Notation`"];*)

Begin["`Private`"]

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

(* Other functions *)
Options[Form2Metric] = {Assumptions->None};
Form2Metric[g_, xx_] := Simplify[Table[If[i === j, 1, 1/2] Coefficient[g, d[xx[[i]]] d[ xx[[j]]]], {i, 
    Length[xx]}, {j, Length[xx]}],OptionValue[Assumptions]];

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

(* Define the custom Star function that will wrap HodgeStar
Star[form_, g_, xx_] := HodgeStar[form, g, xx];

(* Define the custom notation for the \[Star] operator *)
InfixNotation[ParsedBoxWrapper[
  SubsuperscriptBox["\[Star]", "xx_", "g_"]], 
  Star];*)

Options[SelfDualQ] = {Assumptions->None};
SelfDualQ[formT_, gT_, xx_, OptionsPattern[]] := Simplify[formT-HodgeStar[formT,gT,xx],OptionValue[Assumptions]]===0;
Options[AntiSelfDualQ] = {Assumptions->None};
AntiSelfDualQ[formT_, gT_, xx_, OptionsPattern[]] := Simplify[formT+HodgeStar[formT,gT,xx],OptionValue[Assumptions]]===0;


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

Options[Lie] = {Assumptions->None, "Up"->None, "Down"->None};
Lie[vecT_, tensorT_, xx_, OptionsPattern[]] := Module[
  {r, dim, up, down, lie, tensor, vec},

  up = OptionValue["Up"];
  down = OptionValue["Down"];

  If[FormDegree[tensorT] != 0,
    r = FormDegree[tensorT];
    up = {};
    down = Range[r];
    tensor = Form2Tensor[tensorT, xx],

    tensor = tensorT;
    r = TensorRank[tensor];
    If[StringContainsQ[ToString[r], "TensorRank"], r=0];
  ];

  If[up==None && down==None,
    up={};
    down=Range[r];
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
  
  If[rank > 0,
    lie += -Total[TensorTranspose[TensorContract[TensorProduct[Grad[vec, xx], tensor], {2, 2 + #}], Join[{#}, Complement[Range[r], {#}]]] & /@ up];
    lie += Total[TensorTranspose[TensorContract[TensorProduct[Grad[vec, xx], tensor], {1, 2 + #}], Join[{#}, Complement[Range[r], {#}]]] & /@ down];
  ];
  Return[lie]
]

(* TODO:
-inner prod
-Killing vec eq
-Conf Killing vec eq
*)


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

End[]                                                                           
                                                                                
EndPackage[]                                                                    
    
                                                                            
Print["benterre's GR package includes exterior algebra of forms and curvature invariants!

For more information about each, type ?ExteriorAlgebra or ?Invariants for more information."']