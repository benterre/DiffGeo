# DiffGeo
A Mathematica package to handle everything in differential geometry. From calculating curvature invariants to Hodge star operations on forms, this package can do it all! (That's the aim at least)

So far, the package can be organised into two distinct parts:
 - [Standard GR-invariants](#standard-gr-invariants)
 - [Exterior Algebra](#exterior-algebra)

## Future Additions
- [x] Conformal Killing Vector
- [x] Spin connection
- [ ] Fix optional arguments when specifying part of the inputs. eg `Lie[vec, g, xx, "Down" -> {1, 2}]` works but `Lie[vec, g, "Down" -> {1, 2}]` confuses the option for `xx`
- [ ] Fix https://github.com/benterre/DiffGeo/issues/1 AntiSelfDualQ returns False on array input
- [x] Allow Form2Tensor to convert tensor-valued n-forms
- [ ] Lie Bracket
- [x] Codifferential of forms
- [x] Laplace-Beltrami Operator
- [ ] Vielbeins
- [ ] Allow Form2Tensor to convert directly to non-coordinate basis
- [ ] Gamma matrices
- [ ] Killing Spinors
- [ ] Fefferman-Graham Gauge

## Standard GR-invariants
> [!NOTE]
> For these functions, the metric should be given as a $d\times d$ array and the list of coordinates as a $d$ array.

- [Christoffel symbols](#christoffelg-xx)
- [Riemann Tensor](#riemanntensorg-xx)
- [Ricci Tensor](#riccitensorg-xx)
- [Ricci Scalar](#ricciscalarg-xx)
- [Kretschmann Scalar](#kretschmannscalarg-xx)
- [Weyl Tensor](#weyltensorg-xx)
- [Euler Scalar](#eulerscalarg-xx)
- [Cotton Tensor](#cottontensorg-xx)
- [Schouten Tensor](#schoutentensorg-xx)
- [Lovelock Invariants](#lovelockk-ng-xx)

### Christoffel[g, xx]
Calculates the Christoffel symbols for the associated metric `g` and coordinates `xx`.
Conventions are such that its components are given by spacetime indices $`\Gamma^\mu{}_{\nu\rho}`$ in that same order.
In other words,

`Christoffel[g, xx][[\[Mu],\[Nu],\[Rho]]]`=$`\Gamma^{\mu}{}_{\nu\rho}`$.

These are computed from the metric using

$$\Gamma^\mu_{\nu\rho}=\frac{1}{2}g^{\mu\sigma}(\partial_\nu g_{\rho\sigma}+\partial_\rho g_{\nu\sigma}-\partial_\sigma g_{\nu\rho})$$

### RiemannTensor[g, xx]
Calculates the Riemann tensor for the associated metric `g` and coordinates `xx`.
Conventions are such that all spacetime indices are down $R_{\mu\nu\rho\sigma}$.
The components are computed from the metric using

```math
R_{\mu\nu\rho\sigma} = g_{\mu\alpha} \partial_\rho\Gamma^\alpha{}_{\sigma\nu} - \partial_\sigma \Gamma^\alpha{}_{\rho\nu} + \Gamma^{\alpha}{}_{\rho\lambda}\Gamma^\lambda{}_{\sigma\nu} - \Gamma^{\alpha}{}_{\sigma\lambda}\Gamma^{\lambda}{}_{\rho\nu}
```

### RicciTensor[g, xx]
Calculates the Ricci tensor for the associated metric `g` and coordinates `xx`.
Conventions are such that all spacetime indices are down $R_{\mu\nu}$.

The components are computed from the metric using

$$R_{\mu\nu}=g^{\rho\sigma} R_{\rho\mu\sigma\nu}$$

### RicciScalar[g, xx]
Calculates the Ricci scalar for the associated metric `g` and coordinates `xx`.

The components are computed from the metric using

$$R=g^{\mu\nu}R_{\mu\nu}$$

### KretschmannScalar[g, xx]
Calculates the Kretschmann scalar for the associated metric `g` and coordinates `xx`.
Its expression is identical to the third Lovelock invariant at level 2, `I32`.

The components are computed from the metric using

$$K_1 = R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}$$

### WeylTensor[g, xx]
Calculates the Weyl tensor for the associated metric `g` and coordinates `xx`.
Conventions are such that all spacetime indices are down $C_{\mu\nu\rho\sigma}$.

The components are computed from the metric using

$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma} + \frac{1}{d-2} (R_{\mu\sigma}g_{\nu\rho} - R_{\mu\rho}g_{\nu\sigma} + R_{\nu\rho}g_{\mu\sigma} - R_{\nu\sigma}g_{\mu\rho}) + \frac{1}{(d-1)(d-2)} R(g_{\mu\rho}g_{\nu\sigma}- g_{\mu\sigma}g_{\nu\rho})$$

### EulerScalar[g, xx]
Calculates the Euler scalar for the associated metric `g` and coordinates `xx`.

The components are computed from the metric through the Weyl tensor as

$$K_3 = -C_{\mu\nu\rho\sigma}C^{\mu\nu\rho\sigma} + 2R_{\mu\nu}R^{\mu\nu} - \frac{2}{3}R^2$$

### CottonTensor[g, xx]
Calculates the Cotton tensor for the associated metric `g` and coordinates `xx`.
The components are computed from the metric using

$$C_{\mu\nu\rho} = \nabla_\rho R_{\mu\nu} - \nabla_{\nu}R_{\mu\rho}  +\frac{1}{2(d-1)}(\nabla_{\nu} R g_{\mu\rho} - \nabla_{\rho} R g_{\mu\nu})$$

### SchoutenTensor[g, xx]
Calculates the Schouten tensor for the associated metric `g` and coordinates `xx`.

The components are computed from the metric using

$$P_{\mu\nu} =\frac{1}{d-2} \left( R_{\mu\nu} - \frac{R}{2(d-1)} g_{\mu\nu} \right)$$

### Lovelock[k, n][g, xx]
Returns the `n`th Lovelock invariant at order `k` for the associated metric `g` and coordinates `xx`.
For `n` greater than the number of Lovelock invariant at a given order `k`, `Lovelock[n,k]` might return additional invariants built from the Weyl tensor (marked in gray below).

| Order | Command | Expression |
| :---: | :--- | :--- |
| `k=1` | `L[1,1]` | $L^1_1=R$ |
| `k=2` | `L[2,1]` | $L^2_1=R^2$ |
| | `L[2,2]` | $L^2_2=R_{\mu\nu}R^{\mu\nu}$ |
| | `L[2,3]` | $L^2_3=R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}=K_1$ |
| | `L[2,4]` | $\color{gray}{C_{\mu\nu\rho\sigma}C^{\mu\nu\rho\sigma}}$ |
| `k=3` | `L[3,1]` | $L^3_1=R^3$ |
| | `L[3,2]` | $L^3_2=R R_{\mu\nu}R^{\mu\nu}=R L^2_2$ |
| | `L[3,3]` | $L^3_3=R_{\mu\alpha} R^{\nu}{}_\mu R^{\alpha\mu}$ |
| | `L[3,4]` | $L^3_4=R_{\mu\alpha} R_{\nu\beta} R^{\nu\mu\alpha\beta}$ |
| | `L[3,5]` | $L^3_5=R R_{\mu\nu\alpha\beta} R^{\mu\nu\alpha\beta}=R K_1$ |
| | `L[3,6]` | $L^3_6=R_{\nu\alpha} R_{\beta\gamma\epsilon}{}^\nu R^{\beta\gamma\epsilon\alpha}$ |
| | `L[3,7]` | $L^3_7=R_{\mu\nu\alpha\beta} R^{\mu\nu}{}_{\gamma\epsilon} R^{\alpha\beta\gamma\epsilon}$ |
| | `L[3,8]` | $`L^3_8=R_{\mu\nu\alpha\beta} R^{\mu}{}_{\gamma}{}^{\alpha}{}_{\epsilon} R^{\nu\gamma\beta\epsilon}`$ |
| | `L[3,9]` | $`\color{gray}{C_{\mu\nu\alpha\beta}R^{\mu\alpha}R^{\nu\beta}}`$ |
| | `L[3,10]` | $`\color{gray}{C_{\mu\nu\alpha\beta}C^{\alpha\beta\rho\nu}R^\mu{}_\rho}`$ |
| | `L[3,11]` | $`\color{gray}{C_{\mu\nu\alpha\beta}C^{\alpha\beta\rho\sigma}C_{\rho\sigma}{}^{\mu\nu}}`$ |

### LovelockAll[g, xx]
Returns an array containing all Lovelock Invariant available in this package.

- `LovelockAll[g, xx, "Extend"->True]`: includes invariants built from the Weyl tensor (default)
- `LovelockAll[g, xx, "Extend"->False]`: doesn't include invariants built from the Weyl tensor

### GRInvariants[g, xx]
Returns an array containing all the scalar curvature invariant available in this package.

 - `GRInvariants[g, xx, "Minimal"->True]`: only returns a basis of multiplicatively independent invariants
 - `GRInvariants[g, xx, "Minimal"->False]`: returns all invariants (default)

## Exterior Algebra
> [!NOTE]
> For this part of the package to work properly, avoid changing the definition of the exterior derivative `d` as well as the wedge product `Wedge`

Based on the EDC360 package. The exterior algebra part of this package defines the exterior derivative as `d`. The Wedge product can be used on forms as `d[x]\[Wedge]d[y]`.

- [Form2Metric](#form2metric)
- [Form2Tensor](#form2tensor)
- [Tensor2Form](#tensor2form)
- [HodgeStar](#hodgestar)
- [cod](#cod)
- [LaplaceBeltrami](#laplacebeltrami)
- [SelfDualQ](#selfdualq)
- [AntiSelfDualQ](#antiselfdualq)
- [Lie](#lie)
- [KillingQ](#killingq)
- [Cov](#cov)
- [ConfKillingQ](#confkillingq)
- [InteriorProduct](#interiorproduct)
- [SpinConnection](#spinconnection)

### Form2Metric
`Form2Metric[g,xx]`: Converts a metric, (abusively) written using squared one-forms, into its matrix representation, using the set of coordinates `xx`.

> Example
> ```
> g = d[r]^2 + r^2 d[\[Theta]]^2;
> xx = {r,\[Theta]};
>
> Form2Metric[g, xx]
> ```
> `{{1,0},{0,r^2}}`

### Form2Tensor
`Form2Tensor[form, xx]`: Converts a d-form `form`, represented using the exterior derivative `d`, into its tensor components, using `xx` as a list of coordinates.

`Form2Tensor[form]` will try to use a globally defined variable called `xx` as the set of coordinates and then call `Form2Tensor[form, xx]`.

Optional Arguments:
  - `Assumptions->None`: Array of assumptions used in Simplify
  - `"Type"->String`. When equal to 'SymmetrizedArray' returns the tensor as a fully anti-symmetric SymmetrizedArray. When set to 'Normal' returns a standard array. (defaults to 'Normal')

> Example
> ```
> form = d[x]\[Wedge]d[y];
> xx = {x,y};
>
> Form2Tensor[form, xx]
> ```
> `{{0,1},{-1,0}}`

### Tensor2Form
`Tensor2Form[tensor, xx]`: Converts a tensor `tensor` of any rank into its form representation using `xx` as a list of coordinates.

`Tensor2Form[tensor]` will try to use a globally defined variable called `xx` as the set of coordinates and then call `Tensor2Form[tensor, xx]`.

If the original tensor has a symmetric part, it is automatically set to zero.

Optional Arguments:
  - `Assumptions->None`: Array of assumptions used in Simplify

> Example
> ```
> tensor = {{0,1},{-1,0}};
> xx = {w,z};
>
> Tensor2Form[tensor, xx]
> ```
> `d[w]\[Wedge]d[z]`

### HodgeStar
`HodgeStar[tensor, g]`: Returns the Hodge dual of the tensor `tensor`, with respect to the metric `g`.

`HodgeStar[form, g, xx]`: Converts the form `form` to its tensor representation then calls `HodgeStar[tensor, metric]` using the list of coordinates `xx`.

If metric isn't specified, HodgeStar will try to use a globally defined variable called `g` as the metric.

Optional Arguments:
 - `Assumptions->None`: Array of assumptions used in Simplify

> Example
> ```
> xx = {r,\[Theta],\[Phi]};
> g = Form2Metric[d[r]^2 + r^2 (d[\[Theta]]^2 + Sin[\[Theta]]^2 d[\[Phi]]^2), xx];
>
> form = d[r]\[Wedge]d[\[Theta]];
> 
> HodgeStar[form, g, xx]
> HodgeStar[form, g, xx, Assumptions -> {r > 0, \[Theta] > 0, \[Theta] < \[Pi]}]
> ```
> `(Sqrt[Abs[r^4 Sin[\[Theta]]^2]] d[\[Phi]])/r^2`
>
> `d[\[Phi]] Sin[\[Theta]]`

### cod
`cod[form, g, xx]`: Returns the codifferential of the `form`, given a metric `g` and coordinates `xx`.

If $`\omega`$ is a $`k`$-form on a $`d`$-dimensional manifold with metric signature $`s`$, then its codifferential is calculated using

```math
d^\dagger\omega = (-1)^{d(k+1)+1}s\star d\star\omega
```

Optional Arguments:
 - `Assumptions->None`: Array of assumptions used in Simplify


### LaplaceBeltrami
`LaplaceBeltrami[form, g, xx]`: Returns the Laplace-Beltrami operator acted on the `form`, given a metric `g` and coordinates `xx`.

If $`\omega`$ is a $`k`$-form on a manifold with metric $`g`$, then Laplace-Beltrami operator is defined as

```math
\Delta\omega = (d^\dagger d + d d^\dagger)\omega = (d+d^\dagger)^2\omega
```

Optional Arguments:
 - `Assumptions->None`: Array of assumptions used in Simplify

### SelfDualQ
`SelfDualQ[form, g, xx]`: Return `True` if `form` is self-dual wrt the metric `g` and `False` otherwise.

Accepts rank-r tensor and r-form representation for the input. The former doesn't require the coordinates `xx` to be specified.

`SelfDualQ[form]` will try to use globally defined variables called `xx` as the set of coordinates and `g` as the metric and then calls `SelfDualQ[form, g, xx]`.

Optional Arguments:
  - `Assumptions->None`: Array of assumptions used in Simplify

> Example
> ```
> xx = {x1, x2, x3, x4};
> g = IdentityMatrix[4];
> 
> form1 = d[x1]\[Wedge]d[x2] - d[x3]\[Wedge]d[x4];
> form2 = d[x1]\[Wedge]d[x2] + d[x3]\[Wedge]d[x4];
> 
> SelfDualQ[form1, g, xx]
> SelfDualQ[form2, g, xx]
> ```
> `False`
>
> `True`

### AntiSelfDualQ
`AntiSelfDualQ[form, g, xx]`: Return `True` if `form` is anti-self-dual wrt the metric `g` and `False` otherwise.

Accepts rank-r tensor and r-form representation for the input. The former doesn't require the coordinates `xx` to be specified.

`AntiSelfDualQ[form]` will try to use globally defined variables called `xx` as the set of coordinates and `g` as the metric and then calls `AntiSelfDualQ[form, g, xx]`.

Optional Arguments:
  - `Assumptions->None`: Array of assumptions used in Simplify

> Example
> ```
> xx = {x1, x2, x3, x4};
> g = IdentityMatrix[4];
> 
> form1 = d[x1]\[Wedge]d[x2] - d[x3]\[Wedge]d[x4];
> form2 = d[x1]\[Wedge]d[x2] + d[x3]\[Wedge]d[x4];
> 
> AntiSelfDualQ[form1, g, xx]
> AntiSelfDualQ[form2, g, xx]
> ```
> `True`
>
> `False`

### Lie
`Lie[vec, tensor, xx]`: Returns the Lie derivative of `tensor` with respect to the vector `vec` using the set of coordinates `xx`.
A unit vector in the coordinates can be used as vector input by specifying an integer (in the range 1 to `Length[xx]`) instead of `vec`.
A form of degree r can be used in lieu of `tensor`, after which it is converted to a tensor of rank r.

`Lie[vec, tensor]` will try to use a globally defined variable called `xx` as the set of coordinates and then call `Lie[vec, tensor, xx]`.

Let `tensor` be a $(p,q)$-tensor, whose components are denoted by $`T^{\mu_1\cdots \mu_p}{}_{\nu_1\cdots\nu_q}`$, and let's denote the components of `vec` by $`\xi^\mu`$. Then `Lie[vec, tensor, xx]` is also a $(p,q)$-tensor with components

```math
(\mathcal{L}_{\xi}T)^{\mu_1\cdots \mu_p}{}_{\nu_1\cdots\nu_q} = \xi^\rho\partial_\rho T^{\mu_1\cdots \mu_p}{}_{\nu_1\cdots\nu_q} - \partial_\rho\xi^{\mu_1}T^{\rho\mu_2\cdots \mu_p}{}_{\nu_1\cdots\nu_q} - \cdots - \partial_\rho\xi^{\mu_p}T^{\mu_\cdots \mu_{p-1}\rho}{}_{\nu_1\cdots\nu_q}
+ \partial_{\nu_1}\xi^{\rho}T^{\mu_1\cdots \mu_p}{}_{\rho\nu_2\cdots\nu_q} + \cdots + \partial_{\nu_q}\xi^{\rho}T^{\mu_1\cdots \mu_p}{}_{\nu_1\cdots\nu_{q-1}\rho}
```

Optional Arguments:
  - `Assumptions->None`: Array of assumptions used in Simplify
  - `"Up"->None`: Array containing the set of indices in 'tensor' which are contravariant (or up)
  - `"Down"->Range[rank]`: Array containing the set of indices in 'tensor' which are covariant (or down)
  - `Type->"Other"`: Used to specify the output format. `Type->"Form"` will output the tensor as a rank-r form (defaults to that when using a form as input)

> Example
> ```
> xx = {r, \[Theta], \[Phi]};
> g = DiagonalMatrix[{1, r^2, r^2 Sin[\[Theta]]^2}];
> 
> vec = {r^2, \[Theta], 0};
> tensor = {{r, \[Theta]}, {Sin[\[Theta]], Cos[\[Theta]]}};
> 
> Lie[vec, tensor]
> ```
> `{{r^2, \[Theta]}, {\[Theta] Cos[\[Theta]], -\[Theta] Sin[\[Theta]]}}`

> Example
> ```
> xx = {r, \[Theta], \[Phi]};
> g = DiagonalMatrix[{1, r^2, r^2 Sin[\[Theta]]^2}];
> 
> Lie[1, g]
> Lie[3, g]
> ```
> `{{0, 0, 0}, {0, 2 r, 0}, {0, 0, 2 r Sin[\[Theta]]^2}}`
> 
> `{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}`

### KillingQ
`KillingQ[vector, g, xx]`: Returns `True` if the vector is Killing with respect to the metric `g` given the set of coordinates `xx`, and `False` otherwise.

`KillingQ[vector, g]` will try to use a globally defined variable called `xx` as the set of coordinates and then calls `KillingQ[vector, g, xx]`.

`KillingQ[vector]` will try to use a globally defined variable called `g` as the metric and then calls `KillingQ[vector, g]`.

A vector field $`V=V^\mu\partial_\mu`$ is Killing if the Lie derivative of the metric along it vanishes, $`\mathcal{L}_{V}g=0`$.

Optional Arguments:
  - `Assumptions->None`: Array of assumptions used in Simplify

### Cov
`Cov[vec, tensor, g, xx]`: Returns the covariant derivative of `tensor` with respect to the vector `vec` using the Christoffel symbols built from the metric `g` and set of coordinates `xx`.
A unit vector in the coordinates can be used as vector input by specifying an integer (in the range 1 to `Length[xx]`) instead of `vec`.
A form of degree r can be used in lieu of `tensor`, after which it is converted to a tensor of rank r.

`Cov[vec, tensor, g]` will try to use a globally defined variable called `xx` as the set of coordinates and then calls `Lie[vec, tensor, g, xx]`.
`Cov[vec, tensor]` will try to use a globally defined variable called `g` as the metric and then calls `Lie[vec, tensor, g]`.

Let `tensor` be a $(p,q)$-tensor, whose components are denoted by $`T^{\mu_1\cdots \mu_p}{}_{\nu_1\cdots\nu_q}`$, and let's denote the components of `vec` by $`\xi^\mu`$. Then `Cov[vec, tensor, xx]` is also a $(p,q)$-tensor with components

```math
(\nabla_{\xi}T)^{\mu_1\cdots \mu_p}{}_{\nu_1\cdots\nu_q} = \xi^\rho\partial_\rho T^{\mu_1\cdots \mu_p}{}_{\nu_1\cdots\nu_q} + \xi^{\rho}\Gamma^{\mu_1}{}_{\rho\sigma}T^{\sigma\mu_2\cdots \mu_p}{}_{\nu_1\cdots\nu_q} + \cdots + \xi^{\rho}\Gamma^{\mu_p}{}_{\rho\sigma}T^{\mu_\cdots \mu_{p-1}\sigma}{}_{\nu_1\cdots\nu_q}
- \xi^{\rho}\Gamma^{\sigma}{}_{\rho\nu_1}T^{\mu_1\cdots \mu_p}{}_{\sigma\nu_2\cdots\nu_q} - \cdots - \xi^{\rho}\Gamma^{\sigma}{}_{\rho\nu_q}T^{\mu_1\cdots \mu_p}{}_{\nu_1\cdots\nu_{q-1}\sigma}
```

Optional Arguments:
  - `Assumptions->None`: Array of assumptions used in Simplify
  - `"Up"->None`: Array containing the set of indices in 'tensor' which are contravariant (or up)
  - `"Down"->Range[rank]`: Array containing the set of indices in 'tensor' which are covariant (or down)
  - `Type->"Other"`: Used to specify the output format. `Type->"Form"` will output the tensor as a rank-r form (defaults to that when using a form as input)

> Example
> ```
> 
> ```

### ConfKillingQ
`ConfKillingQ[vector, g, xx]`: Returns `True` if the vector is conformal Killing with respect to the metric `g` given the set of coordinates `xx`, and `False` otherwise.

`ConfKillingQ[vector, g]` will try to use a globally defined variable called `xx` as the set of coordinates and then calls `ConfKillingQ[vector, g, xx]`.

`ConfKillingQ[vector]` will try to use a globally defined variable called `g` as the metric and then calls `ConfKillingQ[vector, g]`.

A vector field $`V=V^\mu\partial_\mu`$ is conformal Killing if it obeys the equation
```math
\mathcal{L}_V g=\frac{2}{d}\nabla_\mu V^\mu g
```

Optional Arguments:
  - `Assumptions->None`: Array of assumptions used in Simplify

### InteriorProduct
`InteriorProduct[vec, form, xx]`: Returns the interior product of the form with respect to the vector `vec` given the set of coordinates `xx`.

`InteriorProduct[vec, tensor, xx]`: Returns the interior product of the tensor with respect to the vector `vec` given the set of coordinates `xx`.

Optional Arguments:
  - `Assumptions->None`: Array of assumptions used in Simplify
  - `Type->"Other"`: Used to specify the ouput format. `Type->"Form"` will output a rank r-1 form (defaults to that when using rank-r form as input)


### SpinConnection
`SpinConnection[e, xx]`: Returns the spin connection one-form $`\omega^{ab}`$ with respect to the array of vielbein one-forms `e` and coordinate list `xx`.
The vielbeins `e` can also be specified using a matrix with non-coordinate components first, i.e. $`e^a{}_\mu`$.

Its components are defined via
```math
\omega_\mu{}^{ab}=e^a{}_\nu(\partial_\mu e^{\nu b}+\Gamma^{\nu}{}_{\mu\rho}e^{\rho b})
```

Optional Arguments:
  - `Assumptions->None`: Array of assumptions used in Simplify

> Example
> ```
> xx = {\[Theta], \[CurlyPhi], \[Tau], t};
> e = {Sin[\[CurlyPhi] + \[Tau]] d[\[Theta]] + Cos[\[Theta]] Sin[\[Theta]] Cos[\[CurlyPhi] + \[Tau]] d[\[CurlyPhi]] \[Minus] Cos[\[Theta]] Sin[\[Theta]] Cos[\[CurlyPhi] + \[Tau]] d[\[Tau]], -Cos[\[CurlyPhi] + \[Tau]] d[\[Theta]] + Cos[\[Theta]] Sin[\[Theta]] Sin[\[CurlyPhi] + \[Tau]] d[\[CurlyPhi]] - Cos[\[Theta]] Sin[\[Theta]] Sin[\[CurlyPhi] + \[Tau]] d[\[Tau]], Sin[\[Theta]]^2 d[\[CurlyPhi]] + Cos[\[Theta]]^2 d[\[Tau]], \[Beta] d[t]};
>
> SpinConnection[e, xx]
> ```
> `{{0, Cos[\[Theta]]^2 d[\[Tau]] + d[\[CurlyPhi]] Sin[\[Theta]]^2, Cos[\[Tau] + \[CurlyPhi]] d[\[Theta]] + Cos[\[Theta]] (d[\[Tau]] - d[\[CurlyPhi]]) Sin[\[Theta]] Sin[\[Tau] + \[CurlyPhi]], 0}, {-Cos[\[Theta]]^2 d[\[Tau]] - d[\[CurlyPhi]] Sin[\[Theta]]^2, 0, -Cos[\[Theta]] Cos[\[Tau] + \[CurlyPhi]] d[\[Tau]] Sin[\[Theta]] \+ Cos[\[Theta]] Cos[\[Tau] + \[CurlyPhi]] d[\[CurlyPhi]] Sin[\[Theta]] + d[\[Theta]] Sin[\[Tau] + \[CurlyPhi]], 0}, {-Cos[\[Tau] + \[CurlyPhi]] d[\[Theta]] + Cos[\[Theta]] (-d[\[Tau]] + d[\[CurlyPhi]]) Sin[\[Theta]] Sin[\[Tau] + \[CurlyPhi]], Cos[\[Theta]] Cos[\[Tau] + \[CurlyPhi]] d[\[Tau]] Sin[\[Theta]] - 1/2 Cos[\[Tau] + \[CurlyPhi]] d[\[CurlyPhi]] Sin[2 \[Theta]] - d[\[Theta]] Sin[\[Tau] + \[CurlyPhi]], 0, 0}, {0, 0, 0, 0}}`
