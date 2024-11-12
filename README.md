# DiffGeo
A Mathematica package to handle everything in differential geometry. From calculating curvature invariants to Hodge star operations on forms, this package can do it all! (That's the aim at least)

## Standard GR-invariants
> [!NOTE]
> For these functions, the metric should be given as a $d\times d$ array and the list of coordinates as a $d$ array.

- [Christoffel symbols](#christoffelg-xx)
- [Riemann Tensor](#riemanntensorg-xx)
- [Ricci Tensor](#riccitensorg-xx)
- [Ricci Scalar](#ricciscalarg-xx)

### Christoffel[g, xx]
Calculates the Christoffel symbols for the associated metric `g` and coordinates `xx`.
Conventions are such that its components are given by spacetime indices $\Gamma^\mu{}_{\nu\rho}$ in that same order.
In other words,

`Christoffel[g, xx][[\[Mu],\[Nu],\[Rho]]]`=$\Gamma^{\mu}_{\nu\rho}$.

These are computed from the metric using

$$\Gamma^\mu_{\nu\rho}=\frac{1}{2}g^{\mu\sigma}(\partial_\nu g_{\rho\sigma}+\partial_\rho g_{\nu\sigma}-\partial_\sigma g_{\nu\rho})$$

### RiemannTensor[g, xx]
Calculates the Riemann tensor for the associated metric `g` and coordinates `xx`.
Conventions are such that all spacetime indices are down $R_{\mu\nu\rho\sigma}$.
The components are computed from the metric using

$$R_{\mu\nu\rho\sigma} = g_{\mu\alpha} \partial_\rho\Gamma^\alpha{}_{\sigma\nu} - \partial_\sigma \Gamma^\alpha{}_{\rho\nu} + \Gamma^{\alpha}{}_{\rho\lambda}\Gamma^\lambda{}_{\sigma\nu} - \Gamma^{\alpha}{}_{\sigma\lambda}\Gamma^{\lambda}{}_{\rho\nu}$$

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
For `n` greater than the number of Lovelock invariant at a given order `k`, `Lovelock[n,k]` might return additional invariants built from the Weyl tensor.
