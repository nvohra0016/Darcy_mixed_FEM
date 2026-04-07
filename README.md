# 1. Mixed Finite Element Method for Darcy Equations

We solve the system of equations

$$
    q = -\frac{\kappa}{\mu} \nabla P,
$$

$$
    \nabla \cdot q = 0 \text{ in } \Omega.
$$

Here $\kappa$ [m $^2$ ] is the permeability of the medium (which may be heterogenous), and $\mu$ [Pa s] is the viscosity of the liquid.

The solver is implemented in [deal.II](https://dealii.org/). For more instructions on installation of deal.II, please follow the steps here: [https://dealii.org/current/readme.html](https://dealii.org/current/readme.html).

# 2. Examples

## 2.1. Homogeneous Media

## 2.2. Heterogeneous Media

<div align="center">
<img src='/Images/flux_heterogeneous.png' width='350' height='350'>
<img src='/Images/pressure_heterogeneous.png' width='350' height='350'>
</div>
