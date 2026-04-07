# 1. Mixed Finite Element Method for Darcy Equations

We solve the system of equations

$$
    q = -\frac{\kappa}{\mu} \nabla P,
$$

$$
    \nabla \cdot q = f \text{ in } \Omega.
$$

Here $\kappa$ [m $^2$ ] is the permeability of the medium (which may be heterogenous), and $\mu$ [Pa s] is the viscosity of the liquid.

The solver is implemented in [deal.II](https://dealii.org/). For more instructions on installation of deal.II, please follow the steps here: [https://dealii.org/current/readme.html](https://dealii.org/current/readme.html).

# 2. Examples
We consider $\Omega = (0, 1)^2$ [m $^2$ ]. The boundary pressure values are 1 [MPa] on $x = 0$ and 0 [Pa] on $x = 1$. On $y = 0$ and $y = 1$ we consider no flux conditions $q \cdot n = 0$. We also choose $f = 0$.

## 2.1. Homogeneous Media

We consider a homogeneous permeability profile such that $\frac{\kappa}{\mu} = 10^{-4}$. The results are shown below. The pressure profile is shown below.

<div align="center">
<img src='/Images/pressure_homogeneous.png' width='350' height='350'>
</div>

For this simple case, it can be verified that the pressure profile is linear and is given by $p = 10^6(1 - x)$. 

## 2.2. Heterogeneous Media

We now choose a heterogeneous permeability such that $\frac{\kappa}{\mu} = 10^{-8}$ and $\frac{\kappa}{\mu} = 10^{-2}$. The results are shown below.

<div align="center">
<img src='/Images/pressure_heterogeneous.png' width='350' height='350'>
<img src='/Images/flux_heterogeneous.png' width='350' height='350'>
</div>

The white arrows in the flux profile show the flux direction.
