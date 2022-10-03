# `internal-fluid-flow` Toolbox for GNU-Octave

[![DOI](https://zenodo.org/badge/509427410.svg)](https://zenodo.org/badge/latestdoi/509427410)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/aumpierre-unb/Internal-Fluid-Flow-for-GNU-Octave)

## Installing and Loading `internal-fluid-flow`

```dotnetcli
# use this call to install version 0.2.0, or modify the command line for match the version
pkg install https://github.com/aumpierre-unb/Internal-Fluid-Flow-for-GNU-Octave/archive/refs/tags/v0.2.0.tar.gz
pkg load internal-fluid-flow
```

## Citation of `internal-fluid-flow`

You can cite all versions (both released and pre-released), by using
[DOI 105281/zenodo.6960263](https://doi.org/10.5281/zenodo.6960263).

This DOI represents all versions, and will always resolve to the latest one.

For citation of the last released version of `internal-fluid-flow`, please check CITATION file at the [maintainer's repository](https://github.com/aumpierre-unb/Internal-Fluid-Flow-for-GNU-Octave).

---

The following is a very short introduction to the steady internal flow of an incompressible and inviscid fluid and to the `internal-fluid-flow` Toolbox for GNU Octave.

Internal flow is a pretty extensive topic in fluid mechanics and there are a lot of important and interesting observations related to it that are not taken into account in this text, because they have no direct impact the computation performed by the functions in this package. Our focus here is a small set of equations that described the phenomenon and are required to solve problems on internal fluid flow.

This text is divided in two main sections: The Theory and The `internal-fluid-flow` Toolbox.

## The Theory

### The Bernoulli Equation

The Bernoulli equation is an expression of the mechanical energy balance for a very particular situation:

- internal steady flow of an
- incompressible inviscid fluid, where
- friction effects and tube fittings can be neglected.

For such a case, the mechanical energy is conserved, and for any two points 1 and 2 we have

$$
{\rho v_2^2 \over 2} + \rho g z_2 + p_2 =
{\rho v_1^2 \over 2} + \rho g z_1 + p_1
$$

or

$$
{v_2^2 \over 2g}+z_2+{p_2 \over \rho g}=
{v_1^2 \over 2g}+z_1+{p_1 \over \rho g}
$$

where

- $\rho$ is the fluid's density,
- *v* is the speed flow,
- *g* is the gravitational acceleration,
- *z* is the elevation, and
- *p* is the static pressure.

### Head Loss

The flow of viscous fluids is accompanied of energy dispersion, which can be measured as pressure drop or, equivalently, as head loss *h*, by the Darcy-Weisbach equation,

$$
h=f{v^2 \over 2g} {L \over D}
$$

where *f* is the Darcy friction factor, *L* is the pipe's length and *D* is the pipe's hydraulic diameter,

$$
D={4A \over P}
$$

where *A* is the cross-sectional area of the flow and *P* is the wet perimeter of the cross-section. *f* is described as a function of the Reynolds number,

$$
Re={\rho vg \over \mu}
$$

and the pipe's relative roughness,

$$
\varepsilon={k \over D}
$$

where

- $\mu$ is the fluid's dynamic viscosity and
- *k* is the pipe's[ internal surface] roughness.

The Reynolds number *Re*, the Darcy friction factor *f*, and the relative roughness $\varepsilon$ completely describe the internal flow of incompressible viscous fluids, for both laminar and turbulent regimes. Usually, *f* is given as a function of *Re* and $\varepsilon$.

The simplest problems on internal fluid flow consist on computing one of them given the two other. More complex situations arise when only one or none of those variables is known. Instead, dimensional variables involved are given. However not always, in most cases iterative computation is required.

### Laminar Flow and Turbulent Flow

For laminar flow, *Re* < 2,500 (typically), the Darcy friction factor is given by the Poiseuille condition,

$$
f={64 \over Re}
$$

For turbulent flow, *Re* > 2,500 (typically), the Darcy friction factor is given implicitly by the Colebrook-White equation,

$$
{1 \over \sqrt{f}}=2 \mathrm{log} {1 \over\displaystyle {3.7 \over \varepsilon} + {2.51 \over {Re \sqrt{f}}}}
$$

## The `internal-fluid-flow` Toolbox

`internal-fluid-flow` provides the following functions:

- Re2f
- f2Re
- hDeps2fDRe
- hveps2fDRe
- hQeps2fDRe
- hvthk2fDRe
- hQthk2fDRe

### `Re2f`

`Re2f` computes the Darcy friction factor *f* given the relative roughness $\varepsilon$ and the Reynolds number *Re*. If given *Re* < 2,500, then flow is assumed to be laminar and *f* is computed using of the Poiseuille condition. Otherwise, flow is assumed to be turbulent and *f* is computed using the Colebrook-White equation.

**Syntax:**

```dotnetcli
[f]=Re2f(Re,[eps[,fig]])
```

**Examples:**

Compute the Darcy friction factor f given
the Reynolds number Re = 120,000 and
the relative roughness eps = 0.001:

```dotnetcli
Re=1.2e5;eps=1e-3;
f=Re2f(Re,eps)
```

Compute f and plot a schematic Moody diagram:

```dotnetcli
f=Re2f(1.2e5,1e-3,true)
```

Compute the Darcy friction factor f given
the Reynolds number Re = 120,000
for a smooth tube and plot
a schematic Moody diagram
with the solution:

```dotnetcli
f=Re2f(1.2e5,:,true)
```

### `f2Re`

`espfD2Re` computes the Reynolds number *Re* given the relative roughness $\varepsilon$ and the Darcy friction factor *f*. Depending on the inputs, solution may be laminar or turbulent flow, or either for smooth pipes with higher friction, or none for lower friction and rough pipes. If the Poiseuille condition produces *Re* < 2,500, laminar solution is accepted. If given *f* is possible for turbulent flow,

$$
{1 \over \sqrt f} < 2 \mathrm{log} {1 \over\displaystyle {3.7 \over \varepsilon}}
$$

(which is Colebrook-White equation for for elevated *Re*) the turbulent solution is accepted. If both solutions are accepted, espfD2Re returns both answers. If neither laminar or turbulent solutions are accepted, espfD2Re returns an empty matrix. If given $\varepsilon$ > 0.05, execution is aborted.

**Syntax:**

```dotnetcli
[Re]=f2Re(f,[eps[,fig]])
```

**Examples:**

Compute the Reynolds number Re given
the Darcy friction factor f = 0.028 and
the relative roughness eps = 0.001.
In this case, both laminar and turbulent
solutions are possible:

```dotnetcli
f=2.8e-2;eps=1e-3;
Re=f2Re(f,eps)
```

Compute Re and plot a schematic Moody diagram:

```dotnetcli
Re=f2Re(2.8e-2,1e-3,true)
```

Compute the Reynolds number Re given
the Darcy friction factor f = 0.028
for a smooth tube and plot
a schematic Moody diagram
with the solution:

```dotnetcli
Re=f2Re(2.8e-2)
```

### `hDeps2fDRe`

`hDeps2fDRe` computes both the Darcy friction factor *f* and the Reynolds number *Re* given the head loss *h*, the pipe's length *L*, relative roughness $\varepsilon$ and hydraulic diameter *D*, the gravitational acceleration *g*, and the fluid's density $\rho$ and dynamic viscosity $\mu$. Replacing speed flow *v* in the Darcy-Weisbach equation by the Reynolds number *Re*,

$$
Re^2 f={2gh\rho^2D^3 \over {\mu^2 L}}
$$

Along with the Colebrook-White equation, this version of the Darcy-Weisbach equation produces a system of two equations with two variables. Solution is computed iteratively, however an analytic solution is possible in this case.

**Syntax:**

```dotnetcli
[Re,f]=hDeps2fRe(h,D,L[,eps[,rho[,mu[,g[,fig]]]]])
```

**Examples:**

Compute the Reynolds number Re and
the Darcy friction factor f, given
the head loss h = 0.40 m,
the pipe's hydraulic diameter D = 10 cm,
length L = 25 m and
relative roughness eps = 0.0027,
for water flow:

```dotnetcli
h=40;Q=1e2;L=2.5e3;eps=2.7e-3; # inputs in cgs units
[Re,f]=hDeps2fRe(h,D,L,eps)
```

Compute the Reynolds number Re and
the Darcy friction factor f, given
in addition
the fluid's density rho = 0.989 kg/L and
dynamic viscosity mu = 0.89 cP:

```dotnetcli
h=40;D=10;L=2.5e3;eps=2.7e-3;rho=0.989;mu=8.9e-3; # inputs in cgs units
[Re,f]=hDeps2fRe(h,D,L,eps,rho,mu)
```

Compute Re and f and plot a schematic Moody diagram:

```dotnetcli
[Re,f]=hDeps2fRe(0.40,0.10,25,2.7e-3,989,8.9e-4,9.81,true) # inputs in a consistent system of units
```

### `hveps2fDRe`

`hveps2fDRe` computes both the Darcy friction factor *f* and the Reynolds number *Re* given the head loss *h*, the pipe's length *L* and relative roughness $\varepsilon$, the speed flow *v*, the gravitational acceleration *g*, and the fluid's density $\rho$ and dynamic viscosity $\mu$. Replacing hydraulic diameter *D* in the Darcy-Weisbach equation by the Reynolds number *Re*,

$$
{f \over Re}={2gh\mu \over {v^3\rho L}}
$$

Along with the Colebrook-White equation, this version of the Darcy-Weisbach equation produces a system of two equations with two variables. Solution is computed iteratively.

**Syntax:**

```dotnetcli
[Re,f]=hveps2fRe(h,v,L[,eps[,rho[,mu[,g[,fig]]]]])
```

**Examples:**

Compute the Reynolds number Re and
the Darcy friction factor f, given
the head loss h = 0.40 m,
the flow speed v = 1.1 m/s,
the pipe's length L = 25 m and
relative roughness eps = 0.0027,
for water flow:

```dotnetcli
h=40;v=1.1e2;L=2.5e3;eps=2.7e-3; # inputs in cgs units
[Re,f]=hveps2fRe(h,v,L,eps)
```

Compute the Reynolds number Re and
the Darcy friction factor f, given
in addition
the fluid's density rho = 0.989 kg/L and
dynamic viscosity mu = 0.89 cP:

```dotnetcli
h=40;v=1.1e2;L=2.5e3;eps=2.7e-3;rho=0.989;mu=8.9e-3; # inputs in cgs units
[Re,f]=hveps2fRe(h,v,L,eps,rho,mu)
```

Compute Re and f and plot a schematic Moody diagram:

```dotnetcli
[Re,f]=hveps2fRe(0.40,1.1,25,2.7e-3,989,8.9e-4,9.81,true) # inputs in a consistent system of units
```

### `hQeps2fDRe`

`hQeps2fDRe` computes both the Darcy friction factor *f* and the Reynolds number *Re* given the head loss *h*, the pipe's length *L* and relative roughness $\varepsilon$, the volumetric flow rate *Q*, the gravitational acceleration *g*, and the fluid's density $\rho$ and dynamic viscosity $\mu$. Replacing hydraulic diameter *D* in the Darcy-Weisbach equation by the Reynolds number *Re*,

$$
{Re^5 f}={2ghQ^3 \over\displaystyle {{\left[ {\pi \over 4} \right]}^3 {\left[ {\mu \over \rho} \right]}^5 L}}
$$

Along with the Colebrook-White equation, this version of the Darcy-Weisbach equation produces a system of two equations with two variables. Solution is computed iteratively.

**Syntax:**

```dotnetcli
[Re,f]=hQeps2fRe(h,Q,L[,eps[,rho[,mu[,g[,fig]]]]])
```

**Examples:**

Compute the Reynolds number Re and
the Darcy friction factor f, given
the head loss h = 0.40 m,
the volumetric flow rate Q = 8.6 L/s,
the pipe's length L = 25 m and
relative roughness eps = 0.0027,
for water flow:

```dotnetcli
h=40;Q=8.6e3;L=2.5e3;eps=2.7e-3; # inputs in cgs units
[Re,f]=hQeps2fRe(h,Q,L,eps)
```

Compute the Reynolds number Re and
the Darcy friction factor f, given
in addition
the fluid's density rho = 0.989 kg/L and
dynamic viscosity mu = 0.89 cP:

```dotnetcli
h=40;Q=8.6e3;L=2.5e3;eps=2.7e-3;rho=0.989;mu=8.9e-3; # inputs in cgs units
[Re,f]=hQeps2fRe(h,Q,L,eps,rho,mu)
```

Compute Re and f and plot a schematic Moody diagram:

```dotnetcli
[Re,f]=hQeps2fRe(0.40,8.6e-3,25,2.7e-3,989,8.9e-4,9.81,true) # inputs in a consistent system of units
```

### `hvthk2fDRe`

`hvthk2fDRe` computes both the Darcy friction factor *f* and the Reynolds number *Re* given the head loss *h*, the pipe's length *L* and roughness *k*, the speed flow *v*, the gravitational acceleration *g*, and the fluid's density $\rho$ and dynamic viscosity $\mu$. Replacing hydraulic diameter *D* in the Darcy-Weisbach equation by the Reynolds number *Re*,

$$
{f \over Re}={2gh\mu \over {v^3\rho L}}
$$

Along with the Colebrook-White equation, this version of the Darcy-Weisbach equation produces a system of two equations with two variables. Solution is computed iteratively.

**Syntax:**

```dotnetcli
[Re,f]=hvthk2fRe(h,v,L[,thk[,rho[,mu[,g[,fig]]]]])
```

**Examples:**

Compute the Reynolds number Re and
the Darcy friction factor f, given
the head loss h = 0.40 m,
the flow speed v = 1.1 m/s,
the pipe's length L = 25 m and
roughness thk = 0.27 mm,
for water flow:

```dotnetcli
h=40;v=1.1e2;L=2.5e3;thk=2.7e-2; # inputs in cgs units
[Re,f]=hvthk2fRe(h,v,L,thk)
```

Compute the Reynolds number Re and
the Darcy friction factor f, given
in addition
the fluid's density rho = 0.989 kg/L and
dynamic viscosity mu = 0.89 cP:

```dotnetcli
h=40;v=1.1e2;L=2.5e3;thk=2.7e-2;rho=0.989;mu=8.9e-3; # inputs in cgs units
[Re,f]=hvthk2fRe(h,v,L,thk,rho,mu)
```

Compute Re and f and plot a schematic Moody diagram:

```dotnetcli
[Re,f]=hvthk2fRe(0.40,1.1,25,2.7e-4,989,8.9e-4,9.81,true) # inputs in a consistent system of units
```

### `hQthk2fDRe`

`hQthk2fDRe` computes both the Darcy friction factor *f* and the Reynolds number *Re* given the head loss *h*, the pipe's length *L* and roughness *k*, the volumetric flow rate *Q*, the gravitational acceleration *g*, and the fluid's density $\rho$ and dynamic viscosity $\mu$. Replacing hydraulic diameter *D* in the Darcy-Weisbach equation by the Reynolds number *Re*,

$$
{Re^5 f}={2ghQ^3 \over\displaystyle {{\left[ {\pi \over 4} \right]}^3 {\left[ {\mu \over \rho} \right]}^5 L}}
$$

Along with the Colebrook-White equation, this version of the Darcy-Weisbach equation produces a system of two equations with two variables. Solution is computed iteratively.

**Syntax:**

```dotnetcli
[Re,f]=hQthk2fRe(h,Q,L[,thk[,rho[,mu[,g[,fig]]]]])
```

**Examples:**

Compute the Reynolds number Re and
the Darcy friction factor f, given
the head loss h = 0.40 m,
the volumetric flow rate Q = 8.6 L/s,
the pipe's length L = 25 m and
roughness thk = 0.27 mm
for water flow:

```dotnetcli
h=40;Q=8.6e3;L=2.5e3;thk=2.7e-2; # inputs in cgs units
[Re,f]=hQthk2fRe(h,Q,L,thk)
```

Compute the Reynolds number Re and
the Darcy friction factor f, given
in addition
the fluid's density rho = 0.989 kg/L and
dynamic viscosity mu = 0.89 cP:

```dotnetcli
h=40;Q=8.6e3;L=2.5e3;thk=2.7e-2;rho=0.989;mu=8.9e-3; # inputs in cgs units
[Re,f]=hQthk2fRe(h,Q,L,thk,rho,mu)
```

Compute Re and f and plot a schematic Moody diagram:

```dotnetcli
[Re,f]=hQthk2fRe(0.40,8.6e-3,25,2.7e-4,989,8.9e-4,9.81,true) # inputs in a consistent system of units
```

Copyright &copy; 2022 Alexandre Umpierre

email: aumpierre@gmail.com
