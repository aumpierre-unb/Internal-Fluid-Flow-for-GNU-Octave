# `internal-fluid-flow` Toolbox for GNU-Octave

[![DOI](https://zenodo.org/badge/509427410.svg)](https://zenodo.org/badge/latestdoi/509427410)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/aumpierre-unb/Internal-Fluid-Flow-for-GNU-Octave)

![Illustrative graphical output](https://github.com/aumpierre-unb/Internal-Fluid-Flow-for-GNU-Octave/blob/main/pics/D2fRe.png "Example of graphical output")

## Installing and Loading `internal-fluid-flow`

```dotnetcli
# e.g. this call installs version 0.3.2
pkg install https://github.com/aumpierre-unb/Internal-Fluid-Flow-for-GNU-Octave/archive/refs/tags/v0.3.2.tar.gz
pkg load internal-fluid-flow
```

## Citation of `internal-fluid-flow`

You can cite all versions (both released and pre-released), by using
[DOI 105281/zenodo.6960263](https://doi.org/10.5281/zenodo.6960263).
This DOI represents all versions, and will always resolve to the latest one.

<!--To cite the last released version, please check
https://zenodo.org/account/settings/github/repository/aumpierre-unb/Internal-Fluid-Flow-for-GNU-Octave.-->

---

The following is a very short introduction to the steady internal flow of an incompressible and inviscid fluid and to the `internal-fluid-flow` toolbox for GNU Octave.

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
{1 \over \sqrt{f}}=2 \mathrm{log} {1 \over\displaystyle {\varepsilon \over 3.7} + {2.51 \over {Re \sqrt{f}}}}
$$

## The `internal-fluid-flow` Toolbox

`internal-fluid-flow` provides the following functions:

- `Re2f`
- `f2Re`
- `h2fDRe`

### `Re2f`

`Re2f` computes the Darcy friction factor *f* given the relative roughness $\varepsilon$ and the Reynolds number *Re*. If given *Re* < 2,500, then flow is assumed to be laminar and *f* is computed using of the Poiseuille condition. Otherwise, flow is assumed to be turbulent and *f* is computed using the Colebrook-White equation.

**Syntax:**

```dotnetcli
 -- f=Re2f(Re[,eps][,fig])
```

By default, pipe is assumed to be smooth, eps = 0. If eps > 0.05, eps is reset to eps = 0.05.

If fig = true is given, a schematic Moody diagram is plotted as a graphical representation of the solution.

**Examples:**

Compute the Darcy friction factor *f* given the Reynolds number *Re* = 120,000 and the relative roughness $\varepsilon$ = 0.001:

```dotnetcli
f=Re2f(120e3,1e-3)
```

Compute the Darcy friction factor *f* given the Reynolds number *Re* = 120,000 for a smooth tube and displays a schematic Moody diagram:

```dotnetcli
f=Re2f(120e3,:,true)
```

### `f2Re`

`f2Re` computes the Reynolds number *Re* given the relative roughness $\varepsilon$ and the Darcy friction factor *f*. Depending on the inputs, solution may be laminar or turbulent flow, or either for smooth pipes with higher friction, or none for lower friction and rough pipes. If the Poiseuille condition produces *Re* < 2,500, laminar solution is accepted. If given *f* is possible for turbulent flow,

$$
{1 \over \sqrt f} < 2 \mathrm{log} {1 \over\displaystyle {3.7 \over \varepsilon}}
$$

(which is Colebrook-White equation for for elevated *Re*) the turbulent solution is accepted. If both solutions are accepted, `f2Re` returns both answers. If neither laminar or turbulent solutions are accepted, `f2Re` returns an empty matrix.

**Syntax:**

```dotnetcli
 -- Re=f2Re(f[,eps][,fig])
```

By default, pipe is assumed to be smooth, eps = 0. If eps > 0.05, eps is reset to eps = 0.05.

If fig = true is given, a schematic Moody diagram is plotted as a graphical representation of the solution.

**Examples:**

Compute the Reynolds number *Re* given the Darcy friction factor *f* = 0.028 and the relative roughness $\varepsilon$ = 0.001. In this case, both laminar and turbulent solutions are possible:

```dotnetcli
Re=f2Re(2.8e-2,1e-3)
```

Compute the Reynolds number Re given the Darcy friction factor f = 0.028 for a smooth pipe and displays a schematic Moody diagram. In this case, both turbulent and laminar solutions are possible:

```dotnetcli
Re=f2Re(2.8e-2,:,true)
```

### `h2fDRe`

`h2fDRe` computes both the Darcy friction factor *f* and the Reynolds number *Re* given the head loss *h*, the pipe's hydraulic diameter *D* or the flow speed *v* or the volumetric flow rate *Q*, the pipe's length *L*, the pipe's roughness *k* or the pipe's relative roughness $\varepsilon$, the fluid's density $\rho$, the fluid's dynamic viscosity $\mu$, and the gravitational accelaration *g*.

**Syntax:**

```dotnetcli
 -- [Re,f]=h2fRe(h,D,:,:,eps,:[,L][,rho][,mu][,g][,fig])
 -- [Re,f]=h2fRe(h,:,v,:,eps,:[,L][,rho][,mu][,g][,fig])
 -- [Re,f]=h2fRe(h,:,:,Q,eps,:[,L][,rho][,mu][,g][,fig])
 -- [Re,f]=h2fRe(h,D,:,:,:,k[,L][,rho][,mu][,g][,fig])
 -- [Re,f]=h2fRe(h,:,v,:,:,k[,L][,rho][,mu][,g][,fig])
 -- [Re,f]=h2fRe(h,:,:,Q,:,k[,L][,rho][,mu][,g][,fig])
```

By default, pipe is assumed to be 1 m long, L = 100 (in cm).

By default, pipe is assumed to be smooth. Relative roughness is reset to eps = 0.05, if eps > 0.05.

By default, fluid is assumed to be water at 25 °C, and
rho = 0.997 (in g/cc) and mu = 0.0091 (in g/cm/s).

By default, gravitational acceleration is that of Earth,
g = 981 (in cm/s/s).

Notice that default values are given in the cgs unit system and, if taken, all other parameters must as well be given in cgs units.

If parameter fig = true is given a schematic Moody diagram is plotted as a graphical representation of the solution.

**Examples:**

Compute the Reynolds number *Re* and
the Darcy friction factor *f* given
the head loss *h* = 0.40 m,
the pipe's hydraulic diameter *D* = 10 cm,
length *L* = 25 m and
relative roughness $\varepsilon$ = 0.0027,
for water flow:

```dotnetcli
[Re,f]=h2fRe(40,10,:,:,2.7e-3,:,2.5e3)
```

Compute the Reynolds number *Re* and the Darcy friction factor *f* given the head loss per meter *h*/*L* = 1.6 cm/m, the volumetric flow rate *Q* = 8.6 L/s, the fluid's density $\rho$ = 0.989 g/cc and dynamic viscosity $\mu$ = 0.89 cP for a smooth pipe and show results on a schematic Moody diagram:

```dotnetcli
[Re,f]=h2fRe(1.6,:,:,8.6e3,0,:,1,0.989,8.9e-3,:,true)
```

Compute the Reynolds number *Re* and the Darcy friction factor *f*, given the head loss *h* = 0.40 m, the flow speed *v* = 1.1 m/s, the pipe's length *L* = 25 m for water flow for a smooth pipe:

```dotnetcli
[Re,f]=h2fRe(40,:,1.1e2,:,:,0,2.5e3)
```

### See Also

[Psychrometrics-for-GNU-Octave](https://github.com/aumpierre-unb/Psychrometrics-for-GNU-Octave),
[McCabe-Thiele-for-GNU-Octave](https://github.com/aumpierre-unb/McCabe-Thiele-for-GNU-Octave),
[Ponchon-Savarit-for-GNU-Octave](https://github.com/aumpierre-unb/Ponchon-Savarit-for-GNU-Octave).

Copyright &copy; 2022 2023 Alexandre Umpierre

email: aumpierre@gmail.com
