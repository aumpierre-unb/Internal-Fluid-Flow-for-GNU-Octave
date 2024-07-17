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

## The `internal-fluid-flow` Toolbox

`internal-fluid-flow` provides the following functions:

- `Re2f`
- `f2Re`
- `h2fDRe`

### `Re2f`

`Re2f` computes the Darcy friction factor f given the relative roughness eps and the Reynolds number Re. If given Re < 2,500, then flow is assumed to be laminar and f is computed using of the Poiseuille condition. Otherwise, flow is assumed to be turbulent and f is computed using the Colebrook-White equation.

**Syntax:**

```dotnetcli
 -- f=Re2f(Re[,eps][,fig])
```

By default, pipe is assumed to be smooth, eps = 0. If eps > 0.05, eps is reset to eps = 0.05.

If fig = true is given, a schematic Moody diagram is plotted as a graphical representation of the solution.

**Examples:**

Compute the Darcy friction factor given the Reynolds number is 120,000 and the relative roughness is 0.001:

```dotnetcli
f=Re2f(Re=120e3,eps=1e-3)
```

Compute the Darcy friction factor given the Reynolds number is 120,000 for a smooth tube and displays a schematic Moody diagram:

```dotnetcli
f=Re2f(Re=120e3,:,true)
```

### `f2Re`

`f2Re` computes the Reynolds number Re given the relative roughness eps and the Darcy friction factor f. Depending on the inputs, solution may be laminar or turbulent flow, or either for smooth pipes with higher friction, or none for lower friction and rough pipes. If the Poiseuille condition produces Re < 2,500, laminar solution is accepted. If given f is possible for turbulent flow,

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

Compute the Reynolds number given the Darcy friction factor is 0.028 and the relative roughness is 0.001. In this case, both laminar and turbulent solutions are possible:

```dotnetcli
Re=f2Re(f=2.8e-2,eps=1e-3)
```

Compute the Reynolds number Re given the Darcy friction factor is 0.028 for a smooth pipe and displays a schematic Moody diagram. In this case, both turbulent and laminar solutions are possible:

```dotnetcli
Re=f2Re(f=2.8e-2,:,true)
```

### `h2fDRe`

`h2fDRe` computes both the Darcy friction factor f and the Reynolds number Re given the head loss h, the pipe's hydraulic diameter D or the flow speed v or the volumetric flow rate Q, the pipe's length L, the pipe's roughness k or the pipe's relative roughness eps, the fluid's density rho, the fluid's dynamic viscosity mu, and the gravitational accelaration g.

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

Compute the Reynolds number and the Darcy friction factor given the head loss is 0.40 m, the pipe's hydraulic diameter is 10 cm, the pipe's length is 25 m and the pipe's relative roughness is 0.0027, for water flow:

```dotnetcli
[Re,f]=h2fRe(h=40,D=10,:,:,eps=2.7e-3,:,L=2.5e3)
```

Compute the Reynolds number and the Darcy friction factor given the head loss per meter is 1.6 cm/m, the volumetric flow rate is 8.6 L/s, the fluid's density is 0.989 g/cc and the fluid's dynamic viscosity is 0.89 cP for a smooth pipe and show results on a schematic Moody diagram:

```dotnetcli
[Re,f]=h2fRe(h=1.6,:,:,Q=8.6e3,eps=0,:,L=1,rho=0.989,mu=8.9e-3,:,true)
```

Compute the Reynolds number and the Darcy friction factor given the head loss is 0.40 m, the flow speed is 1.1 m/s, the pipe's length is 25 m for water flow for a smooth pipe:

```dotnetcli
[Re,f]=h2fRe(h=40,:,v=1.1e2,:,:,k=0,L=2.5e3)
```

### See Also

[Psychrometrics-for-GNU-Octave](https://github.com/aumpierre-unb/Psychrometrics-for-GNU-Octave),
[McCabe-Thiele-for-GNU-Octave](https://github.com/aumpierre-unb/McCabe-Thiele-for-GNU-Octave),
[Ponchon-Savarit-for-GNU-Octave](https://github.com/aumpierre-unb/Ponchon-Savarit-for-GNU-Octave).

Copyright &copy; 2022 2023 Alexandre Umpierre

email: aumpierre@gmail.com
