# Internal Fluid Flow Toolbox for GNU-Octave

[![DOI](https://zenodo.org/badge/509427410.svg)](https://zenodo.org/badge/latestdoi/509427410)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Install and load Internal Fluid Flow

```dotnetcli
>> pkg install https://github.com/aumpierre-unb/Internal-Fluid-Flow-for-GNU-Octave/archive/refs/tags/v0.1.2.tar.gz``

>> pkg load internal-fluid-flow
```

---

The following is a very short introduction to the steady internal flow of an incompressible and inviscid fluid and to the Internal Fluid Flow Toolbox for GNU Octave.

Our focus here is a small set of equations that described the phenomenon and are required to solve problems on internal fluid flow. Internal flow is a pretty extensive topic in fluid mechanics and there are a lot of important and interesting observations related to it that are not taken into account in this text, because they have no direct impact the computation performed by the functions in this package.

This text is divided in two main sections: The Theory and The Internal Fluid Flow Toolbox.

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

For laminar flow, *Re* < 2500 (typically), the Darcy friction factor is given by the Poiseuille condition,

$$
f={64 \over Re}
$$

For turbulent flow, *Re* > 2500 (typically), the Darcy friction factor is given implicitly by the Colebrook-White equation,

$$
{1 \over \sqrt{f}}=2 \mathrm{log} {1 \over\displaystyle {3.7 \over \varepsilon} + {2.51 \over {Re \sqrt{f}}}}
$$

## The Internal Fluid Flow Toolbox

<!-- This package provides a set of functions designed to solve problems of internal fluid flow. All functions are based on the Poiseuille condition for laminar flow, the Colebrook-White equation for turbulent flow, and the Darcy-Weisbach equation for head loss. The simplest problems on internal flow consist in computing either the Reynolds number or the Darcy friction factor given the other and the relative roughness. For those cases, this package provides functions f2Re and Re2f, respectively. More elaborated problems consist in computing both the Reynolds number and the Darcy friction factor given the head loss, the tube length, the fluid's density and dynamic viscosity, the gravitational acceleration, the relative roughness and either the dynamic diameter or the linear velocity or the volumetric flow. For those cases, this package provides functions hDeps2fRe, hveps2fRe and hQeps2fRe, respectively. A slightly more elaborate situation arises when roughness is given instead of relative roughness along with the linear velocity or the volumetric flow. For those cases, this package provides functions hvthk2fRe and hQthk2fRe, respectively. All function in this package offer the option of plotting the solution on a schematic Moody diagram. -->

Internal Fluid Flow provides the following functions:

- Re2f
- f2Re
- hDeps2fDRe
- hveps2fDRe
- hQeps2fDRe
- hvthk2fDRe
- hQthk2fDRe

### Re2f

Re2f computes the Darcy friction factor *f* given the relative roughness $\varepsilon$ and the Reynolds number *Re*. If given *Re* < 2500, then flow is assumed to be laminar and *f* is computed using of the Poiseuille condition. Otherwise, flow is assumed to be turbulent and *f* is computed using the Colebrook-White equation.

**Syntax:**

``[f]=Re2f(Re,[eps[,fig]])``

*e.g.* this call computes *f* and shows no plot:

``>> Re=1.2e5;eps=0.002;``

``>> f=Re2f(Re,eps)``

*e.g.* this call computes *f* for the default relative roughness $\varepsilon$ = 0.002 and shows plot:

``>> f=Re2f(1.2e5,:,true)``

*e.g.* this call computes *f* and shows plot:

``>> f=Re2f(1.2e5,0.002,true)``

### f2Re

espfD2Re computes the Reynolds number *Re* given the relative roughness $\varepsilon$ and the Darcy friction factor *f*. Depending on the inputs, solution may be laminar or turbulent flow, or either for smooth pipes with higher friction, or none for lower friction and rough pipes. If the Poiseuille condition produces Re < 2500, laminar solution is accepted. If given *f* is possible for turbulent flow,

$$
{1 \over \sqrt f} < 2 \mathrm{log} {1 \over\displaystyle {3.7 \over \varepsilon}}
$$

(which is Colebrook-White equation for for elevated *Re*) the turbulent solution is accepted. If both solutions are accepted, espfD2Re returns both answers. If neither laminar or turbulent solutions are accepted, espfD2Re returns an empty matrix. If given $\varepsilon$ > 0.05, execution is aborted.

**Syntax:**

``[Re]=f2Re(f,[eps[,fig]])``

*e.g.* this call computes *Re* for both laminar and turbulent regimes (if possible) and shows no plot:

``>> f=0.025;eps=0.002;``

``>> Re=f2Re(f,eps)``

*e.g.* this call computes *Re* for the default relative roughness $\varepsilon$ = 0.002 for both laminar and turbulent regimes (if possible) and shows plot:

``>> Re=f2Re(0.025,:,true)``

*e.g.* e.g. this call computes *Re* for both laminar and turbulent regimes (if possible) and shows plot:

``>> Re=f2Re(0.025,0.002,true)``

### hDeps2fDRe

hDeps2fDRe computes both the Darcy friction factor *f* and the Reynolds number *Re* given the head loss *h*, the pipe's length *L*, relative roughness $\varepsilon$ and hydraulic diameter *D*, the gravitational acceleration *g*, and the fluid's density $\rho$ and dynamic viscosity $\mu$. Replacing speed flow *v* in the Darcy-Weisbach equation by the Reynolds number *Re*,

$$
Re^2 f={2gh\rho^2D^3 \over {\mu^2 L}}
$$

Along with the Colebrook-White equation, this version of the Darcy-Weisbach equation produces a system of two equations with two variables. Solution is computed iteratively, however an analytic solution is possible in this case.

**Syntax:**

``[Re,f]=hDeps2fRe(h,D,L[,eps[,rho[,mu[,g[,fig]]]]])``

*e.g.* this call computes *Re* and *f* and shows no plot:

``>> h=40;D=10;L=2500;eps=0.0025;g=981;mu=0.0089;rho=0.989;``

``>> thk=eps*D``

``>> [Re,f]=hDeps2fRe(h,D,L,eps,rho,mu,g)``

``>> v=Re*mu/rho/D``

``>> Q=v*(pi/4*D^2)``

*e.g.* this call computes *Re* and *f* and shows plot:

``>> [Re,f]=hDeps2fRe(40,10,2500,0.0025,0.989,0.0089,981,true)``

### hveps2fDRe

hveps2fDRe computes both the Darcy friction factor *f* and the Reynolds number *Re* given the head loss *h*, the pipe's length *L* and relative roughness $\varepsilon$, the speed flow *v*, the gravitational acceleration *g*, and the fluid's density $\rho$ and dynamic viscosity $\mu$. Replacing hydraulic diameter *D* in the Darcy-Weisbach equation by the Reynolds number *Re*,

$$
{f \over Re}={2gh\mu \over {v^3\rho L}}
$$

Along with the Colebrook-White equation, this version of the Darcy-Weisbach equation produces a system of two equations with two variables. Solution is computed iteratively.

**Syntax:**

``[Re,f]=hveps2fRe(h,v,L[,eps[,rho[,mu[,g[,fig]]]]])``

*e.g.* this call computes *Re* and *f* and shows no plot:

``>> h=40;v=110;L=2500;eps=0.0025;g=981;mu=0.0089;rho=0.989;``

``>> [Re,f]=hveps2fRe(h,v,L,eps,rho,mu,g)``

``>> D=Re*mu/rho/v``

``>> thk=eps*D``

``>> Q=v*(pi/4*D^2)``

*e.g.* this call computes *Re* and *f* and shows plot:

``>> [Re,f]=hveps2fRe(40,110,2500,0.0025,0.989,0.0089,981,true)``

### hQeps2fDRe

hQeps2fDRe computes both the Darcy friction factor *f* and the Reynolds number *Re* given the head loss *h*, the pipe's length *L* and relative roughness $\varepsilon$, the volumetric flow rate *Q*, the gravitational acceleration *g*, and the fluid's density $\rho$ and dynamic viscosity $\mu$. Replacing hydraulic diameter *D* in the Darcy-Weisbach equation by the Reynolds number *Re*,

$$
{Re^5 f}={2ghQ^3 \over\displaystyle {{\left[ {\pi \over 4} \right]}^3 {\left[ {\mu \over \rho} \right]}^5 L}}
$$

Along with the Colebrook-White equation, this version of the Darcy-Weisbach equation produces a system of two equations with two variables. Solution is computed iteratively.

**Syntax:**

``[Re,f]=hQeps2fRe(h,Q,L[,eps[,rho[,mu[,g[,fig]]]]])``

*e.g.* this call computes *Re* and *f* and shows no plot:

``>> h=40;Q=8600;L=2500;eps=0.0025;g=981;mu=0.0089;rho=0.989;``

``>> [Re,f]=hQeps2fRe(h,Q,L,eps,rho,mu,g)``

``>> D=Q*rho/(pi/4)/Re/mu``

``>> thk=eps*D``

``>> v=Q/(pi/4*D^2)``

*e.g.* this call computes *Re* and *f* and shows plot:

``>> [Re,f]=hQeps2fRe(40,8600,2500,0.0025,0.989,0.0089,981,true)``

### hvthk2fDRe

hvthk2fDRe computes both the Darcy friction factor *f* and the Reynolds number *Re* given the head loss *h*, the pipe's length *L* and roughness *k*, the speed flow *v*, the gravitational acceleration *g*, and the fluid's density $\rho$ and dynamic viscosity $\mu$. Replacing hydraulic diameter *D* in the Darcy-Weisbach equation by the Reynolds number *Re*,

$$
{f \over Re}={2gh\mu \over {v^3\rho L}}
$$

Along with the Colebrook-White equation, this version of the Darcy-Weisbach equation produces a system of two equations with two variables. Solution is computed iteratively.

**Syntax:**

``[Re,f]=hvthk2fRe(h,v,L[,thk[,rho[,mu[,g[,fig]]]]])``

*e.g.* this call computes *Re* and *f* and shows no plot:

``>> h=40;v=110;L=2500;thk=0.025;g=981;mu=0.0089;rho=0.989;``

``>> [Re,f]=hvthk2fRe(h,v,L,thk,rho,mu,g)``

``>> D=Re*mu/rho/v``

``>> eps=thk/D``

``>> Q=v*(pi/4*D^2)``

*e.g.* this call computes *Re* and *f* and shows plot:

``>> [Re,f]=hvthk2fRe(40,110,2500,0.025,0.989,0.0089,981,true)``

### hQthk2fDRe

hQthk2fDRe computes both the Darcy friction factor *f* and the Reynolds number *Re* given the head loss *h*, the pipe's length *L* and roughness *k*, the volumetric flow rate *Q*, the gravitational acceleration *g*, and the fluid's density $\rho$ and dynamic viscosity $\mu$. Replacing hydraulic diameter *D* in the Darcy-Weisbach equation by the Reynolds number *Re*,

$$
{Re^5 f}={2ghQ^3 \over\displaystyle {{\left[ {\pi \over 4} \right]}^3 {\left[ {\mu \over \rho} \right]}^5 L}}
$$

Along with the Colebrook-White equation, this version of the Darcy-Weisbach equation produces a system of two equations with two variables. Solution is computed iteratively.

**Syntax:**

``[Re,f]=hQthk2fRe(h,Q,L[,thk[,rho[,mu[,g[,fig]]]]])``

*e.g.* this call computes *Re* and *f* and shows no plot:

``>> h=40;Q=8600;L=2500;thk=0.025;g=981;mu=0.0089;rho=0.989;``

``>> [Re,f]=hQthk2fRe(h,Q,L,thk,rho,mu,g)``

``>> D=Q*rho/(pi/4)/Re/mu``

``>> eps=thk/D``

``>> v=Q/(pi/4*D^2)``

*e.g.* this call computes *Re* and *f* and shows plot:

``>> [Re,f]=hQthk2fRe(40,8600,2500,0.025,0.989,0.0089,981,true)``

Copyright &copy; 2022 Alexandre Umpierre

email: aumpierre@gmail.com
