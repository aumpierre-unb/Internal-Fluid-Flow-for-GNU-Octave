## Using the InternalFluidFlow module for Julia

InernalFluidFlow can be added and loaded either 
from the JuliaHub repository or from the source repository.

<hr/>

### Add and load InternalFluidFlow from JuliaHub repository (last released version)

``julia> using Pkg``

``julia> Pkg.add("InternalFluidFlow")``

``julia> using InternalFluidFlow``

### Add and load InternalFluidFlow from source repository (pre-realease version)

``(@v1.8) pkg> add "https://github.com/aumpierre-unb/InternalFluidFlow.jl"``

### The InternalFluidFlow

InternalFluidFlow module has seven functions. 

All functions have an optional input argument with 
default value *false*. If it is given as *true*, 
a schematic Moody diagram is plot with the computed solution.

The first two are basic functions, to compute either 
the Reynolds number *Re* or the Darcy friction factor *f*, 
provided the other is given along with the relative roughness *eps*:

``julia> f=Re2f(Re,[eps[,fig]])``

``julia> Re=f2Re(f,[eps[,fig]])``

The other five functions compute both 
the Reynolds number *Re* and the Darcy friction factor *f*, given 
the head loss *h*, 
the pipe's length *L*, 
the gravitational acceleration *g*, 
the fluid's density *rho* and dynamic viscosity *mu*, 
the pipe's roughness *thk* or relative roughness *eps*, and 
the pipe's hydraulic diameter *D* or the speed flow *v* or the volumetric rate flow *Q*:

``julia> Re,f=hDeps2fRe(h,D,L,eps,g,mu,rho[,fig])``

``julia> Re,f=hveps2fRe(h,v,L,eps,g,mu,rho[,fig])``

``julia> Re,f=hvthk2fRe(h,v,L,thk,g,mu,rho[,fig])``

``julia> Re,f=hQeps2fRe(h,Q,L,eps,g,mu,rho[,fig])``

``julia> Re,f=hQthk2fRe(h,Q,L,thk,g,mu,rho[,fig])``

For a short description of InternalFluidFlow module, type

``help?> InternalFluidFlow``

For further details on syntax, check on the functions descriptions. *e.g.*, for the description of function hQthk2fRe, type

``help?> hQthk2fRe``