# Internal-Fluid-Flow-for-GNU-Octave
This package provides a set of functions designed 
to solve problems of internal fluid flow. 
All functions are based on the Colebrooke-White 
equation for turbulent flow, the Poiseuille 
condition for laminar flow, and the Darcy-Weissbach 
equation for head loss. The simplest problems 
on internal flow consist in computing either 
the Reynolds number or the Darcy friction factor 
given the other and the relative roughness. 
For those cases, this package provides functions 
f2Re and Re2f, respectively. More elaborated 
problems consist in computing both the Reynolds 
number and the Darcy friction factor given the 
head loss, the tube length, the fluid's density and 
dynamic viscosity, the gravitational acceleration, 
the relative roughness and either the dynamic diameter 
or the linear velocity or the volumetric flow. 
For those cases, this package provides functions 
hDeps2fRe, hveps2fRe and hQeps2fRe, respectively. 
A slightly more elaborate situation arises when 
roughness is given instead of relative roughness 
along with the linear velocity or the volumetric flow. 
For those cases, this package provides functions 
hvthk2fRe and hQthk2fRe, respectively. All function 
in this package offer the option of plotting 
the solution on a schematic Moody diagram. 
Enjoy!
