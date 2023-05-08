# Copyright (C) 2022 2023 Alexandre Umpierre
#
# This file is part of internal-fluid-flow toolbox for GNU Octave.
# internal-fluid-flow toolbox for GNU Octave is free software:
# you can redistribute it and/or modify it under the terms
# of the GNU General Public License (GPL) version 3
# as published by the Free Software Foundation.
#
# internal-fluid-flow toolbox for GNU Octave is distributed in the hope
# that it will be useful,but WITHOUT ANY WARRANTY;
# without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the
# GNU General Public License along with this program
# (license GNU GPLv3.txt).
# It is also available at https://www.gnu.org/licenses/.

function [Re,f]=h2fRe(h,D=NaN,v=NaN,Q=NaN,eps=NaN,k=NaN,L=100,rho=0.997,mu=9.1e-3,g=981,fig=false)
    # Syntax:
    # -- [Re,f]=h2fRe(h,D,:,:,eps,:[,L][,rho][,mu][,g][,fig])
    # -- [Re,f]=h2fRe(h,:,v,:,eps,:[,L][,rho][,mu][,g][,fig])
    # -- [Re,f]=h2fRe(h,:,:,Q,eps,:[,L][,rho][,mu][,g][,fig])
    # -- [Re,f]=h2fRe(h,D,:,:,:,k[,L][,rho][,mu][,g][,fig])
    # -- [Re,f]=h2fRe(h,:,v,:,:,k[,L][,rho][,mu][,g][,fig])
    # -- [Re,f]=h2fRe(h,:,:,Q,:,k[,L][,rho][,mu][,g][,fig])
    #
    # h2fRe computes the Reynolds number Re and
    #  the Darcy friction factor f given
    #  the head loss h,
    #  the pipe's hydraulic diameter D or
    #  the flow speed v or
    #  the volumetric flow rate Q,
    #  the pipe's length L (default L = 100),
    #  the pipe's roughness k (default k = 0) or
    #  the pipe's relative roughness eps (default eps = 0),
    #  the fluid's density rho (default rho = 0.997),
    #  the fluid's dynamic viscosity mu (default mu = 0.0091),and
    #  the gravitational accelaration g (default g = 981).
    # By default,pipe is assumed to be smooth.
    #  Relative roughness is reset to eps = 0.05,if eps > 0.05.
    # Notice that default values are given in the cgs unit system and,
    #  if taken,all other parameters must as well
    #  be given in cgs units.
    # If parameter fig = true is given
    #  a schematic Moody diagram
    #  is plotted as a graphical representation
    #  of the solution.
    # h2fRe is a main function of
    #  the internal-fluid-flow toolbox for GNU Octave.
    #
    # Examples:
    # # Compute the Reynolds number and
    # # the Darcy friction factor given
    # # the head loss is 40 cm,
    # # the pipe's hydraulic diameter is 10 cm,
    # # the pipe's length is 25 m and
    # # the pipe's relative roughness is 0.0027
    # # for water flow:
    # [Re,f]=h2fRe(h=40,D=10,:,:,eps=2.7e-3,:,L=2.5e3)
    #
    # # Compute the Reynolds number and
    # # the Darcy friction factor given
    # # the head loss per meter is 1.6 cm/m,
    # # the volumetric flow rate is 8.6 L/s,
    # # the fluid's density is 0.989 g/cc and
    # # the fluid's dynamic viscosity is 0.89 cP
    # # for a smooth pipe and
    # # show results on a schematic Moody diagram:
    # [Re,f]=h2fRe(h=1.6,:,:,Q=8.6e3,eps=0,:,L=1,rho=0.989,mu=8.9e-3,:,true)
    #
    # # Compute the Reynolds number and
    # # the Darcy friction factor given
    # # the head loss is 0.40 m,
    # # the flow speed is 1.1 m/s,
    # # the pipe's length is 25 m
    # # for water flow for a smooth pipe:
    # [Re,f]=h2fRe(h=40,:,v=1.1e2,:,:,k=0,L=2.5e3)
    #
    # See also: Re2f,f2Re.
    a=isnan([D,v,Q])!=1;
    if sum(a)!=1
        error("h2fRe requires that either\nthe hydraulic diameter,\nthe flow speed or\nthe flow rate\nbe given alone.")
    end
    b=isnan([eps,k])!=1;
    if sum(b)!=1
        error("h2fRe requires that either\nthe pipe's roughness or\nthe pipe's relative roughness\nbe given alone.")
    end
    if a==[1,0,0] && b==[1,0]
        [Re,f]=hDeps2fRe(h,D,L,eps,rho,mu,g,fig);
    elseif a==[1,0,0] && b==[0,1]
        [Re,f]=hDeps2fRe(h,D,L,k / D,rho,mu,g,fig);
    elseif a==[0,1,0] && b==[1,0]
        [Re,f]=hveps2fRe(h,v,L,eps,rho,mu,g,fig);
    elseif a==[0,1,0] && b==[0,1]
        [Re,f]=hvthk2fRe(h,v,L,k,rho,mu,g,fig);
    elseif a==[0,0,1] && b==[1,0]
        [Re,f]=hQeps2fRe(h,Q,L,eps,rho,mu,g,fig);
    elseif a==[0,0,1] && b==[0,1]
        [Re,f]=hQthk2fRe(h,Q,L,k,rho,mu,g,fig);
    end
end

