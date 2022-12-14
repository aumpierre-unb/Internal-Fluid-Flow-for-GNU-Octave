# Copyright (C) 2022 Alexandre Umpierre
#
# This file is part of internal-fluid-flow toolbox for GNU Octave.
# internal-fluid-flow toolbox for GNU Octave is free software:
# you can redistribute it and/or modify it under the terms
# of the GNU General Public License (GPL) version 3
# as published by the Free Software Foundation.
#
# internal-fluid-flow toolbox for GNU Octave is distributed in the hope
# that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the
# GNU General Public License along with this program
# (license GNU GPLv3.txt).
# It is also available at https://www.gnu.org/licenses/.

function [Re,f]=hQthk2fRe(h,Q,L,thk=0,rho=0.997,mu=9.1e-3,g=981,fig=false)
    # Syntax:
    # [Re,f]=hQthk2fRe(h,Q,L,thk,g,mu,rho[,fig=false])
    #
    # hQthk2fRe computes
    #  the Reynolds number Re and
    #  the Darcy friction factor f given
    #  the head loss h,
    #  the volumetric flow Q,
    #  the pipe"s length L,
    #  the pipe"s roughness thk,
    #  the fluid"s density rho,
    #  the fluid"s dynamic viscosity mu, and
    #  the gravitational accelaration g.
    # By default, pipe is assumed to be smooth, thk = 0.
    # By default, the fluid is assumed to be water at 25 degC,
    #  rho = 0.997 kg/L and mu = 0.91 cP,
    #  and gravitational acceleration is assumed to be
    #  g = 9.81 m/s/s.
    # Please, notice that these default values are given
    #  in the cgs unit system and, if taken,
    #  all other inputs must as well be given in cgs units.
    # If fig = true is given, a schematic Moody diagram
    #  is plotted as a graphical representation
    #  of the solution.
    # hQthk2fRe is a main function of
    #  the internal-fluid-flow toolbox for GNU Octave.
    #
    # Examples:
    # # Compute the Reynolds number Re and
    # # the Darcy friction factor f given
    # # the head loss h = 0.40 m,
    # # the volumetric flow rate Q = 8.6 L/s,
    # # length L = 25 m and
    # # roughness thk = 0.27 mm,
    # # for water:
    # h=40;Q=8.6e3;L=2.5e3;thk=2.7e-2; # inputs in cgs units
    # [Re,f]=hQthk2fRe(h,Q,L,thk)
    #
    # # Compute the Reynolds number Re and
    # # the Darcy friction factor f given
    # # the head loss h = 0.40 m,
    # # the volumetric flow rate Q = 8.6 L/s,
    # # length L = 25 m and
    # # the fluid"s density rho = 0.989 kg/L and
    # # dynamic viscosity mu = 0.89 cP, and
    # # in a smooth pipe:
    # h=40;Q=8.6e3;L=2.5e3;rho=0.989;mu=8.9e-3; # inputs in cgs units
    # [Re,f]=hQthk2fRe(h,Q,L,:,rho,mu)
    #
    # # Compute the Reynolds number Re and
    # # the Darcy friction factor f given
    # # the head loss h = 0.40 m,
    # # the volumetric flow rate Q = 8.6 L/s,
    # # length L = 25 m and
    # # roughness thk = 0.27 mm,
    # # the fluid"s dynamic viscosity mu = 0.89 cP and
    # # density rho = 0.989 kg/L, and
    # # display a schematic Moody Diagram:
    # [Re,f]=hQthk2fRe(0.40,8.6e-3,25,2.7e-4,989,8.9e-4,9.81,true)
    #
    # See also: Re2f, f2Re, hDeps2fRe, hveps2fRe, hvthk2fRe, hQeps2fRe.
    P=2*g*h*Q^3/(pi/4)^3/(mu/rho)^5/L;
    foo=@(f) (1/f^(1/2)+2*log10(...
    thk/(rho*Q/(pi/4)/mu/((P/f)^(1/5)))... # this is eps
    /3.7+2.51/(P/f)^(1/5)/f^(1/2)));
    f=newtonraphson(foo,1e-2,1e-4);
    Re=(P/f)^(1/5);
    if Re>2.3e3
        islam=false;
    else
        Re=(P/64)^(1/4);
        f=64/Re;
        islam=true;
    end
    D=rho*Q/Re/mu/(pi/4);
    eps=thk/D;
    if fig
        figure;
        hold on;
        x=[5e-2 2.5e-2 1e-2 3e-3 1e-3 3e-4 1e-4];
        for i=1:length(x)
            if abs(x(i)-eps)>eps/10
                turbulent(x(i),"k",1);
                feps=(-2*log10(x(i)/3.7))^-2;
                text(8e7,feps*1.07,num2str(x(i),4),"color","k","fontsize",11,"horizontalalignment","right");
            end
        end
        rough("-.b",1.5);
        if thk~=0
            smooth("-.b",1.5);
            text(7e6,8e-3,"Smooth pipe","color","b","fontsize",11,"horizontalalignment","right");
            text(4e4,7.6e-2,"Fully rough flow","color","b","fontsize",11);
        else
            text(7e6,8e-3,"Smooth pipe","color","r","fontsize",11,"horizontalalignment","right");
            text(4e4,7.5e-2,"Fully rough flow","color","b","fontsize",11);
        end
        if islam
            laminar("r",2);
        else
            laminar("k",1);
            turbulent(eps,"r",2);
            feps=(-2*log10(eps/3.7))^-2;
            text(8e6,feps*1.07,num2str(eps,4),"color","r","fontsize",11,"horizontalalignment","right");
        end
        loglog(Re,f,"or","markersize",8,"markerfacecolor","r");
        line("xdata",[(P/6e-3)^(1/5) (P/1e-1)^(1/5)],...
             "ydata",[6e-3 1e-1],...
             "linewidth",1,...
             "linestyle","--",...
             "color","r");
        grid on;
        axis([1e2 1e8 6e-3 1e-1]);
        xlabel("Reynolds number \itRe");
        ylabel("Darcy friction factor \itf");
        set(gca,...
           "fontsize",14,...
           "box","on",...
           "ytick",[6e-3,8e-3,1e-2,2e-2,4e-2,6e-2,8e-2,1e-1],...
           "xtick",[1e2,1e3,1e4,1e5,1e6,1e7,1e8]);
        hold off;
    end
end

