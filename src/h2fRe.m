# Copyright (C) 2022 2023 Alexandre Umpierre
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

function [Re,f]=h2fRe(h,D=NaN,v=NaN,Q=NaN,eps=NaN,k=NaN,L=100,rho=0.997,mu=9.1e-3,g=981,fig=false)
    # Syntax:
    # -- [Re,f]=h2fRe(h,D=double,eps=double[,L=double][,rho=double][,mu=double][,g=double][,fig=logical])
    # -- [Re,f]=h2fRe(h,D=double,k=double[,L=double][,rho=double][,mu=double][,g=double][,fig=logical])
    # -- [Re,f]=h2fRe(h,v=double,eps=double[,L=double][,rho=double][,mu=double][,g=double][,fig=logical])
    # -- [Re,f]=h2fRe(h,v=double,k=double[,L=double][,rho=double][,mu=double][,g=double][,fig=logical])
    # -- [Re,f]=h2fRe(h,Q=double,eps=double[,L=double][,rho=double][,mu=double][,g=double][,fig=logical])
    # -- [Re,f]=h2fRe(h,Q=double,k=double[,L=double][,rho=double][,mu=double][,g=double][,fig=logical])
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
    #  the fluid's dynamic viscosity mu (default mu = 0.0091), and
    #  the gravitational accelaration g (default g = 981).
    # By default, pipe is assumed to be smooth.
    #  Relative roughness is reset to eps = 0.05, if eps > 0.05.
    # Notice that default values are given in the cgs unit system and,
    #  if taken, all other parameters must as well
    #  be given in cgs units.
    # If parameter fig = true is given
    #  a schematic Moody diagram
    #  is plotted as a graphical representation
    #  of the solution.
    # h2fRe is a main function of
    #  the internal-fluid-flow toolbox for GNU Octave.
    #
    # Examples:
    # # Compute the Reynolds number Re and
    # # the Darcy friction factor f given
    # # the head loss h = 40 cm,
    # # the pipe's hydraulic diameter D = 10 cm,
    # # length L = 25 m and
    # # relative roughness eps = 0.0027
    # # for water flow:
    # [Re,f]=h2fRe(40,D=10,L=2.5e3,eps=2.7e-3)
    #
    # # Compute the Reynolds number Re and
    # # the Darcy friction factor f given
    # # the head loss per meter h/L = 1.6 cm/m,
    # # the volumetric flow rate Q = 8.6 L/s,
    # # the fluid's density rho = 0.989 g/cc and
    # # dynamic viscosity mu = 0.89 cP
    # # for a smooth pipe and
    # # show results on a schematic Moody diagram:
    # [Re,f]=h2fRe(1.6,Q=8.6e3,eps=0,rho=0.989,mu=8.9e-3,fig=true)
    #
    # # Compute the Reynolds number Re and
    # # the Darcy friction factor f, given
    # # the head loss h = 0.40 m,
    # # the flow speed v = 1.1 m/s,
    # # the pipe's length L = 25 m
    # # for water flow for a smooth pipe:
    # [Re,f]=h2fRe(40,v=1.1e2,L=2.5e3,k=0)
    #
    # See also: Re2f, f2Re.
    a = isnan([D, v, Q]) != 1
    if sum(a) != 1
        error("h2fRe requires that either\nthe hydraulic diameter,\nthe flow speed or\nthe flow rate\nbe given alone.")
    end
    b = isnan([eps, k]) != 1
    if sum(b) != 1
        error("h2fRe requires that either\nthe pipe's roughness or\nthe pipe's relative roughness\nbe given alone.")
    end
    K=2*g*h*rho^2*D^3/mu^2/L;
    foo=@(f) (1/f^(1/2)+2*log10(eps/3.7+2.51/(K/f)^(1/2)/f^(1/2)));
    f=newtonraphson(foo,1e-2,1e-4);
    Re=(K/f)^(1/2);
    if Re>2.3e3
        islam=false;
    else
        Re=K/64;
        f=64/Re;
        islam=true;
    end
    if fig
        figure;
        hold on;
        x=[5e-2 2.5e-2 1e-2 3e-3 1e-3 3e-4 1e-4];
        for i=1:length(x)
            if abs(x(i)-eps)>eps/10
                turbulent(x(i),'k',1);
                feps=(-2*log10(x(i)/3.7))^-2;
                text(8e7,feps*1.07,num2str(x(i),4),...
                     'color','k',...
                     'fontsize',11,...
                     'horizontalalignment','right');
            end
        end
        rough('-.b',1.5);
        if eps~=0
            smooth('-.b',1.5);
            text(7e6,8e-3,'Smooth pipe',...
                 'color','b',...
                 'fontsize',11,...
                 'horizontalalignment','right');
            text(4e4,7.6e-2,'Fully rough flow',...
                 'color','b',...
                 'fontsize',11);
        else
            text(7e6,8e-3,'Smooth pipe',...
                 'color','r',...
                 'fontsize',11,...
                 'horizontalalignment','right');
            text(4e4,7.5e-2,'Fully rough flow',...
                 'color','b',...
                 'fontsize',11);
        end
        if islam
            laminar('r',2);
        else
            laminar('k',1);
            turbulent(eps,'r',2);
            feps=(-2*log10(eps/3.7))^-2;
            text(8e6,feps*1.07,num2str(eps,4),...
                 'color','r',...
                 'fontsize',11,...
                 'horizontalalignment','right');
        end
        loglog(Re,f,'or',...
               'markersize',8,...
               'markerfacecolor','r');
        line('xdata',[(K/6e-3)^(1/2) (K/1e-1)^(1/2)],...
             'ydata',[6e-3 1e-1],...
             'linewidth',1,...
             'linestyle','--',...
             'color','r');
        grid on;
        axis([1e2 1e8 6e-3 1e-1]);
        xlabel('Reynolds Number \itRe');
        ylabel('Darcy friction factor \itf');
        set(gca,...
           'fontsize',14,...
           'box','on',...
           'ytick',[6e-3,8e-3,1e-2,2e-2,4e-2,6e-2,8e-2,1e-1],...
           'xtick',[1e2,1e3,1e4,1e5,1e6,1e7,1e8]);
        hold off;
    end
end

