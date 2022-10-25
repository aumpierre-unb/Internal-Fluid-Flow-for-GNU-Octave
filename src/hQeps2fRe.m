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

function [Re,f]=hQeps2fRe(h,Q,L,eps=0,rho=0.997,mu=9.1e-3,g=981,fig=false)
    # Syntax:
    # [Re,f]=hQeps2fRe(h,Q,L,eps,g,mu,rho[,fig])
    #
    # hQeps2fRe computes
    #  the Reynolds number Re and
    #  the Darcy friction factor f given
    #  the head loss h,
    #  the volumetric flow Q,
    #  the pipe's length L,
    #  the pipe's relative roughness eps,
    #  the fluid's density rho,
    #  the fluid's dynamic viscosity mu, and
    #  the gravitational accelaration g.
    # By default, pipe is assumed to be smooth, eps = 0.
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
    # hQeps2fRe is a main function of
    #  the internal-fluid-flow toolbox for GNU Octave.
    #
    # Examples:
    # # Compute the Reynolds number Re and
    # # the Darcy friction factor f given
    # # the head loss h = 0.40 m,
    # # the volumetric flow rate Q = 8.6 L/s,
    # # length L = 25 m and
    # # relative roughness eps = 0.0027,
    # # for water:
    # h=40;Q=8.6e3;L=2.5e3;eps=2.7e-3; # inputs in cgs units
    # [Re,f]=hQeps2fRe(h,Q,L,eps)
    # D=rho/mu*Q/(%pi/4)/Re # pipe's hydraulic diameter in cm
    # thk=eps*D # pipe's roughness in cm
    # v=Re*mu/rho/D # flow speed in cm/s
    #
    # # Compute the Reynolds number Re and
    # # the Darcy friction factor f given
    # # the head loss h = 0.40 m,
    # # the volumetric flow rate Q = 8.6 L/s,
    # # length L = 25 m and
    # # the fluid's density rho = 0.989 kg/L and
    # # dynamic viscosity mu = 0.89 cP, and
    # # in a smooth pipe:
    # h=40;Q=8.6e3;L=2.5e3;rho=0.989;mu=8.9e-3; # inputs in cgs units
    # [Re,f]=hQeps2fRe(h,Q,L,:,rho,mu)
    #
    # # Compute the Reynolds number Re and
    # # the Darcy friction factor f given
    # # the head loss h = 0.40 m,
    # # the volumetric flow rate Q = 8.6 L/s,
    # # length L = 25 m and
    # # relative roughness eps = 0.0027,
    # # the fluid's dynamic viscosity mu = 0.89 cP and
    # # density rho = 0.989 kg/L, and
    # # display a schematic Moody Diagram:
    # [Re,f]=hQeps2fRe(0.40,8.6e-3,25,2.7e-3,989,8.9e-4,9.81,true)
    #
    # See also: Re2f, f2Re, hDeps2fRe, hveps2fRe, hvthk2fRe, hQthk2fRe.
    P=2*g*h*Q^3/(pi/4)^3/(mu/rho)^5/L;
    islam=true;
    Re=(P/64)^(1/4);
    f=64/Re;
    if Re>2.3e3
        islam=false;
        Re=1e4;
        f=Re2f(Re,eps);
        while abs(f-P/Re^5)/f>5e-3
            if f-P/Re^5<0 Re=Re*1.02;
            else
                Re=Re*0.98;
                if Re<2.3e3
                    islam=true;
                    Re=(P/64)^(1/4);
                    f=64/Re;
                    printf("Solution found in extended laminar range.\n");
                    break
                end
            end
            f=Re2f(Re,eps);
        end
    end
    if fig
        figure;
        if islam
          laminar('r');
          hold on;turb(eps,'k');
        else
          laminar('k');
          hold on;turb(eps,'r');
        end
        if eps<1e-4, hold on;turb(1e-5,'k');
        else hold on;turb(eps/3,'k');
        end
        if eps<1e-4, hold on;turb(1e-4,'k');
        else hold on;turb(eps/10,'k');
        end
        if eps<1e-4, hold on;turb(1e-3,'k');
        elseif eps*3>5e-2, hold on;turb(5e-2,'k');
        else hold on;turb(eps*3,'k');
        end
        if eps<1e-4, hold on;turb(5e-3,'k');
        elseif eps*10>5e-2, hold on;turb(eps/6,'k');
        else hold on;turb(eps*10,'k');
        end
        hold on;rough('b');
        if ~eps==0, hold on;smooth('b'); end
        hold on;loglog(Re,f,'dr');
        line('xdata',[(P/6e-3)^(1/5) (P/1e-1)^(1/5)],...
             'ydata',[6e-3 1e-1],...
             'linewidth',1,...
             'linestyle','--',...
             'color','r');
        grid on;
        axis([1e2 1e8 6e-3 1e-1]);
        xlabel('{\itRe} = {(4/\it\pi)} \times {\it\rhoQ} / {\it\muD}');
        ylabel('{\itf} = ({4}/{\it\pi})^3 \times 2{\itghQ}^3{\it\rho}^5 / {\it\mu}^5{\itL} \times {\itRe}^{-5}');
        set(gca,'fontsize',14);
    end
end

function laminar(t)
    line('xdata',[5e2 4e3],...
         'ydata',[64/5e2 64/4e3],...
         'linewidth',1,...
         'color',t);
end

function turb(eps,t)
    Re=[];
    f=[];
    N=51;
    for i=1:N
        w=log10(2e3)+(i-1)*(log10(1e8)-log10(2e3))/(N-1);
        Re=[Re;10^w];
        foo=@(f) 1/sqrt(f)+2*log10(eps/3.7+2.51/Re(end)/sqrt(f));
        f=[f;bissecao(foo,6e-4,1e-1,1e-4)];
    end
    loglog(Re,f,t);
end

function smooth(t)
    Re=[];
    f=[];
    N=31;
    for i=1:N
        w=log10(2e3)+(i-1)*(log10(1e7)-log10(2e3))/(N-1);
        Re=[Re;10^w];
        foo=@(f) 1/sqrt(f)+2*log10(2.51/Re(end)/sqrt(f));
        f=[f;bissecao(foo,6e-3,1e-1,1e-4)];
    end
    loglog(Re,f,t);
end

function rough(t)
    eps=[];
    f=[];
    Re=[];
    N=31;
    for i=1:N
        w=log10(4e-5)+(i-1)*(log10(5e-2)-log10(4e-5))/(N-1);
        eps=[eps;10^w];
        f=[f;1.01*(2*log10(3.7/eps(end)))^-2];
        z=f2Re(f(end),eps(end));
        Re=[Re;z(end)];
    end
    loglog(Re,f,t);
end

function x2=bissecao(f,x1,x2,tol)
  while abs(f(x2))>tol
    x=(x1+x2)/2;
    if f(x)*f(x1)>0 x1=x;
    else x2=x;
    end
  end
end

