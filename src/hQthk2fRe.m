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
    #  the pipe's length L,
    #  the pipe's roughness thk,
    #  the fluid's density rho,
    #  the fluid's dynamic viscosity mu, and
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
    # [Re,f]=hQthk2fRe(h,Q,L,:,rho,mu)
    #
    # # Compute the Reynolds number Re and
    # # the Darcy friction factor f given
    # # the head loss h = 0.40 m,
    # # the volumetric flow rate Q = 8.6 L/s,
    # # length L = 25 m and
    # # roughness thk = 0.27 mm,
    # # the fluid's dynamic viscosity mu = 0.89 cP and
    # # density rho = 0.989 kg/L, and
    # # display a schematic Moody Diagram:
    # [Re,f]=hQthk2fRe(0.40,8.6e-3,25,2.7e-4,989,8.9e-4,9.81,true)
    #
    # See also: Re2f, f2Re, hDeps2fRe, hveps2fRe, hvthk2fRe, hQeps2fRe.
    P=2*g*h*Q^3/(pi/4)^3/(mu/rho)^5/L;
    islam=true;
    Re=(P/64)^(1/4);
    f=64/Re;
    if Re>2.3e3
        islam=false;
        Re=1e4;
        f=P/Re^5;
        D=rho*Q/Re/mu/(pi/4);
        eps=thk/D;
        f=Re2f(Re,eps);
        while abs(f-P/Re^5)/f>5e-3
            if f-P/Re^5<0 Re=Re*1.02;
            else
                Re=Re*0.98;
                if Re<2.3e3
                    islam=true;
                    Re=(P/64)^(1/4);
                    f=64/Re;
                    D=rho*Q/Re/mu/(pi/4);
                    eps=thk/D;
                    printf("Solution found in extended laminar range.\n");
                    break
                end
            end
            f=P/Re^5;
            D=rho*Q/Re/mu/(pi/4);
            eps=thk/D;
            f=Re2f(Re,eps);
        end
    end
    if fig
        figure;
        if islam
          laminar('r',2);
          hold on;turb(eps,'k',1);
        else
          laminar('k',1);
          hold on;turb(eps,'r',2);
        end
        if eps<1e-4, hold on;turb(1e-5,'k',1);
        else hold on;turb(eps/3,'k',1);
        end
        if eps<1e-4, hold on;turb(1e-4,'k',1);
        else hold on;turb(eps/10,'k',1);
        end
        if eps<1e-4, hold on;turb(1e-3,'k',1);
        elseif eps*3>5e-2, hold on;turb(5e-2,'k',1);
        else hold on;turb(eps*3,'k',1);
        end
        if eps<1e-4, hold on;turb(5e-3,'k',1);
        elseif eps*10>5e-2, hold on;turb(eps/6,'k',1);
        else hold on;turb(eps*10,'k',1);
        end
        hold on;rough('b',1.5);
        if ~eps==0, hold on;smooth('b',1.5); end
        hold on;loglog(Re,f,'or','markersize',10","markerfacecolor","r");
        line('xdata',[(P/6e-3)^(1/5) (P/1e-1)^(1/5)],...
             'ydata',[6e-3 1e-1],...
             'linewidth',1,...
             'linestyle','--',...
             'color','r');
        grid on;
        axis([1e2 1e8 6e-3 1e-1]);
        xlabel('{\itRe} = {(4/\it\pi)} \times {\it\rhoQ} / {\it\muD}');
        ylabel('{\itf} = ({4}/{\it\pi})^3 \times 2{\itghQ}^3{\it\rho}^5 / {\it\mu}^5{\itL} \times {\itRe}^{-5}');
        set(gca,'fontsize',14,'box','on');
    end
end

function laminar(t,w)
    line('xdata',[5e2 4e3],...
         'ydata',[64/5e2 64/4e3],...
         'linewidth',w,...
         'color',t);
end

function turb(eps,t,w)
    Re=[];
    f=[];
    N=51;
    for i=1:N
        u=log10(2e3)+(i-1)*(log10(1e8)-log10(2e3))/(N-1);
        Re=[Re;10^u];
        foo=@(f) 1/sqrt(f)+2*log10(eps/3.7+2.51/Re(end)/sqrt(f));
        f=[f;bissection(foo,6e-4,1e-1,1e-4)];
    end
    loglog(Re,f,t,'linewidth',w);
end

function smooth(t,w)
    Re=[];
    f=[];
    N=31;
    for i=1:N
        u=log10(2e3)+(i-1)*(log10(1e7)-log10(2e3))/(N-1);
        Re=[Re;10^u];
        foo=@(f) 1/sqrt(f)+2*log10(2.51/Re(end)/sqrt(f));
        f=[f;bissection(foo,6e-3,1e-1,1e-4)];
    end
    loglog(Re,f,t,'linewidth',w);
end

function rough(t,w)
    eps=[];
    f=[];
    Re=[];
    N=31;
    for i=1:N
        u=log10(4e-5)+(i-1)*(log10(5e-2)-log10(4e-5))/(N-1);
        eps=[eps;10^u];
        f=[f;1.01*(2*log10(3.7/eps(end)))^-2];
        z=f2Re(f(end),eps(end));
        Re=[Re;z(end)];
    end
    loglog(Re,f,t,'linewidth',w);
end

function x2=bissection(f,x1,x2,tol)
  while abs(f(x2))>tol
    x=(x1+x2)/2;
    if f(x)*f(x1)>0 x1=x;
    else x2=x;
    end
  end
end

