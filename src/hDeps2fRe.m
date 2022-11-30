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

function [Re,f]=hDeps2fRe(h,D,L,eps=0,rho=0.997,mu=9.1e-3,g=981,fig=false)
    # Syntax:
    # [Re,f]=hDeps2fRe(h,D,L,eps,g,mu,rho[,fig=false])
    #
    # hDeps2fRe computes
    #  the Reynolds number Re and
    #  the Darcy friction factor f given
    #  the head loss h,
    #  the pipe's hydraulic diameter D,
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
    # hDeps2fRe is a main function of
    #  the internal-fluid-flow toolbox for GNU Octave.
    #
    # Examples:
    # # Compute the Reynolds number Re and
    # # the Darcy friction factor f given
    # # the head loss h = 0.40 m,
    # # the pipe's hydraulic diameter D = 10 cm,
    # # length L = 25 m and
    # # relative roughness eps = 0.0027,
    # # for water:
    # h=40;D=10;L=2.5e3;eps=2.7e-3; # inputs in cgs units
    # [Re,f]=hDeps2fRe(h,D,L,eps)
    # thk=eps*D # pipe's roughness in cm
    # v=Re*mu/rho/D # flow speed in cm/s
    # Q=v*(pi/4*D^2) # volumetric flow rate in cc/s
    #
    # # Compute the Reynolds number Re and
    # # the Darcy friction factor f given
    # # the head loss h = 0.40 m,
    # # the pipe's hydraulic diameter D = 10 cm,
    # # length L = 25 m and
    # # the fluid's density rho = 0.989 kg/L and
    # # dynamic viscosity mu = 0.89 cP,
    # # in a smooth pipe:
    # h=40;D=10;L=2.5e3;rho=0.989;mu=8.9e-3; # inputs in cgs units
    # [Re,f]=hDeps2fRe(h,D,L,:,rho,mu)
    #
    # # Compute the Reynolds number Re and
    # # the Darcy friction factor f given
    # # the head loss h = 0.40 m,
    # # the pipe's hydraulic diameter D = 10 cm,
    # # length L = 25 m and
    # # relative roughness eps = 0.0027,
    # # the fluid's density rho = 0.989 kg/L and
    # # dynamic viscosity mu = 0.89 cP, and
    # # display a schematic Moody Diagram:
    # [Re,f]=hDeps2fRe(0.40,0.10,25,2.7e-3,989,8.9e-4,9.81,true)
    #
    # See also: Re2f, f2Re, hveps2fRe, hvthk2fRe, hQeps2fRe, hQthk2fRe.
    K=2*g*h*rho^2*D^3/mu^2/L;
    islam=true;
    Re=K/64;
    f=64/Re;
    if Re>2.3e3
        islam=false;
        Re=1e4;
        f=Re2f(Re,eps);
        while abs(f-K/Re^2)/f>5e-3
            if f-K/Re^2<0 Re=Re*1.02;
            else
                Re=Re*0.98;
                if Re<2.3e3
                    islam=true;
                    Re=K/64;
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
        if eps<1e-4
            hold on;turb(1e-5,'k',1);
            feps=(-2*log10(1e-5/3.7))^-2;
            text(2e7,feps*1.07,num2str(1e-5,4),'color','k','fontsize',11);
        else
            hold on;turb(eps/3,'k',1);
            feps=(-2*log10(eps/3/3.7))^-2;
            text(2e7,feps*1.07,num2str(eps/3,4),'color','k','fontsize',11);
        end
        if eps<1e-4
            hold on;turb(1e-4,'k',1);
            feps=(-2*log10(1e-4/3.7))^-2;
            text(2e7,feps*1.07,num2str(1e-4,4),'color','k','fontsize',11);
        else
            hold on;turb(eps/10,'k',1);
            feps=(-2*log10(eps/10/3.7))^-2;
            text(2e7,feps*1.07,num2str(eps/10,4),'color','k','fontsize',11);
        end
        if eps<1e-4
            hold on;turb(1e-3,'k',1);
            feps=(-2*log10(1e-3/3.7))^-2;
            text(2e7,feps*1.07,num2str(1e-3,4),'color','k','fontsize',11);
        elseif eps*3>5e-2
            hold on;turb(5e-2,'k',1);
            feps=(-2*log10(5e-2/3.7))^-2;
            text(2e7,feps*1.07,num2str(5e-2,4),'color','k','fontsize',11);
        else
            hold on;turb(eps*3,'k',1);
            feps=(-2*log10(eps*3/3.7))^-2;
            text(2e7,feps*1.07,num2str(eps*3,4),'color','k','fontsize',11);
        end
        if eps<1e-4
            hold on;turb(5e-3,'k',1);
            feps=(-2*log10(5e-3/3.7))^-2;
            text(2e7,feps*1.07,num2str(5e-3,4),'color','k','fontsize',11);
        elseif eps*10>5e-2
            hold on;turb(eps/1.5,'k',1);
            feps=(-2*log10(eps/1.5/3.7))^-2;
            text(2e7,feps*1.07,num2str(eps/1.5,4),'color','k','fontsize',11);
        else
            hold on;turb(eps*10,'k',1);
            feps=(-2*log10(eps*10/3.7))^-2;
            text(2e7,feps*1.07,num2str(eps*10,4),'color','k','fontsize',11);
        end
        hold on;rough('-.b',1.5);
        if ~eps==0
            hold on;smooth('-.b',1.5);
        end
        if islam
            laminar('r',2);
            hold on;turb(eps,'k',1);
            feps=(-2*log10(eps/3.7))^-2;
            text(2e7,feps*1.07,num2str(eps,4),'color','k','fontsize',11);
        else
            laminar('k',1);
            hold on;turb(eps,'r',2);
            feps=(-2*log10(eps/3.7))^-2;
            text(2e7,feps*1.07,num2str(eps,4),'color','r','fontsize',11);
        end
        hold on;loglog(Re,f,'or','markersize',8,'markerfacecolor','r');
        line('xdata',[(K/6e-3)^(1/2) (K/1e-1)^(1/2)],...
        'ydata',[6e-3 1e-1],...
        'linewidth',1.5,...
        'linestyle','--',...
        'color','r');
        grid on;
        axis([1e2 1e8 6e-3 1e-1]);
        xlabel('{\itRe} = {\it\rhouD} / {\it\mu}');
        ylabel('{\itf} = 2{\itghD}^3{\it\rho}^2 / {\it\mu}^2{\itL} \times {\itRe}^{-2}');
        set(gca,...
        'fontsize',14,...
        'box','on',...
        'ytick',[6e-3,8e-3,1e-2,2e-2,4e-2,6e-2,8e-2,1e-1],...
        'xtick',[1e2,1e3,1e4,1e5,1e6,1e7,1e8]);
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

