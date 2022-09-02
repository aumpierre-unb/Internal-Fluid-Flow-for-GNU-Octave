#{
Copyright (C) 2022 Alexandre Umpierre

This file is part of Internal Fluid Flow Toolbox.
Internal Fluid Flow Toolbox is free software:
you can redistribute it and/or modify it under the terms
of the GNU General Public License (GPL) version 3
as published by the Free Software Foundation.

Internal Fluid Flow Toolbox is distributed in the hope
that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the
GNU General Public License along with this program
(license GNU GPLv3.txt).
It is also available at https://www.gnu.org/licenses/.
#}

function [Re,f]=hvthk2fRe(h,v,L,thk,g,mu,rho,fig=false)
    # [Re,f]=hvthk2fRe(h,v,L,thk,g,mu,rho[,fig]) computes
    # the Reynolds number Re and
    # the Darcy friction factor f, given
    # the head loss h,
    # the flow speed v,
    # the pipe's length L,
    # the pipe's roughness thk,
    # the gravitational accelaration g,
    # the fluid's dynamic viscosity mu and
    # the fluid's density rho.
    # If fig=true is given,a schematic Moody diagram
    #  is plotted as a graphical representation
    #  of the solution.
    #
    # # e.g. Compute the Reynolds number Re and
    # # the Darcy friction factor f, given
    # # the head loss h=40 cm,
    # # the flow speed v=100 cm/s,
    # # the pipe's length L=2500 cm and
    # # roughness thk=0.0025 cm,
    # # the gravitational acceleration g=981 cm/s/s, and
    # # the fluid's dynamic viscosity mu=0.0089 g/cm/s and
    # # density rho=0.989 g/cc:
    # h=40;v=110;L=2500;eps=0.0025;g=981;mu=0.0089;rho=0.989;
    # [Re,f]=hvthk2fRe(h,v,L,eps,g,mu,rho)
    # D=Re*mu/rho/v #pipe's hydraulic diameter
    # eps=thk/D #pipe's relative roughness
    # Q=v*(pi/4*D^2) #volumetric flow rate
    # # Alternatively:
    # [Re,f]=hvthk2fRe(40,100,2500,0.005,981,0.0089,0.989)
    # # This call computes Re and f for the given inputs and
    # # displays a schematic Moody Diagram:
    # [Re,f]=hvthk2fRe(40,100,2500,0.0025,981,0.0089,0.989,true)
    #
    # # e.g. Compute the Reynolds number Re and
    # # the Darcy friction factor f, given
    # # the head loss h=40 cm,
    # # the flow speed v=20 cm/s,
    # # the pipe's length L=2500 cm and
    # # roughness thk=0.0025 cm,
    # # the gravitational acceleration g=981 cm/s/s, and
    # # the fluid's dynamic viscosity mu=0.0089 g/cm/s and
    # # density rho=0.989 g/cc, and
    # # display a schematic Moody Diagram:
    # [Re,f]=hvthk2fRe(40,20,2500,0.0025,981,0.0089,0.989,true)
    #
    # # e.g. Compute the Reynolds number Re and
    # # the Darcy friction factor f, given
    # # the head loss h=40 cm,
    # # the flow speed v=27 cm/s,
    # # the pipe's length L=2500 cm and
    # # roughness thk=0.0025 cm,
    # # the gravitational acceleration g=981 cm/s/s, and
    # # the fluid's dynamic viscosity mu=0.0089 g/cm/s and
    # # density rho=0.989 g/cc, and
    # # display a schematic Moody Diagram:
    # [Re,f]=hvthk2fRe(40,27,2500,0.0025,981,0.0089,0.989,true)
    #
    # See also: Re2f, f2Re, hDeps2fRe, hveps2fRe, hQeps2fRe, hQthk2fRe
    M=2*g*mu*h/v^3/rho/L;
    isturb=true;
    Re=1e4;
    f=M*Re;
    D=Re*mu/v/rho;
    eps=thk/D;
    f=Re2f(Re,eps);
    while abs(f-Re*M)/f>5e-3
        if f-Re*M>0 Re=Re*1.02;
        else
            Re=Re*0.98;
            if Re<2.3e3
                isturb=false;
                Re=sqrt(64/M);
                f=64/Re;
                D=Re*mu/v/rho;
                eps=thk/D;
                break
            end
        end
        f=M*Re;
        D=Re*mu/v/rho;
        eps=thk/D;
        f=Re2f(Re,eps);
    end
    if isturb && sqrt(64/M)<2.3e3
        Re=[sqrt(64/M);Re];
        f=[64/sqrt(64/M);f];
    end
    if size(Re,1)==2
        warning("Laminar and turbulent solutions found.");
    end
    if fig
        figure;
        if min(Re)<2.3e3 laminar('r');
        else laminar('k');
        end
        if max(Re)>2.3e3, hold on;turb(eps,'r');
        else hold on;turb(eps,'k');
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
        line('xdata',[6e-3/M 1e-1/M],...
             'ydata',[6e-3 1e-1],...
             'linewidth',1,...
             'linestyle','--',...
             'color','r');
        grid on;
        axis([1e2 1e8 6e-3 1e-1]);
        xlabel('{\itRe} = {\it\rhouD} / {\it\mu}');
        ylabel('{\itf} = 2{\itghD}^3{\it\rho}^2 / {\it\mu}^2{\itL} \times {\itRe}');
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
        f=[f;bissecao(foo,6e-3,1e-1,1e-4)];
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

