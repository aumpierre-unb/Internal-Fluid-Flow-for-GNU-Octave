#{
Copyright @copyright{} 2022 Alexandre Umpierre

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

function [Re,f]=hQeps2fRe(h,Q,L,rr,g,mu,rho,s=0)
    # [Re,f]=hQeps2fRe(h,Q,L,rr,g,mu,rho[,s]) computes
    # the Reynolds number Re and
    # the Darcy friction factor f, given
    # the head loss h,
    # the volumetric flow Q,
    # the tube lenght L,
    # the relative roughness rr,
    # the gravitational accelaration g,
    # the fluid's dynamic viscosity mu and
    # the fluid's density rho.
    # If s=1 is given,a schematic Moody diagram
    #  is plotted as a graphical representation
    #  of the computation.
    #
    # # e.g. (computes Re and f, pops no plot)
    # h=40;Q=8500;L=2500;rr=0.0025;g=981;mu=.0089;rho=.989;
    # [Re,f]=hQeps2fRe(h,Q,L,rr,g,mu,rho)
    # D=Q*rho/(pi/4)/Re/mu
    # thk=rr*D
    # v=Q/(pi/4*D^2)
    #
    # # e.g. (computes Re and f and displays plot)
    # [Re,f]=hQeps2fRe(40,8600,2500,0.0025,981,0.0089,0.989,1)
    #
    # See also: Re2f, f2Re, hDeps2fRe, hveps2fRe, hvthk2fRe, hQthk2fRe
    P=2*g*h*Q^3/(pi/4)^3/(mu/rho)^5/L;
    Re=(P/64)^(1/4);
    f=64/Re;
    if Re>25e2
        Re=1e4;
        f=Re2f(Re,rr);
        while abs(f-P/Re^5)/f>5e-3
            if f-P/Re^5<0 Re=Re*1.02;
            else Re=Re*0.98;end
            f=Re2f(Re,rr);
        end
    end
    if s==1
        figure
        laminar('k')
        hold on,turb(rr,'k')
        hold on,turb(rr*3,'k')
        hold on,turb(rr*10,'k')
        hold on,turb(rr/3,'k')
        hold on,turb(rr/10,'k')
        hold on,rough('b')
        hold on,loglog(Re,f,'dr')
        hold on,loglog([Re/10 Re*10],[P/(Re/10)^5 P/(Re*10)^5],'--r')
        grid on
        axis([1e2 1e7 1e-2 1e-1])
        xlabel('{\itRe} = {\it\rho}{\ituD}/{\it\mu}')
        ylabel('{\itf} = {\ith} / ({\itv}^2/{\itg} {\itL}/{\itD})')
        set(gca,'fontsize',14)
    end
end

function laminar(t)
    Re=[5e-2 4e3];
    f=64 ./ Re;
    loglog(Re,f,t);
end

function turb(rr,t)
    Re=[];
    f=[];
    N=50;
    for i=1:N
        w=log10(2e3)+i*(log10(1e8)-log10(2e3))/N;
        Re=[Re;10^w];
        foo=@(f) 1/sqrt(f)+2*log10(rr/3.7+2.51/Re(end)/sqrt(f));
        f=[f;bissecao(foo,1e-2,1e-1,1e-4)];
    end
    loglog(Re,f,t);
end

function rough(t)
    rr=[];
    f=[];
    Re=[];
    N=30;
    for i=1:N
        w=log10(4e-5)+i*(log10(5e-2)-log10(4e-5))/N;
        rr=[rr;10^w];
        f=[f;1.02*(2*log10(3.7/rr(end)))^-2];
        z=f2Re(f(end),rr(end));
        Re=[Re;z(2)];
    end
    loglog(Re,f,t);
end

function x2=bissecao(f,x1,x2,tol)
  while abs(f(x2))>tol
    x=(x1+x2)/2;
    if f(x)*f(x1)>0
      x1=x;
    else
      x2=x;
    end
  end
end

