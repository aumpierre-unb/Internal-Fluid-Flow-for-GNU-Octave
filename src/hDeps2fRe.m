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

function [Re,f]=hDeps2fRe(h,D,L,eps,g,mu,rho,s=0)
    # [Re,f]=hDeps2fRe(h,D,L,eps,g,mu,rho[,s]) computes
    # the Reynolds number Re and
    # the Darcy friction factor f, given
    # thehead loss h,
    # the hydraullic diameter D,
    # the tube lenght L,
    # the relative roughness eps,
    # the gravitational accelaration g,
    # the fluid's dynamic viscosity mu and
    # the fluid's density rho.
    # If s=1 is given,a schematic Moody diagram
    #  is plotted as a graphical representation
    #  of the computation.
    #
    # e.g.
    # h=40;D=10;L=2500;thk=0.025;g=981;mu=.0089;rho=.989;
    # eps=thk/D
    # [Re,f]=hDeps2fRe(h,D,L,eps,g,mu,rho)
    # v=Re*mu/rho/D
    # Q=v*(pi/4*D^2)
    #
    # See also: Re2f, f2Re, hveps2fRe, hvthk2fRe, hQeps2fRe, hQthk2fRe
    K=2*g*h*rho^2*D^3/mu^2/L;
    Re=K/64;
    f=64/Re;
    if Re>25e2
        Re=1e4;
        f=Re2f(Re,eps);
        while abs(f-K/Re^2)/f>5e-3
            if f-K/Re^2<0 Re=Re*1.02;
            else Re=Re*0.98;end
            f=Re2f(Re,eps);
        end
    end
    if s==1
      figure#clf
      laminar()
      hold on,turb(eps,'k')
      hold on,turb(eps*3,'k')
      hold on,turb(eps*10,'k')
      hold on,turb(eps/3,'k')
      hold on,turb(eps/10,'k')
      hold on,loglog(Re,f,'dk')
      hold on,loglog([Re/10 Re*10],[K/(Re/10)^2 K/(Re*10)^2],'--k')
      grid on
      axis([1e2 1e7 1e-2 1e-1])
      xlabel('{\itRe} = {\it\rho}{\ituD}/{\it\mu}')
      ylabel('{\itf} = {\ith} / ({\itv}^2/{\itg} {\itL}/{\itD})')
      set(gca,'fontsize',14)
    end
end

function laminar()
    Re=[5e-2 4e3];
    f=64 ./ Re;
    loglog(Re,f,'k');
end

function turb(eps,t)
    Re=[];
    f=[];
    N=50;
    for i=1:N
        w=log10(2e3)+i*(log10(1e7)-log10(2e3))/N;
        Re=[Re;10^w];
        foo=@(f) 1/sqrt(f)+2*log10(eps/3.7+2.51/Re(end)/sqrt(f));
        f=[f;bissecao(foo,1e-2,1e-1,1e-4)];
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
