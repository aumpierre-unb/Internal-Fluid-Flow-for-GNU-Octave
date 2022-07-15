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

function [Re]=f2Re(f,eps=2e-3+1e-10,s=0)
    # [Re]=f2Re(f,[eps[,s]]) computes
    # the Reynolds number Re, given
    # the Darcy friction factor f and
    # the relative roughness eps for
    # for laminar regime and,
    # when possible, also
    # for turbulent regime.
    # By default eps=2e-3.
    # If eps>5e-2, execution is aborted.
    # If s=1 is given,a schematic Moody diagram
    #  is plotted as a graphical representation
    #  of the computation.
    #
    # # e.g. (computes f for both laminar and
    # # turbulent regimes if possible, pops no plot.
    # f=0.025;eps=0.002;
    # Re=f2Re(f,eps)
    #
    # # e.g. (computes f for laminar regime
    # # and displays plot)
    # Re=f2Re(0.025,:,1)
    #
    # # e.g. (computes f for both laminar and
    # # turbulent regimes if possible,
    # # and displays plot)
    # Re=f2Re(0.025,0.002,1)
    #
    # See also: Re2f, hDeps2fRe, hveps2fRe, hvthk2fRe, hQeps2fRe, hQthk2fRe
    if eps>5e-2 abort end
    Re=[];
    fD=[];
    if 64/f<3000
        Re=[Re;64/f];
        fD=[fD;f];
    end
    if f>(2*log10(3.7/eps))^-2 && eps~=2e-3+1e-10
        foo=@(Re) 1/sqrt(f)+...
                  2*log10(eps/3.7+2.51/Re/sqrt(f));
        Re=[Re;bissecao(foo,1e3,1e8,1e-4)];
        fD=[fD;f];
    end
    if s==1
        figure
        laminar('k')
        hold on,turb(eps,'k')
        hold on,turb(eps*3,'k')
        hold on,turb(eps*10,'k')
        hold on,turb(eps/3,'k')
        hold on,turb(eps/10,'k')
        hold on,rough('b')
        hold on,loglog(Re,fD,'dr')
        grid on
        axis([1e2 1e8 6e-3 1e-1])
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

function turb(eps,t)
    Re=[];
    f=[];
    N=50;
    for i=1:N
        w=log10(2e3)+i*(log10(1e8)-log10(2e3))/N;
        Re=[Re;10^w];
        foo=@(f) 1/sqrt(f)+2*log10(eps/3.7+2.51/Re(end)/sqrt(f));
        f=[f;bissecao(foo,6e-3,1e-1,1e-4)];
    end
    loglog(Re,f,t);
end

function rough(t)
    eps=[];
    f=[];
    Re=[];
    N=30;
    for i=1:N
        w=log10(4e-5)+i*(log10(5e-2)-log10(4e-5))/N;
        eps=[eps;10^w];
        f=[f;1.02*(2*log10(3.7/eps(end)))^-2];
        z=f2Re(f(end),eps(end));
        Re=[Re;z(end)];
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

