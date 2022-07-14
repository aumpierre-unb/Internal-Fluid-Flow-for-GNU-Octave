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

function Re=f2Re(f,eps,s=0)
    # Re=f2Re(f,eps[,s]) computes
    # the Reynolds number Re, given
    # the Darcy friction factor f and
    # the relative roughness eps.
    # If s=1 is given,a schematic Moody diagram
    #  is plotted as a graphical representation
    #  of the computation.
    #
    # e.g.
    # f=0.025;eps=0.002;
    # Re=f2Re(f,eps)
    #
    # See also: Re2f, hDeps2fRe, hveps2fRe, hvthk2fRe, hQeps2fRe, hQthk2fRe
    foo=@(Re) 1/sqrt(f)+...
              2*log10(eps/3.7+2.51/Re/sqrt(f));
    Re=bissecao(foo,1e3,1e8,1e-3);
    if s==1
      figure#clf
      laminar()
      hold on,turb(eps,'k')
      hold on,turb(eps*3,'k')
      hold on,turb(eps*10,'k')
      hold on,turb(eps/3,'k')
      hold on,turb(eps/10,'k')
      hold on,loglog(Re,f,'dk')
      grid on
      axis([1e2 1e7 1e-2 1e-1])
      xlabel('{\itRe} = {\it\rho}{\ituD}/{\it\mu}')
      ylabel('{\itf} = {\ith} / ({\itv}^2/{\itg} {\itL}/{\itD})')
      set(gca,'fontsize',14)
    end
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

