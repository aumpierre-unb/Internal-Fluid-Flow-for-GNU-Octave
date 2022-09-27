# Copyright (C) 2022 Alexandre Umpierre
#
# This file is part of Internal Fluid Flow Toolbox.
# Internal Fluid Flow Toolbox is free software:
# you can redistribute it and/or modify it under the terms
# of the GNU General Public License (GPL) version 3
# as published by the Free Software Foundation.
#
# Internal Fluid Flow Toolbox is distributed in the hope
# that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the
# GNU General Public License along with this program
# (license GNU GPLv3.txt).
# It is also available at https://www.gnu.org/licenses/.

function [Re]=f2Re(f,eps=0,fig=false)
    # Syntax:
    # [Re]=f2Re(f,[eps[,fig]])
    #
    # f2Re computes
    #  the Reynolds number Re, given
    #  the Darcy friction factor f and
    #  the pipe's relative roughness eps for
    #  for laminar regime and,
    #  when possible, also
    #  for turbulent regime.
    # By default, pipe is assumed to be smooth, eps=0.
    # If eps>5e-2, execution is aborted.
    # If fig=true is given,a schematic Moody diagram
    #  is plotted as a graphical representation
    #  of the solution.
    #
    # Examples:
    # # Compute Reynolds number Re for
    # # Darcy friction factor f = 0.028 and
    # # relative roughness eps = 0.001.
    # # In this case, both laminar and turbulent
    # # solutions are possible:
    # f=2.8e-2;eps=1e-3;
    # Re=f2Re(f,eps)
    #
    # # Compute the Reynolds number Re given
    # # the Darcy friction factor f = 0.028
    # # in a smooth pipe and
    # # displays a schematic Moody Diagram:
    # # In this case, both turbulent and laminar
    # # solutions are possible:
    # Re=f2Re(2.8e-2,:,true)
    #
    # See also: Re2f, hDeps2fRe, hveps2fRe, hvthk2fRe, hQeps2fRe, hQthk2fRe.
    Re=[];
    fD=[];
    if 64/f<2.3e3
        Re=[Re;64/f];
        fD=[fD;f];
    end
    if f>(2*log10(3.7/eps))^-2
        foo=@(Re) 1/sqrt(f)+...
                  2*log10(eps/3.7+2.51/Re/sqrt(f));
        r=bissecao(foo,1e3,1e8,1e-4);
        if r>2.3e3
            Re=[Re;r];
            fD=[fD;f];
        end
    end
    if ~isempty(fD) && fig
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
        hold on;loglog(Re,fD,'dr');
        line('xdata',[1e2 1e8],...
             'ydata',[f f],...
             'linewidth',1,...
             'linestyle','--',...
             'color','r');
        grid on;
        axis([1e2 1e8 6e-3 1e-1]);
        xlabel('{\itRe}');
        ylabel('{\itf}');
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

