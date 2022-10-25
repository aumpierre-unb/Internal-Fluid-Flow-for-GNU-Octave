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

function [f]=Re2f(Re,eps=0,fig=false)
    # Syntax:
    # [f]=Re2f(Re,[eps[,fig]])
    #
    # Re2f computes
    #  the Darcy friction f factor given
    #  the Reynolds number Re and
    #  the pipe's relative roughness eps.
    # By default, pipe is assumed to be smooth, eps = 0.
    # If eps > 5e-2, eps is reset to eps = 5e-2.
    # If fig = true is given, a schematic Moody diagram
    #  is plotted as a graphical representation
    #  of the solution.
    # Re2f is a main function of
    #  the internal-fluid-flow toolbox for GNU Octave.
    #
    # Examples:
    # # Compute the Darcy friction factor f given
    # # the Reynolds number Re = 120,000 and
    # # the relative roughness eps = 0.001:
    # Re=1.2e5;eps=0.001;
    # f=Re2f(Re,eps)
    #
    # # Compute the Darcy friction factor f given
    # # the Reynolds number Re = 120,000
    # # in a smooth pipe and
    # # displays a schematic Moody Diagram:
    # f=Re2f(1.2e5,:,true)
    #
    # See also: f2Re, hDeps2fRe, hveps2fRe, hvthk2fRe, hQeps2fRe, hQthk2fRe.
    if eps>5e-2
        eps=5e-2;
    end
    if Re<2.3e3
        f=64/Re;
    else
        foo=@(f) 1/sqrt(f)+...
                 2*log10(eps/3.7+2.51/Re/sqrt(f));
        f=bissection(foo,1e-2,1e-1,1e-4);
    end
    if fig
        figure;
        if Re<2.3e3
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
        line('xdata',[Re Re],...
             'ydata',[6e-3 1e-1],...
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
        f=[f;bissection(foo,6e-4,1e-1,1e-4)];
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
        f=[f;bissection(foo,6e-3,1e-1,1e-4)];
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

function x2=bissection(f,x1,x2,tol)
  while abs(f(x2))>tol
    x=(x1+x2)/2;
    if f(x)*f(x1)>0 x1=x;
    else x2=x;
    end
  end
end

