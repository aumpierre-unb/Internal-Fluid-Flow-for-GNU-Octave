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
    # [f]=Re2f(Re,[eps[,fig=false]])
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
        hold on;
        if eps<1e-4
            turb(1e-5,'k',1);
            feps=(-2*log10(1e-5/3.7))^-2;
            text(2e7,feps*1.07,num2str(1e-5,4),'color','k','fontsize',11);
        else
            turb(eps/3,'k',1);
            feps=(-2*log10(eps/3/3.7))^-2;
            text(2e7,feps*1.07,num2str(eps/3,4),'color','k','fontsize',11);
        end
        if eps<1e-4
            turb(1e-4,'k',1);
            feps=(-2*log10(1e-4/3.7))^-2;
            text(2e7,feps*1.07,num2str(1e-4,4),'color','k','fontsize',11);
        else
            turb(eps/10,'k',1);
            feps=(-2*log10(eps/10/3.7))^-2;
            text(2e7,feps*1.07,num2str(eps/10,4),'color','k','fontsize',11);
        end
        if eps<1e-4
            turb(1e-3,'k',1);
            feps=(-2*log10(1e-3/3.7))^-2;
            text(2e7,feps*1.07,num2str(1e-3,4),'color','k','fontsize',11);
        elseif eps*3>5e-2
            turb(5e-2,'k',1);
            feps=(-2*log10(5e-2/3.7))^-2;
            text(2e7,feps*1.07,num2str(5e-2,4),'color','k','fontsize',11);
        else
            turb(eps*3,'k',1);
            feps=(-2*log10(eps*3/3.7))^-2;
            text(2e7,feps*1.07,num2str(eps*3,4),'color','k','fontsize',11);
        end
        if eps<1e-4
            turb(5e-3,'k',1);
            feps=(-2*log10(5e-3/3.7))^-2;
            text(2e7,feps*1.07,num2str(5e-3,4),'color','k','fontsize',11);
        elseif eps*10>5e-2
            turb(eps/1.5,'k',1);
            feps=(-2*log10(eps/1.5/3.7))^-2;
            text(2e7,feps*1.07,num2str(eps/1.5,4),'color','k','fontsize',11);
        else
            turb(eps*10,'k',1);
            feps=(-2*log10(eps*10/3.7))^-2;
            text(2e7,feps*1.07,num2str(eps*10,4),'color','k','fontsize',11);
        end
        hold on;rough('-.b',1.5);
        if ~eps==0
            smooth('-.b',1.5);
            text(7e6,8e-3,'Smooth pipe','color','b','fontsize',11,'horizontalalignment','right');
            text(4e4,7.5e-2,'Fully rough flow','color','b','fontsize',11);
        else
            text(7e6,8e-3,'Smooth pipe','color','r','fontsize',11,'horizontalalignment','right');
            text(4e4,7.5e-2,'Fully rough flow','color','b','fontsize',11);
        end
        if Re<2.3e3
            laminar('r',2);
            turb(eps,'k',1);
            feps=(-2*log10(eps/3.7))^-2;
            text(2e7,feps*1.07,num2str(eps,4),'color','k','fontsize',11);
        else
            laminar('k',1);
            turb(eps,'r',2);
            feps=(-2*log10(eps/3.7))^-2;
            text(2e7,feps*1.07,num2str(eps,4),'color','r','fontsize',11);
        end
        loglog(Re,f,'or','markersize',8,'markerfacecolor','r');
        line('xdata',[Re Re],...
             'ydata',[6e-3 1e-1],...
             'linewidth',1,...
             'linestyle','--',...
             'color','r');
        grid on;
        axis([1e2 1e8 6e-3 1e-1]);
        xlabel('{\itRe}');
        ylabel('{\itf}');
        set(gca,...
           'fontsize',14,...
           'box','on',...
           'ytick',[6e-3,8e-3,1e-2,2e-2,4e-2,6e-2,8e-2,1e-1],...
           'xtick',[1e2,1e3,1e4,1e5,1e6,1e7,1e8]);
        hold off;
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

