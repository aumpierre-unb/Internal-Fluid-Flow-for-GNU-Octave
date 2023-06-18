# Copyright (C) 2022 2023 Alexandre Umpierre
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

function f=Re2f(Re,eps=0,fig=false)
    # Syntax:
    # -- f=Re2f(Re[,eps][,fig])
    #
    # Re2f computes
    #  the Darcy friction f factor given
    #  the Reynolds number Re and
    #  the pipe's relative roughness eps.
    # By default, pipe is assumed to be smooth.
    #  Relative roughness is reset to eps = 0.05, if eps > 0.05.
    # If fig = true is given, a schematic Moody diagram
    #  is plotted as a graphical representation
    #  of the solution.
    # Re2f is a main function of
    #  the internal-fluid-flow toolbox for GNU Octave.
    #
    # Examples:
    # # Compute the Darcy friction factor given
    # # the Reynolds number is 120,000 and
    # # the relative roughness is 0.001:
    # f=Re2f(Re=120e3,eps=1e-3)
    #
    # # Compute the Darcy friction factor given
    # # the Reynolds number is 120,000
    # # for a smooth pipe and
    # # displays a schematic Moody diagram:
    # f=Re2f(Re=120e3,:,true)
    #
    # See also: f2Re, h2fRe.
    if eps>5e-2
        eps=5e-2;
    end
    if Re<2.3e3
        f=64/Re;
    else
        foo=@(f) (1/f^.5+2*log10(eps/3.7+2.51/Re/f^.5));
        f=newtonraphson(foo,1e-2,1e-4);
    end
    if fig
        figure;
        hold on;
        x=[5e-2 2.5e-2 1e-2 3e-3 1e-3 3e-4 1e-4];
        for i=1:length(x)
            if abs(x(i)-eps)>eps/10
                turbulent(x(i),'k',1);
                feps=(-2*log10(x(i)/3.7))^-2;
                text(8e7,feps*1.07,num2str(x(i),4),...
                     'color','k',...
                     'fontsize',11,...
                     'horizontalalignment','right');
            end
        end
        rough('.-b',1.5);
        if eps~=0
            smooth('.-b',1.5);
            text(7e6,8e-3,'Smooth pipe',...
                 'color','b',...
                 'fontsize',11,...
                 'horizontalalignment','right');
            text(4e4,7.6e-2,'Fully rough flow',...
                 'color','b',...
                 'fontsize',11);
        else
            text(7e6,8e-3,'Smooth pipe',...
                 'color','r',...
                 'fontsize',11,...
                 'horizontalalignment','right');
            text(4e4,7.5e-2,'Fully rough flow',...
                 'color','b',...
                 'fontsize',11);
        end
        if Re<2.3e3
            laminar('r',2);
            turbulent(eps,'k',1);
            feps=(-2*log10(eps/3.7))^-2;
            text(9e6,feps*1.07,num2str(eps,4),...
                 'color','k',...
                 'fontsize',11,...
                 'horizontalalignment','right');
        else
            laminar('k',1);
            turbulent(eps,'r',2);
            feps=(-2*log10(eps/3.7))^-2;
            text(9e6,feps*1.07,num2str(eps,4),...
                 'color','r',...
                 'fontsize',11,...
                 'horizontalalignment','right');
        end
        loglog(Re,f,'or',...
               'markersize',8,...
               'markerfacecolor','r');
        line('xdata',[Re Re],...
             'ydata',[6e-3 1e-1],...
             'linewidth',1,...
             'linestyle','--',...
             'color','r');
        grid on;
        axis([1e2 1e8 6e-3 1e-1]);
        xlabel('Reynolds Number \itRe');
        ylabel('Darcy friction factor \itf');
        set(gca,...
           'fontsize',14,...
           'box','on',...
           'ytick',[6e-3,8e-3,1e-2,2e-2,4e-2,6e-2,8e-2,1e-1],...
           'xtick',[1e2,1e3,1e4,1e5,1e6,1e7,1e8]);
        hold off;
    end
end

