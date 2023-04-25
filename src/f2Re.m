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

function Re=f2Re(f,eps=0,fig=false)
    # Syntax:
    # Re=f2Re(f,eps=0,fig=false)
    #
    # f2Re computes
    #  the Reynolds number Re given
    #  the Darcy friction factor f and
    #  the pipe's relative roughness eps for
    #  for laminar regime and,
    #  when possible, also
    #  for turbulent regime.
    # By default, pipe is assumed to be smooth.
    #  Relative roughness is reset to eps = 0.05, if eps > 0.05.
    # If fig = true is given, a schematic Moody diagram
    #  is plotted as a graphical representation
    #  of the solution.
    # f2Re is a main function of
    #  the internal-fluid-flow toolbox for GNU Octave.
    #
    # Examples:
    # # Compute Reynolds number Re given
    # # the Darcy friction factor f = 0.028 and
    # # the relative roughness eps = 0.001.
    # # In this case, both laminar and turbulent
    # # solutions are possible:
    # f=2.8e-2;eps=1e-3;
    # Re=f2Re(f,eps)
    #
    # # Compute the Reynolds number Re given
    # # the Darcy friction factor f = 0.028
    # # for a smooth pipe and
    # # plot a schematic Moody diagram
    # # with the solution.
    # # In this case, both turbulent and laminar
    # # solutions are possible:
    # Re=f2Re(2.8e-2,:,true)
    #
    # See also: Re2f, h2fRe.
    if eps>5e-2
        eps=5e-2;
    end
    Re=[];
    fD=[];
    Re_=64/f;
    if Re_<2.3e3
        Re=[Re;Re_];
        fD=[fD;f];
    end
    if f>(2*log10(3.7/eps))^-2
        Re_=2.51/(10^(1/f^.5/-2)-eps/3.7)/f^.5;
        if Re_>2.3e3
            Re=[Re;Re_];
            fD=[fD;f];
        end
    end
    if ~isempty(fD) && fig
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
        rough('-.b',1.5);
        if eps~=0
            smooth('-.b',1.5);
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
        if min(Re)<2.3e3
            laminar('r',2);
        else
            laminar('k',1);
        end
        if max(Re)>2.3e3
            turbulent(eps,'r',2);
            feps=(-2*log10(eps/3.7))^-2;
            text(2e7,feps*1.07,num2str(eps,4),...
                 'color','r',...
                 'fontsize',11,...
                 'horizontalalignment','right');
        else
            turbulent(eps,'k',1);
        end
        loglog(Re,f,'or',...
               'markersize',8,...
               'markerfacecolor','r');
        line('xdata',[1e2 1e8],...
             'ydata',[f f],...
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

