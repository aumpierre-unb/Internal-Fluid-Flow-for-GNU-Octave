#  Copyright (C) 2022 2023 Alexandre Umpierre
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

function smooth(t,w)
    # Syntax:
    # smooth(t,w)
    #
    # smooth produces a line
    #  of color t and width w
    #  that represents the limit for
    #  the smooth pipes
    #  on a schematic Moody diagram.
    # smooth is an internal function of
    #  the internal-fluid-flow toolbox for GNU Octave.
    Re=[];
    f=[];
    N=30;
    for n=1:N
        u=log10(2.3e3)+(n-1)*(log10(1e7)-log10(2.3e3))/(N-1);
        Re=[Re;10^u];
        f=[f;Re2f(Re(end),eps)];
    end
    loglog(Re,f,t,...
           'linewidth',w);
end
