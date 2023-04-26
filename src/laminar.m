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

function laminar(t,w)
    # Syntax:
    # -- laminar(t,w)
    #
    # laminar produces a line
    #  of color t and width w
    #  that relates the Reynolds number and
    #  the Darcy friction factor according to
    #  the Poiseuille condition for laminar flow
    #  on a schematic Moody diagram.
    # laminar is an internal function of
    #  the internal-fluid-flow toolbox for GNU Octave.
    line('xdata',[5e2 4e3],...
         'ydata',[64/5e2 64/4e3],...
         'linewidth',w,...
         'color',t);
end
