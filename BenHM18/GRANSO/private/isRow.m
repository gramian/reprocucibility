function tf = isRow(r)
%   isRow:
%       Checks whether r is a nonempty row vector, that is, for 
%       [m,n] = size(r), both m = 1 and n > 0 must hold.
%
%
%   For comments/bug reports, please visit the GRANSO GitLab webpage:
%   https://gitlab.com/timmitchell/GRANSO
%
%   isRow.m introduced in GRANSO Version 1.0.
%
% =========================================================================
% |  isRow.m                                                              |
% |  Copyright (C) 2016 Tim Mitchell                                      |
% |                                                                       |
% |  This file is originally from URTM.                                   |
% |                                                                       |
% |  URTM is free software: you can redistribute it and/or modify         |
% |  it under the terms of the GNU Affero General Public License as       |
% |  published by the Free Software Foundation, either version 3 of       |
% |  the License, or (at your option) any later version.                  |
% |                                                                       |
% |  URTM is distributed in the hope that it will be useful,              |
% |  but WITHOUT ANY WARRANTY; without even the implied warranty of       |
% |  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        |
% |  GNU Affero General Public License for more details.                  |
% |                                                                       |
% |  You should have received a copy of the GNU Affero General Public     |
% |  License along with this program.  If not, see                        |
% |  <http://www.gnu.org/licenses/>.                                      |
% =========================================================================
%
% =========================================================================
% |  GRANSO: GRadient-based Algorithm for Non-Smooth Optimization         |
% |  Copyright (C) 2016 Tim Mitchell                                      |
% |                                                                       |
% |  This file is part of GRANSO.                                         |
% |                                                                       |
% |  GRANSO is free software: you can redistribute it and/or modify       |
% |  it under the terms of the GNU Affero General Public License as       |
% |  published by the Free Software Foundation, either version 3 of       |
% |  the License, or (at your option) any later version.                  |
% |                                                                       |
% |  GRANSO is distributed in the hope that it will be useful,            |
% |  but WITHOUT ANY WARRANTY; without even the implied warranty of       |
% |  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        |
% |  GNU Affero General Public License for more details.                  |
% |                                                                       |
% |  You should have received a copy of the GNU Affero General Public     |
% |  License along with this program.  If not, see                        |
% |  <http://www.gnu.org/licenses/>.                                      |
% =========================================================================

    [m,n]   = size(r);
    tf      = m == 1 && n > 0;
end