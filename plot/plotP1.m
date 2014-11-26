function plotP1(c4n, n4e, x, OPTtitle)
%% Draw a P1-function.
%   plotP1(c4n, n4e, x, OPTtitle) draws the P1-function defined by the grid
%                                 (c4n, n4e) and the basis coefficients (x). 
%                                 The input argument OPTtitle is
%                                 optional, it sets the title of the figure. 
%                                 The default value is empty.

    if(size(n4e,1)>2000)
        trisurf(n4e,c4n(:,1),c4n(:,2),x,'EdgeColor','none');
    else
        trisurf(n4e,c4n(:,1),c4n(:,2),x);
    end

    if nargin == 4
            title(OPTtitle);
        else
            title('');
    end
    drawnow;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2009 
% Numerical Analysis Group
% Prof. Dr. Carsten Carstensen
% Humboldt-University  
% Departement of Mathematics
% 10099 Berlin
% Germany
%
% This file is part of AFEM.
%
% AFEM is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% AFEM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
