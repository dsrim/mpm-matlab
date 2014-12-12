function e4n = computeE4n(n4e)
%% computeE4n - Elements for nodes.
%   computeE4n(n4e) returns a sparse matrix in which the entry (j,k) contains
%               the number of the element that has the nodes j and k in
%               counterclockwise order (or 0 if no such element exists).
%               n4e is as specified in the documentation.
%
%   See also: computeE4s

    if isempty(n4e)
        e4n = [];
        return;
    end

    %% Compute e4n.
    % Create a list of all sides in the decomposition and build a sparse
    % matrix such that each side computes its proper element number.
    allSides = [n4e(:,[1 2]); n4e(:,[2 3]); n4e(:,[3 1])];
    nrElems = size(n4e,1);
    N = max(max(n4e));
    elemNumbers = [1:nrElems 1:nrElems 1:nrElems];
    e4n = sparse(allSides(:,1),allSides(:,2),elemNumbers,N,N);
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
