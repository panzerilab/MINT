function atoms = fit_core(p_S)
% Copyright (C) 2024 Gabriel Matias Lorenz, Nicola Marie Engel
% This file is part of MINT.
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

node_location = {1,2};
pidLattice = pid_lattice(3);
p_Yt = permute(p_S, [4 2 3 1]);
atoms = zeros(2,1);
atoms(1) = pidLattice.calculate_atom(p_S,  4, node_location);
atoms(2) = pidLattice.calculate_atom(p_Yt, 4, node_location);
end