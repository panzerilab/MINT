function atoms = cfit_core(p_S4, opts)
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

pidLattice3 = pid_lattice(3);
pidLattice4 = pid_lattice(4);

node_location3 = {1, 2};
node_location4 = {1, 2, 3};
p_Yt4 = permute(p_S4, [5 2 3 4 1]);
p_S3  = squeeze(sum(p_S4,4));
p_Yt3 = permute(p_S3, [4 2 3 1]);
atoms = zeros(2,2);
atoms(1,1) =  pidLattice3.calculate_atom(p_S3,  4, node_location3);
atoms(1,2) =  pidLattice3.calculate_atom(p_Yt3, 4, node_location3);
           
atoms(2,1) =  pidLattice4.calculate_atom(p_S4,  5, node_location4);
atoms(2,2) =  pidLattice4.calculate_atom(p_Yt4, 5, node_location4);
end