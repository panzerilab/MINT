function atoms = ii_core(p_src, opts)
% Copyright (C) 2024 Gabriel Matias Lorenz, Nicola Marie Engel
% This file is part of MINT.
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses

p_crs = permute(p_src, [3 2 1]);
atoms = zeros(2,1);
atoms_1 = feval(opts.function, p_src);
atoms_2 = feval(opts.function, p_crs);
atoms(1) = atoms_1(1);
atoms(2) = atoms_2(1);
end
