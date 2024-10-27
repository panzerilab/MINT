function pid_v = pidimin(pdf_dirty)
% function pid_v = pidimin(pdf_dirty)
%
% ### Description
% Compute partial information decomposition (PID) values for a given probability distribution.
%
% ### Inputs:
% - pdf_dirty: n X m X k array containing the joint probability distribution of three variables. The values must sum to one and cannot contain negative or NaN values.
%
% ### Outputs:
% - pid_v: 1x4 vector containing the PID values for different informational components:
%     - Synergy: shared information that only the combination of all variables provides.
%     - Unique1: information uniquely contributed by the first variable.
%     - Unique2: information uniquely contributed by the second variable.
%     - Redundancy: information that is redundantly contributed by any of the variables.
%
% #### Dependency
% This function depends on pid_lattice which calculates the PID components based on the provided probability distribution and desired information structure (like synergy, redundancy, etc.).

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

if nargin < 1
    msg = 'not enough input arguments.';
    error('pidimin:notEnoughInput', msg);
end

if sum(pdf_dirty) == 0
    msg = 'sum of pdf cannot be zero';
    error('pidimin:InvalidInput', msg);
end

if any(pdf_dirty < 0, "all")
    msg = 'negative values in pdf are not allowed';
    error('pidimin:InvalidInput', msg);
end

if any(isnan(pdf_dirty), "all")
    msg = 'pdf contains NaNs. Aborting.';
    error('pidimin:NaNInput', msg);
end

if size(pdf_dirty, 1)==1
    warning("Entropy of target is zero.")
    pid_v = [0, 0, 0, 0];
    return
end

prob_xyz = pdf_dirty / sum(pdf_dirty, 'all');

p = prob_xyz .* (prob_xyz > 1e-300);

new_order = [2:ndims(p) 1];
p = permute(p, new_order);

si = {1,2};
u1 = {1};
u2 = {2};
ci = {[1 2]};

pidLattice = pid_lattice(ndims(p)-1);
pid_v = [];
pid_v = [pid_v pidLattice.calculate_atom(p, 3, si)];
pid_v = [pid_v pidLattice.calculate_atom(p, 3, u1)];
pid_v = [pid_v pidLattice.calculate_atom(p, 3, u2)];
pid_v = [pid_v pidLattice.calculate_atom(p, 3, ci)];

% for node=1:size(pidLattice.lat,2)
%     disp(pidLattice.lat{1,node})
% end

end