function pid_v = pidimminsb(inputs)
% *function pid_v = pidimminsb(inputs_1d)*
%
% ### Description
% Compute Partial Information Decomposition (PID) values from a given probability distribution (pdf_dirty).
%
% Partial Information Decomposition separates the information that multiple sources provide about a target into unique, redundant, and synergistic components. 
% This function takes a 3D probability distribution (joint pdf of three variables) and calculates the PID values.
%
% ### Inputs:
% - pdf_dirty: A 3-dimensional array representing the joint probability density function of three variables, typically denoted as X, Y, and Z. 
%
% ### Outputs:
% - pid_v: A vector containing four values representing the components of PID:
% - Redundant information (common information both sources provide about the target)
% - Unique information from the first variable (information only the first source provides about the target)
% - Unique information from the second variable (information only the second source provides about the target)
% - Synergistic information (information that is only available when both sources are considered together)

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

nsb_opts = struct();
nsb_opts.bias = 'nsb';

I1  = cell2mat(MI({inputs{1}, inputs{end}}, {'I(A;B)'}, nsb_opts));
I2  = cell2mat(MI({inputs{2}, inputs{end}}, {'I(A;B)'}, nsb_opts));
I12 = cell2mat(MI({cat(1, inputs{1}, inputs{2}), inputs{end}}, {'I(A;B)'}, nsb_opts));

nTimepoints = size(inputs{1},2);
pid_v = cell(1,nTimepoints);
red = min([I1;I2], [], 1);
u1 = I1 - red;
u2 = I2 - red;
syn = I12 - u1 - u2 - red;
for i=1:nTimepoints
    pid_v{i} = [red(i), u1(i), u2(i), syn(i)];
end

end



