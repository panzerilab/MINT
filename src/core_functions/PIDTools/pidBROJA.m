function [pid_v, table_prob] = pidBROJA(pdf_dirty)
% *function [pid_v, table_prob] = pidBROJA(pdf_dirty)*
%
% ### Description
% Compute the Partial Information Decomposition (PID) according to the BROJA-2 method. The function evaluates the unique, shared,
% and synergistic information contributions between two input variables (X, Y) with respect to a target variable (Z).
%
% ### Inputs:
% - pdf_dirty: A probability distribution function array (nX x nY x nZ), representing the joint probability of the variables X, Y, and Z
%
% ### Outputs:
% - pid_v: A vector containing the information components in the following order: Shared information (SI), unique information of Y (UIY), unique information of Z (UIZ), and complementary information (CI).
% - table_prob: The non-zero probabilities used in the computation
%
% ### Further notes:
% The function uses the BROJA-2 method to optimize the PID measures.

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
    msg = 'Not enough input';
    error('pidbroja:notEnoughInput', msg);
   
end 
if iscell(pdf_dirty)
    pdf_dirty = pdf_dirty{1};
end

if sum(pdf_dirty) == 0
    msg = 'sum of pdf cannot be zero';
    error('pidbroja:InvalidInput', msg);
end

if any(pdf_dirty < 0, "all") 
    msg = 'negative values in pdf are not allowed';
    error('pidbroja:InvalidInput', msg);
end

if any(isnan(pdf_dirty), "all")
    msg = 'pdf contains NaNs. Aborting.';
    error('pidbroja:NaNInput', msg);
end

if size(pdf_dirty, 1)==1   
    warning("Entropy of target is zero.")
    pid_v = [0, 0, 0, 0];
    table_prob = 0;
    return
end

prob_xyz = pdf_dirty / sum(pdf_dirty, 'all');
prob_xyz = prob_xyz .* (prob_xyz > 1e-300);

marg_xy = sum(prob_xyz, 3);
marg_xz = sum(prob_xyz, 2);
marg_xy = reshape(marg_xy, [size(marg_xy, 1) size(marg_xy, 2)]);
marg_xz = reshape(marg_xz, [size(marg_xz, 1) size(marg_xz, 3)]);

conic_solver = Solve_w_ECOSarrays(marg_xy, marg_xz);
conic_solver.create_model();
conic_solver.solve();
conic_solver.addqdistr();
conic_solver.create_qxyz();
entropy_X     = conic_solver.entropy_X2();
condent       = conic_solver.condentropy2();
condent__orig = conic_solver.condentropy_orig2(prob_xyz);
condYmutinf   = conic_solver.condYmutinf2();
condZmutinf   = conic_solver.condZmutinf2();
bits = 1/log(2);

si  = bits * (entropy_X  - condent - condZmutinf - condYmutinf) ;
uiy = bits * condZmutinf;
uiz = bits * condYmutinf;
ci  = bits * (condent - condent__orig);

%joint   = mutualInformationXYZ(prob_xyz);
%qjoint  = mutualInformationXYZ(conic_solver.triplet_nonzeroqpos);
%single1 = mutualInformation(marg_xy);
%single2 = mutualInformation(marg_xz);
%si  = single1 + single2 - qjoint;%bits * (entropy_X  - condent - condZmutinf - condYmutinf) ;
%uiy = single1 - si;%bits * condZmutinf;
%uiz = single2 - si;%bits * condYmutinf;
%ci  = joint - qjoint;%bits * (condent - condent__orig);

pid_v = [si uiy uiz ci];
table_prob = conic_solver.triplet_nonzeroqpos;
end

function MI = mutualInformation(pxy)
    % Calculate marginal probabilities
    px = sum(pxy, 2); % Marginal probability density of X
    py = sum(pxy, 1); % Marginal probability density of Y

    % Initialize mutual information
    MI = 0;

    % Loop through each value of x and y
    for i = 1:size(pxy, 1)
        for j = 1:size(pxy, 2)
            if pxy(i,j) ~= 0 && px(i) ~= 0 && py(j) ~= 0
                MI = MI + pxy(i,j) * log2(pxy(i,j) / (px(i) * py(j)));
            end
        end
    end
end

function MIX = mutualInformationXYZ(pxyz)
    % Calculate marginal probabilities
    px  = sum(sum(pxyz, 2), 3); % Marginal probability density of X
    pyz = squeeze(sum(pxyz, 1)); % Joint probability density of Y and Z

    % Initialize mutual information
    MIX = 0;

    % Loop through each value of x, y, and z
    for i = 1:size(pxyz, 1)
        for j = 1:size(pxyz, 2)
            for k = 1:size(pxyz, 3)
                if pxyz(i,j,k) ~= 0 && px(i) ~= 0 && pyz(j,k) ~= 0
                    MIX = MIX + pxyz(i,j,k) * log2(pxyz(i,j,k) / (px(i) * pyz(j,k)));
                end
            end
        end
    end
end
