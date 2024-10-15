 function pid_v = pidimmi(pdf_dirty)
%%% *function pid_v = pidimmi(pdf_dirty)*
%%%
%%% ### Description
%%% Compute Partial Information Decomposition (PID) values from a given probability distribution (pdf_dirty).
%%%
%%% Partial Information Decomposition separates the information that multiple sources provide about a target into unique, redundant, and synergistic components. 
%%% This function takes a 3D probability distribution (joint pdf of three variables) and calculates the PID values.
%%%
%%% ### Inputs:
%%% - pdf_dirty: A 3-dimensional array representing the joint probability density function of three variables, typically denoted as X, Y, and Z. 
%%%
%%% ### Outputs:
%%% - pid_v: A vector containing four values representing the components of PID:
%%% - Redundant information (common information both sources provide about the target)
%%% - Unique information from the first variable (information only the first source provides about the target)
%%% - Unique information from the second variable (information only the second source provides about the target)
%%% - Synergistic information (information that is only available when both sources are considered together)

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

p = permute(p, [2 3 1]);

MIxyz = mutualInformationXYZ(p);
MIxz  = mutualInformation(squeeze(sum(p,2)));
MIyz  = mutualInformation(squeeze(sum(p,1)));

red = min([MIxz, MIyz]);
u1 = MIxz - red;
u2 = MIyz - red;
syn = MIxyz - u1 - u2 - red;

pid_v = [red u1 u2 syn];
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

function MIZ = mutualInformationXYZ(pxyz)
    % Calculate marginal probabilities
    pz = sum(sum(pxyz, 1), 2); % Marginal probability density of Z
    pxy = sum(pxyz, 3); % Joint probability density of X and Y

    % Initialize mutual information
    MIZ = 0;

    % Loop through each value of x, y, and z
    for i = 1:size(pxyz, 1)
        for j = 1:size(pxyz, 2)
            for k = 1:size(pxyz, 3)
                if pxyz(i,j,k) ~= 0 && pz(k) ~= 0 && pxy(i,j) ~= 0
                    MIZ = MIZ + pxyz(i,j,k) * log2(pxyz(i,j,k) / (pz(k) * pxy(i,j)));
                end
            end
        end
    end
end