function [MI_values, MI_naive, MI_nullDist] = MI(inputs, varargin)
%%% MI - Calculate Mutual Information (MI) and related information-theoretic quantities
%%%
%%% This function calculates mutual information (MI) and other related
%%% measures based on the provided inputs, outputs, and optional parameters.
%%% It computes various forms of MI, such as I(A;B), Ish(A;B), Ilin(A;B),
%%% conditional information measures, and shuffled control variants.
%%%
%%% Inputs:
%%%   - inputs: A cell array containing the input time series data. Each cell represents a time series, where:
%%%             - inputs{1}: First time series (A) with dimensions
%%%                          nDims X nTimepoints X nTrials
%%%             - inputs{2}: Second time series (B) with dimensions
%%%                          nDims X nTimepoints X nTrials
%%%   - outputs: A cell array of strings specifying which entropies to compute.
%%%              Possible outputs include:
%%%               - 'I(A;B)'    : Mutual Information I(A;B)
%%%               - 'Ilin(A;B)' : Linear MI Ilin(A;B)
%%%               - 'coI(A;B)'  : Co-information coI(A;B)
%%%               - 'Iss(A)'    : Sig. Sim Information Iss(A)
%%%               - 'Ic(A;B)'   : Sum of Ici and Icd - Ic(A;B)
%%%               - 'Ici(A;B)'  : Correlation Independent Information Ici(A;B)
%%%               - 'Icd(A;B)'  : Correlation Dependent Information Icd(A;B)
%%%
%%%   - varargin: Optional arguments, passed as a structure. Fields may include:
%%%              - bias:           Specifies the bias correction method to be used.
%%%                                Possible values include:
%%%                                'naive'                      :(default) - No correction applied.
%%%                                'qe', 'le'                   :quadratic/linear extrapolation (need to specify xtrp as number of extrapolations).
%%%                                'ShuffSub'                   :Shuffle Substraction (need to specify shuff as number of shufflings).
%%%                                'qe_ShuffSub', 'le_ShuffSub' :Combination of qe/le and Shuffsub (need to specify shuff and xtrp).
%%%                                'pt'                         :Panzeri-Treves bias correction.
%%%                                'bub'                        :
%%%                                'shuffCorr'                  :
%%%                                Users can also define their own custom bias correction method
%%%                                (see help for correction.m)
%%%
%%%              - bin_method:     Cell array specifying the binning method.
%%%                                It can have one or two entries:
%%%                                If one entry is provided, it will be applied to both A and B.
%%%                                Possible values include:
%%%                                'eqpop'     : Equal population binning.
%%%                                'eqspace'   : Equal space binning.
%%%                                'threshold' : Binning based on a specified threshold.
%%%                                Users can also define their own custom binning method
%%%                                (see help for binning.m).
%%%                                Default is {'none'}.
%%%
%%%              - n_bins:         Specifies the number of bins to use for binning.
%%%                                It can be a single integer or a cell array with one or two entries.
%%%                                If one entry is provided, it will be used for both A and B.
%%%                                This integer defines how the continuous values will be
%%%                                discretized into bins for analysis.
%%%                                Default number of bins is {3}.
%%%
%%%              - suppressWarnings: Boolean (true/false) to suppress warning messages.
%%%                                  Default is false, meaning warnings will be shown.
%%% Outputs:
%%%   - MI_values: A cell array containing the computed MI values as specified in the outputs argument.
%%%   - MI_naive: A cell array containing the naive MI estimates.
%%%   - MI_shuff_all: A value indicating the all results of the shuffling procedure (0 if not performed).
%%%
%%% EXAMPLE
%%% Suppose we have two time series of groups of neurons X1 and X2
%%% and two time series of groups of neurons Y1 and Y2.
%%% (Structur of X1, X2, Y1, Y1 is nNeurons x nTrials)
%%%
%%% We can structure our inputs as follows:
%%% Thus, the total input for A and B would be:
%%% A = cat(1, X1, X2);  % Concatenates X1 and X2 along the first dimension (neurons)
%%% B = cat(1, Y1, Y2);  % Concatenates Y1 and Y2 along the first dimension (neurons)
%%%
%%% To compute the Mutual Information and the Shuffled MI from A (X) to B (Y), the function can be called as:
%%% MI_values = MI({A, B}, {'I(A;B)', Ish(A;B}, opts);
%%%
%%% Here, 'opts' represents additional options you may want to include, such as
%%% specifying the bias correction method, number of bins (n_bins), and other parameters as needed.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Check Inputs, Check OutputList, Fill missing opts with default values %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(varargin) > 1
    opts = varargin{2};
    if isfield(opts, 'isChecked')
        if opts.isChecked
            outputs = varargin{1};
        end
    else
        [outputs, opts] = check_inputs('MI',inputs,varargin{:});
    end
else
    [outputs, opts] = check_inputs('MI',inputs,varargin{:});
end
possibleOutputs = {'I(A;B)', 'Ish(A;B)',  'Ilin(A;B)', 'coI(A;B)', 'coIsh(A;B)', ...
    'Iss(A)', 'Ic(A;B)', 'Icsh(A;B)', 'Ici(A;B)', 'Icd(A;B)', 'Icdsh(A;B)'};
[isMember, indices] = ismember(outputs, possibleOutputs);
if any(~isMember)
    nonMembers = outputs(~isMember);
    msg = sprintf('Invalid Outputs: %s', strjoin(nonMembers, ', '));
    error('MI:invalidOutput', msg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Step 2: Compute required Entropies for the requested Outputs           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
entropy_dependencies = struct( ...
    'I_A_B', {{'H(A)', 'H(A|B)'}}, ...
    'Ish_A_B', {{'H(A)', 'H(A|B)', 'Hsh(A|B)', 'Hind(A|B)'}}, ...
    'Ilin_A_B', {{'Hlin(A)', 'H(A|B)', 'Hind(A|B)'}}, ...
    'coI_A_B', {{'H(A)', 'H(A|B)', 'Hlin(A)', 'Hind(A|B)'}}, ...
    'coIsh_A_B', {{'H(A)', 'H(A|B)', 'Hsh(A|B)', 'Hlin(A)'}}, ...
    'Iss_A', {{'Hind(A)', 'Hlin(A)'}}, ...
    'Ic_A_B', {{'H(A)', 'H(A|B)', 'Hind(A|B)', 'Hind(A)'}}, ...
    'Icsh_A_B', {{'H(A)', 'H(A|B)', 'Hsh(A|B)', 'Hind(A)'}}, ...
    'Ici_A_B', {{'Chi(A)', 'Hind(A)'}}, ...
    'Icd_A_B', {{'H(A)', 'H(A|B)', 'Chi(A)', 'Hind(A|B)'}}, ...
    'Icdsh_A_B', {{'H(A)', 'H(A|B)', 'Hsh(A|B)', 'Chi(A)'}} ...
    );
required_entropies = {};
for ind = 1:length(indices)
    idx = indices(ind);
    switch possibleOutputs{idx}
        case 'I(A;B)'
            required_entropies = [required_entropies, entropy_dependencies.I_A_B{:}];
        case 'Ish(A;B)'
            required_entropies = [required_entropies, entropy_dependencies.Ish_A_B{:}];
        case 'Ilin(A;B)'
            required_entropies = [required_entropies, entropy_dependencies.Ilin_A_B{:}];
        case 'coI(A;B)'
            required_entropies = [required_entropies, entropy_dependencies.coI_A_B{:}];
        case 'coIsh(A;B)'
            required_entropies = [required_entropies, entropy_dependencies.coIsh_A_B{:}];
        case 'Iss(A)'
            required_entropies = [required_entropies, entropy_dependencies.Iss_A{:}];
        case 'Ic(A;B)'
            required_entropies = [required_entropies, entropy_dependencies.Ic_A_B{:}];
        case 'Icsh(A;B)'
            required_entropies = [required_entropies, entropy_dependencies.Icsh_A_B{:}];
        case 'Ici(A;B)'
            required_entropies = [required_entropies, entropy_dependencies.Ici_A_B{:}];
        case 'Icd(A;B)'
            required_entropies = [required_entropies, entropy_dependencies.Icd_A_B{:}];
        case 'Icdsh(A;B)'
            required_entropies = [required_entropies, entropy_dependencies.Icdsh_A_B{:}];
    end
end
required_entropies = unique(required_entropies);
corr = opts.bias;
opts_entropies = opts;
opts_entropies.compute_nulldist = false;
opts_entropies.isBinned = false;
if strcmp(opts.bias, 'shuffCorr')
    [MI_values, MI_naive, MI_nullDist] = correction(inputs, outputs, corr, @MI, opts_entropies);
    return
end
[H_values, H_naive, H_shuff_all] = H(inputs, required_entropies, opts_entropies);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Step 3: Compute requested Output Values                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize cell for MI values
MI_values = cell(1, length(outputs));
MI_naive = cell(1, length(outputs));
MI_shuff_all = cell(1, length(outputs));

for i = 1:length(indices)
    idx = indices(i);
    nOut = nargout;
    switch possibleOutputs{idx}
        case 'I(A;B)'
            % I(A;B) = H(A) - H(A|B)
            %____________________________________________________________________________________________%
            if strcmp(opts.bias, 'shuffSub')
                H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')};
                H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')};
                MI_naive{i} = H_A_naive - H_A_B_naive;
                for shuffIdx = 1:opts.shuff
                    H_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'H(A|B)')}(shuffIdx);
                    MI_shuff_all{i}(shuffIdx) = H_A_naive - H_A_B_shuff;
                end
                MI_values{i} = MI_naive{i} - mean(MI_shuff_all{i});
                nOut = 1;
            else
                H_A = H_values{strcmp(required_entropies, 'H(A)')};
                H_A_B = H_values{strcmp(required_entropies, 'H(A|B)')};
                MI_values{i} = H_A - H_A_B;
            end
            %____________________________________________________________________________________________%
            if nOut > 1
                H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')};
                H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')};
                MI_naive{i} = H_A_naive - H_A_B_naive;
                if nOut > 2 && opts.shuff > 0
                    for shuffIdx = 1:opts.shuff
                        H_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'H(A|B)')}(shuffIdx);
                        MI_shuff_all{i}(shuffIdx) = H_A_naive - H_A_B_shuff;
                    end
                end
            end

        case 'Ish(A;B)'
            % Ish(A;B) = H(A) - Hind(A|B) + Hsh(A|B) - H(A|B)
            %____________________________________________________________________________________________%
            if strcmp(opts.bias, 'shuffSub')
                H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')};
                Hind_A_B_naive = H_naive{strcmp(required_entropies, 'Hind(A|B)')};
                Hsh_A_B_naive = H_naive{strcmp(required_entropies, 'Hsh(A|B)')};
                H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')};
                MI_naive{i} = H_A_naive - Hind_A_B_naive + Hsh_A_B_naive - H_A_B_naive;
                for shuffIdx = 1:opts.shuff
                    Hind_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'Hind(A|B)')}(shuffIdx);
                    Hsh_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'Hsh(A|B)')}(shuffIdx);
                    H_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'H(A|B)')}(shuffIdx);
                    MI_shuff_all{i}(shuffIdx) = H_A_naive - Hind_A_B_shuff + Hsh_A_B_shuff - H_A_B_shuff;
                end
                MI_values{i} = MI_naive{i} - mean(MI_shuff_all{i});
                nOut = 1;
            else
                H_A = H_values{strcmp(required_entropies, 'H(A)')};
                Hind_A_B = H_values{strcmp(required_entropies, 'Hind(A|B)')};
                Hsh_A_B = H_values{strcmp(required_entropies, 'Hsh(A|B)')};
                H_A_B = H_values{strcmp(required_entropies, 'H(A|B)')};
                MI_values{i} = H_A - Hind_A_B + Hsh_A_B - H_A_B;
            end
            %____________________________________________________________________________________________%
            if nOut > 1
                H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')};
                Hind_A_B_naive = H_naive{strcmp(required_entropies, 'Hind(A|B)')};
                Hsh_A_B_naive = H_naive{strcmp(required_entropies, 'Hsh(A|B)')};
                H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')};
                MI_naive{i} = H_A_naive - Hind_A_B_naive + Hsh_A_B_naive - H_A_B_naive;
                if nOut > 2 && opts.shuff > 0
                    for shuffIdx = 1:opts.shuff
                        Hind_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'Hind(A|B)')}(shuffIdx);
                        Hsh_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'Hsh(A|B)')}(shuffIdx);
                        H_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'H(A|B)')}(shuffIdx);
                        MI_shuff_all{i}(shuffIdx) = H_A_naive - Hind_A_B_shuff + Hsh_A_B_shuff - H_A_B_shuff;
                    end
                end
            end

        case 'Ilin(A;B)'
            % I_lin(A;B) = Hlin(A) - Hind(A|B)
            %____________________________________________________________________________________________%
            if strcmp(opts.bias, 'shuffSub')
                Hlin_A_naive = H_naive{strcmp(required_entropies, 'Hlin(A)')};
                Hind_A_B_naive = H_naive{strcmp(required_entropies, 'Hind(A|B)')};
                MI_naive{i} = Hlin_A_naive - Hind_A_B_naive;
                for shuffIdx = 1:opts.shuff
                    Hind_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'Hind(A|B)')}(shuffIdx);
                    MI_shuff_all{i}(shuffIdx) = Hlin_A_naive - Hind_A_B_shuff;
                end
                MI_values{i} = MI_naive{i} - mean(MI_shuff_all{i});
                nOut = 1;
            else
                Hlin_A = H_values{strcmp(required_entropies, 'Hlin(A)')};
                Hind_A_B = H_values{strcmp(required_entropies, 'Hind(A|B)')};
                MI_values{i} = Hlin_A - Hind_A_B;
            end
            %____________________________________________________________________________________________%
            if nOut > 1
                Hlin_A_naive = H_naive{strcmp(required_entropies, 'Hlin(A)')};
                Hind_A_B_naive = H_naive{strcmp(required_entropies, 'Hind(A|B)')};
                MI_naive{i} = Hlin_A_naive - Hind_A_B_naive;
                if nOut > 2 && opts.shuff > 0
                    for shuffIdx = 1:opts.shuff
                        Hind_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'Hind(A|B)')}(shuffIdx);
                        MI_shuff_all{i}(shuffIdx) = Hlin_A_naive - Hind_A_B_shuff;
                    end
                end
            end

        case 'coI(A;B)'
            % coI(A;B) = H(A) - H(A|B) - Hlin(A) + Hind(A|B)
            %____________________________________________________________________________________________%
            if strcmp(opts.bias, 'shuffSub')
                H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')};
                H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')};
                Hlin_A_naive = H_naive{strcmp(required_entropies, 'Hlin(A)')};
                Hind_A_B_naive = H_naive{strcmp(required_entropies, 'Hind(A|B)')};
                MI_naive{i} = H_A_naive - H_A_B_naive - Hlin_A_naive + Hind_A_B_naive;
                for shuffIdx = 1:opts.shuff
                    H_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'H(A|B)')}(shuffIdx);
                    Hind_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'Hind(A|B)')}(shuffIdx);
                    MI_shuff_all{i}(shuffIdx) = H_A_naive - H_A_B_shuff - Hlin_A_naive + Hind_A_B_shuff;
                end
                MI_values{i} = MI_naive{i} - mean(MI_shuff_all{i});
                nOut = 1;
            else
                H_A = H_values{strcmp(required_entropies, 'H(A)')};
                Hlin_A = H_values{strcmp(required_entropies, 'Hlin(A)')};
                H_A_B = H_values{strcmp(required_entropies, 'H(A|B)')};
                Hind_A_B = H_values{strcmp(required_entropies, 'Hind(A|B)')};
                MI_values{i} = H_A - H_A_B - Hlin_A + Hind_A_B;
            end
            %____________________________________________________________________________________________%
            if nOut > 1
                H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')};
                H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')};
                Hlin_A_naive = H_naive{strcmp(required_entropies, 'Hlin(A)')};
                Hind_A_B_naive = H_naive{strcmp(required_entropies, 'Hind(A|B)')};
                MI_naive{i} = H_A_naive - H_A_B_naive - Hlin_A_naive + Hind_A_B_naive;
                if nOut > 2 && opts.shuff > 0
                    for shuffIdx = 1:opts.shuff
                        H_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'H(A|B)')}(shuffIdx);
                        Hind_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'Hind(A|B)')}(shuffIdx);
                        MI_shuff_all{i}(shuffIdx) = H_A_naive - H_A_B_shuff - Hlin_A_naive + Hind_A_B_shuff;
                    end
                end
            end
        case 'coIsh(A;B)'
            % coIsh(A;B) = H(A) + Hsh(A|B) - H(A|B) - Hlin(A)
            %____________________________________________________________________________________________%
            if strcmp(opts.bias, 'shuffSub')
                H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')};
                Hsh_A_B_naive = H_naive{strcmp(required_entropies, 'Hsh(A|B)')};
                H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')};
                Hlin_A_naive = H_naive{strcmp(required_entropies, 'Hlin(A)')};
                MI_naive{i} = H_A_naive + Hsh_A_B_naive - H_A_B_naive - Hlin_A_naive;
                for shuffIdx = 1:opts.shuff
                    Hsh_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'Hsh(A|B)')}(shuffIdx);
                    H_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'H(A|B)')}(shuffIdx);
                    MI_shuff_all{i}(shuffIdx) = H_A_naive + Hsh_A_B_shuff - H_A_B_shuff - Hlin_A_naive;
                end
                MI_values{i} = MI_naive{i} - mean(MI_shuff_all{i});
                nOut = 1;
            else
                H_A = H_values{strcmp(required_entropies, 'H(A)')};
                Hlin_A = H_values{strcmp(required_entropies, 'Hlin(A)')};
                Hsh_A_B = H_values{strcmp(required_entropies, 'Hsh(A|B)')};
                H_A_B = H_values{strcmp(required_entropies, 'H(A|B)')};
                MI_values{i} = H_A + Hsh_A_B - H_A_B - Hlin_A;
            end
            %____________________________________________________________________________________________%
            if nOut > 1
                H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')};
                Hsh_A_B_naive = H_naive{strcmp(required_entropies, 'Hsh(A|B)')};
                H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')};
                Hlin_A_naive = H_naive{strcmp(required_entropies, 'Hlin(A)')};
                MI_naive{i} = H_A_naive + Hsh_A_B_naive - H_A_B_naive - Hlin_A_naive;
                if nOut > 2 && opts.shuff > 0
                    for shuffIdx = 1:opts.shuff
                        Hsh_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'Hsh(A|B)')}(shuffIdx);
                        H_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'H(A|B)')}(shuffIdx);
                        MI_shuff_all{i}(shuffIdx) = H_A_naive + Hsh_A_B_shuff - H_A_B_shuff - Hlin_A_naive;
                    end
                end
            end
        case 'Iss(A)'
            % Iss(A) = Hind(A) - Hlin(A)
            %____________________________________________________________________________________________%
            if strcmp(opts.bias, 'shuffSub')
                Hind_A_naive = H_naive{strcmp(required_entropies, 'Hind(A)')};
                Hlin_A_naive = H_naive{strcmp(required_entropies, 'Hlin(A)')};
                MI_naive{i} = Hind_A_naive - Hlin_A_naive;
                for shuffIdx = 1:opts.shuff
                    Hind_A_shuff= H_shuff_all{strcmp(required_entropies, 'Hind(A)')}(shuffIdx);
                    MI_shuff_all{i}(shuffIdx) = Hind_A_shuff - Hlin_A_naive;
                end
                MI_values{i} = MI_naive{i} - mean(MI_shuff_all{i});
                nOut = 1;
            else
                Hind_A = H_values{strcmp(required_entropies, 'Hind(A)')};
                Hlin_A = H_values{strcmp(required_entropies, 'Hlin(A)')};
                MI_values{i} = Hind_A - Hlin_A;
            end
            %____________________________________________________________________________________________%
            if nOut > 1
                Hind_A_naive = H_naive{strcmp(required_entropies, 'Hind(A)')};
                Hlin_A_naive = H_naive{strcmp(required_entropies, 'Hlin(A)')};
                MI_naive{i} = Hind_A_naive - Hlin_A_naive;
                if nOut > 2 && opts.shuff > 0
                    for shuffIdx = 1:opts.shuff
                        Hind_A_shuff= H_shuff_all{strcmp(required_entropies, 'Hind(A)')}(shuffIdx);
                        MI_shuff_all{i}(shuffIdx) = Hind_A_shuff - Hlin_A_naive;
                    end
                end
            end
        case 'Ic(A;B)'% possibleOutputs = { 'H(A|B)', 'Hind(A|B)', 'Hsh(A|B)'};
            % Ic(A;B) = H(A) - H(A|B) + Hind(A|B) - Hind(A)
            %____________________________________________________________________________________________%
            if strcmp(opts.bias, 'shuffSub')
                H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')};
                H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')};
                Hind_A_B_naive = H_naive{strcmp(required_entropies, 'Hind(A|B)')};
                Hind_A_naive = H_naive{strcmp(required_entropies, 'Hind(A)')};
                MI_naive{i} = H_A_naive - H_A_B_naive + Hind_A_B_naive - Hind_A_naive;
                for shuffIdx = 1:opts.shuff
                    Hind_A_shuff = H_shuff_all{strcmp(required_entropies, 'Hind(A)')}(shuffIdx);
                    H_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'H(A|B)')}(shuffIdx);
                    Hind_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'Hind(A|B)')}(shuffIdx);
                    MI_shuff_all{i}(shuffIdx) = H_A_naive - H_A_B_shuff + Hind_A_B_shuff - Hind_A_shuff;
                end
                MI_values{i} = MI_naive{i} - mean(MI_shuff_all{i});
                nOut = 1;
            else
                H_A = H_values{strcmp(required_entropies, 'H(A)')};
                H_A_B = H_values{strcmp(required_entropies, 'H(A|B)')};
                Hind_A_B = H_values{strcmp(required_entropies, 'Hind(A|B)')};
                Hind_A = H_values{strcmp(required_entropies, 'Hind(A)')};
                MI_values{i} = H_A - H_A_B + Hind_A_B - Hind_A;
            end
            %____________________________________________________________________________________________%
            if nOut > 1
                H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')};
                H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')};
                Hind_A_B_naive = H_naive{strcmp(required_entropies, 'Hind(A|B)')};
                Hind_A_naive = H_naive{strcmp(required_entropies, 'Hind(A)')};
                MI_naive{i} = H_A_naive - H_A_B_naive + Hind_A_B_naive - Hind_A_naive;
                if nOut > 2 && opts.shuff > 0
                    for shuffIdx = 1:opts.shuff
                        Hind_A_shuff = H_shuff_all{strcmp(required_entropies, 'Hind(A)')}(shuffIdx);
                        H_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'H(A|B)')}(shuffIdx);
                        Hind_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'Hind(A|B)')}(shuffIdx);
                        MI_shuff_all{i}(shuffIdx) = H_A_naive - H_A_B_shuff + Hind_A_B_shuff - Hind_A_shuff;
                    end
                end
            end
        case 'Icsh(A;B)'
            % Icsh(A;B) = H(A) + Hsh(A|B) - H(A|B) - Hind(A)
            %____________________________________________________________________________________________%
            if strcmp(opts.bias, 'shuffSub')
                H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')};
                Hsh_A_B_naive = H_naive{strcmp(required_entropies, 'Hsh(A|B)')};
                H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')};
                Hind_A_naive = H_naive{strcmp(required_entropies, 'Hind(A)')};
                MI_naive{i} = H_A_naive + Hsh_A_B_naive - H_A_B_naive - Hind_A_naive;
                for shuffIdx = 1:opts.shuff
                    Hind_A_shuff = H_shuff_all{strcmp(required_entropies, 'Hind(A)')}(shuffIdx);
                    Hsh_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'Hsh(A|B)')}(shuffIdx);
                    H_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'H(A|B)')}(shuffIdx);
                    MI_shuff_all{i}(shuffIdx) = H_A_naive + Hsh_A_B_shuff - H_A_B_shuff - Hind_A_shuff;
                end
                MI_values{i} = MI_naive{i} - mean(MI_shuff_all{i});
                nOut = 1;
            else
                H_A = H_values{strcmp(required_entropies, 'H(A)')};
                Hsh_A_B = H_values{strcmp(required_entropies, 'Hsh(A|B)')};
                H_A_B = H_values{strcmp(required_entropies, 'H(A|B)')};
                Hind_A = H_values{strcmp(required_entropies, 'Hind(A)')};
                MI_values{i} = H_A + Hsh_A_B - H_A_B - Hind_A;
            end
            %____________________________________________________________________________________________%
            if nOut > 1
                H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')};
                Hsh_A_B_naive = H_naive{strcmp(required_entropies, 'Hsh(A|B)')};
                H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')};
                Hind_A_naive = H_naive{strcmp(required_entropies, 'Hind(A)')};
                MI_naive{i} = H_A_naive + Hsh_A_B_naive - H_A_B_naive - Hind_A_naive;
                if nOut > 2 && opts.shuff > 0
                    for shuffIdx = 1:opts.shuff
                        Hind_A_shuff = H_shuff_all{strcmp(required_entropies, 'Hind(A)')}(shuffIdx);
                        Hsh_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'Hsh(A|B)')}(shuffIdx);
                        H_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'H(A|B)')}(shuffIdx);
                        MI_shuff_all{i}(shuffIdx) = H_A_naive + Hsh_A_B_shuff - H_A_B_shuff - Hind_A_shuff;
                    end
                end
            end
        case 'Ici(A;B)'
            % Ici(A;B) = Chi(A) - Hind(A)
            %____________________________________________________________________________________________%
            if strcmp(opts.bias, 'shuffSub')
                Chi_A_naive = H_naive{strcmp(required_entropies, 'Chi(A)')};
                Hind_A_naive = H_naive{strcmp(required_entropies, 'Hind(A)')};
                MI_naive{i} = Chi_A_naive - Hind_A_naive;
                for shuffIdx = 1:opts.shuff
                    Hind_A_shuff = H_shuff_all{strcmp(required_entropies, 'Hind(A)')}(shuffIdx);
                    Chi_A_shuff = H_shuff_all{strcmp(required_entropies, 'Chi(A)')}(shuffIdx);
                    MI_shuff_all{i}(shuffIdx) = Chi_A_shuff - Hind_A_shuff;
                end
                MI_values{i} = MI_naive{i} - mean(MI_shuff_all{i});
                nOut = 1;
            else
                Chi_A = H_values{strcmp(required_entropies, 'Chi(A)')};
                Hind_A = H_values{strcmp(required_entropies, 'Hind(A)')};
                MI_values{i} = Chi_A - Hind_A;
            end
            %____________________________________________________________________________________________%
            if nOut > 1
                Chi_A_naive = H_naive{strcmp(required_entropies, 'Chi(A)')};
                Hind_A_naive = H_naive{strcmp(required_entropies, 'Hind(A)')};
                MI_naive{i} = Chi_A_naive - Hind_A_naive;
                if nOut > 2 && opts.shuff > 0
                    for shuffIdx = 1:opts.shuff
                        Hind_A_shuff = H_shuff_all{strcmp(required_entropies, 'Hind(A)')}(shuffIdx);
                        Chi_A_shuff = H_shuff_all{strcmp(required_entropies, 'Chi(A)')}(shuffIdx);
                        MI_shuff_all{i}(shuffIdx) = Chi_A_shuff - Hind_A_shuff;
                    end
                end
            end
        case 'Icd(A;B)'
            % Icd(A;B) = H(A) - H(A|B) - Chi(A) + Hind(A|B)
            %____________________________________________________________________________________________%
            if strcmp(opts.bias, 'shuffSub')
                H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')};
                H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')};
                Chi_A_naive = H_naive{strcmp(required_entropies, 'Chi(A)')};
                Hind_A_B_naive = H_naive{strcmp(required_entropies, 'Hind(A|B)')};
                MI_naive{i} = H_A_naive - H_A_B_naive - Chi_A_naive + Hind_A_B_naive;
                for shuffIdx = 1:opts.shuff
                    Chi_A_shuff = H_shuff_all{strcmp(required_entropies, 'Chi(A)')}(shuffIdx);    
                    H_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'H(A|B)')}(shuffIdx);
                    Hind_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'Hind(A|B)')}(shuffIdx);
                    MI_shuff_all{i}(shuffIdx)= H_A_naive - H_A_B_shuff - Chi_A_shuff+ Hind_A_B_shuff;
                end
                MI_values{i} = MI_naive{i} - mean(MI_shuff_all{i});
                nOut = 1;
            else
                H_A = H_values{strcmp(required_entropies, 'H(A)')};
                H_A_B = H_values{strcmp(required_entropies, 'H(A|B)')};
                Chi_A = H_values{strcmp(required_entropies, 'Chi(A)')};
                Hind_A_B = H_values{strcmp(required_entropies, 'Hind(A|B)')};
                MI_values{i} = H_A - H_A_B - Chi_A + Hind_A_B;
            end
            %____________________________________________________________________________________________%
            if nOut > 1
                H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')};
                H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')};
                Chi_A_naive = H_naive{strcmp(required_entropies, 'Chi(A)')};
                Hind_A_B_naive = H_naive{strcmp(required_entropies, 'Hind(A|B)')};
                MI_naive{i} = H_A_naive - H_A_B_naive - Chi_A_naive + Hind_A_B_naive;
                if nOut > 2 && opts.shuff > 0
                    for shuffIdx = 1:opts.shuff
                        Chi_A_shuff = H_shuff_all{strcmp(required_entropies, 'Chi(A)')}(shuffIdx);
                        H_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'H(A|B)')}(shuffIdx);
                        Hind_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'Hind(A|B)')}(shuffIdx);
                        MI_shuff_all{i}(shuffIdx)= H_A_naive - H_A_B_shuff - Chi_A_shuff+ Hind_A_B_shuff;
                    end
                end
            end
        case 'Icdsh(A;B)'
            % Icdsh(A;B) = H(A) + Hsh(A|B) - H(A|B) - Chi(A)
            %____________________________________________________________________________________________%
            if strcmp(opts.bias, 'shuffSub')
                H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')};
                Hsh_A_B_naive = H_naive{strcmp(required_entropies, 'Hsh(A|B)')};
                H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')};
                Chi_A_naive = H_naive{strcmp(required_entropies, 'Chi(A)')};
                MI_naive{i} = H_A_naive + Hsh_A_B_naive - H_A_B_naive - Chi_A_naive;
                for shuffIdx = 1:opts.shuff
                    Hsh_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'Hsh(A|B)')}(shuffIdx);
                    H_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'H(A|B)')}(shuffIdx);
                    Chi_A_shuff = H_shuff_all{strcmp(required_entropies, 'Chi(A)')}(shuffIdx);
                    MI_shuff_all{i}(shuffIdx) = H_A_naive + Hsh_A_B_shuff - H_A_B_shuff - Chi_A_shuff;
                end
                MI_values{i} = MI_naive{i} - mean(MI_shuff_all{i});
                nOut = 1;
            else
                H_A = H_values{strcmp(required_entropies, 'H(A)')};
                Hsh_A_B = H_values{strcmp(required_entropies, 'Hsh(A|B)')};
                H_A_B = H_values{strcmp(required_entropies, 'H(A|B)')};
                Chi_A = H_values{strcmp(required_entropies, 'Chi(A)')};
                MI_values{i} = H_A + Hsh_A_B - H_A_B - Chi_A;
            end
            %____________________________________________________________________________________________%
            if nOut > 1
                H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')};
                Hsh_A_B_naive = H_naive{strcmp(required_entropies, 'Hsh(A|B)')};
                H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')};
                Chi_A_naive = H_naive{strcmp(required_entropies, 'Chi(A)')};
                MI_naive{i} = H_A_naive + Hsh_A_B_naive - H_A_B_naive - Chi_A_naive;
                if nOut > 2 && opts.shuff > 0
                    for shuffIdx = 1:opts.shuff
                        Hsh_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'Hsh(A|B)')}(shuffIdx);
                        H_A_B_shuff = H_shuff_all{strcmp(required_entropies, 'H(A|B)')}(shuffIdx);
                        Chi_A_shuff = H_shuff_all{strcmp(required_entropies, 'Chi(A)')}(shuffIdx);
                        MI_shuff_all{i}(shuffIdx) = H_A_naive + Hsh_A_B_shuff - H_A_B_shuff - Chi_A_shuff;
                    end
                end
            end
    end

    if opts.computeNulldist
        nullDist_opts = opts;
        nullDist_opts.computeNulldist = false;
        MI_nullDist = create_NullDistribution(inputs, outputs, @MI, nullDist_opts);
    else
        MI_nullDist = MI_shuff_all;
    end 
end


