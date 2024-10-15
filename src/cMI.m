function cMI_values =  cMI(inputs, varargin)
%%% *function cMI_values =  cMI(inputs, outputs, opts)*
%%%
%%% This function calculates conditional mutual information (cMI) and
%%% other related measures based on the provided inputs, and requested outputs.
%%%
%%% Inputs:
%%%   - inputs: A cell array containing the input time series data. Each cell
%%%             represents a time series, where:
%%%             - inputs{1}: First time series (A) with dimensions
%%%                          nDims X [nTimepoints X] nTrials
%%%             - inputs{2}: Second time series (B) with dimensions
%%%                          nDims X [nTimepoints X] nTrials
%%%             - inputs{3}: Third time series (C) with dimensions
%%%                          nDims X [nTimepoints X] nTrials
%%%
%%%   - outputs: A cell array of strings specifying which entropies or cMI measures to compute.
%%%              Possible outputs include:
%%%               - 'I(A;B|C)'    : Conditional Mutual Information I(A;B|C)
%%%               - 'Ish(A;B|C)'  : Shuffled Conditional MI Ish(A;B|C)
%%%               - 'Ilin(A;B|C)' : Linear Conditional MI Ilin(A;B|C)
%%%
%%%   - varargin: Optional arguments, passed as a structure. Fields may include:
%%%              - bias:           Specifies the bias correction method to be used.
%%%                                Possible values include:
%%%                                'naive'                      :(default) - No correction applied.
%%%                                'qe', 'le'                   :quadratic/linear extrapolation.
%%%                                'ShuffSub'                   :Shuffle Subtraction (requires shuff parameter).
%%%                                'qe_ShuffSub', 'le_ShuffSub' :Combination of qe/le and ShuffSub.
%%%                                'pt'                         :Panzeri-Treves bias correction.
%%%
%%%              - bin_method:     Cell array specifying the binning method.
%%%                                It can have one or two entries:
%%%                                Possible values include:
%%%                                'eqpop'     : Equal population binning.
%%%                                'eqspace'   : Equal space binning.
%%%                                'threshold' : Binning based on a specified threshold.
%%%                                Default is {'none'}.
%%%
%%%              - n_bins:         Specifies the number of bins to use for binning.
%%%                                It can be a single integer or a cell array with one or two entries.
%%%                                Default number of bins is {3}.
%%%
%%%              - suppressWarnings: Boolean (true/false) to suppress warning messages.
%%%                                  Default is false.
%%%
%%% Outputs:
%%%   - cMI_values: A cell array containing the computed conditional MI values for
%%%                 each requested output type.
%%%
%%% Example:
%%% Suppose we have two time series of neural activity X1 and X2 and a
%%% time series Y. To compute the conditional mutual information I(X1;X2|Y) and the
%%% shuffled version, the function can be called as:
%%%
%%%   cMI_values = cMI({X1, X2, Y}, {'I(A;B|C)', 'Ish(A;B|C)'}, opts);
%%%
%%% Here, 'opts' represents additional options you may want to include, such as
%%% specifying the bias correction method, number of bins (n_bins), and other parameters.

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
        [outputs, opts] = check_inputs('cMI',inputs,varargin{:});
    end
else
    [outputs, opts] = check_inputs('cMI',inputs,varargin{:});
end

possibleOutputs = {'I(A;B|C)', 'Ish(A;B|C)', 'Ilin(A;B|C)'};
[isMember, indices] = ismember(outputs, possibleOutputs);
if any(~isMember)
    nonMembers = outputs(~isMember);
    msg = sprintf('Invalid Outputs: %s', strjoin(nonMembers, ', '));
    error('cMI:invalidOutput', msg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Step 2: Compute required Entropies for the requested Outputs           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
entropy_dependencies = struct( ...
    'I_AB_C', {{'H(A|C)', 'H(A|B,C)'}} ...
    );

required_entropies = {};
for ind = 1:length(indices)
    idx = indices(ind);
    switch possibleOutputs{idx}
        case 'I(A;B|C)'
            required_entropies = [required_entropies, entropy_dependencies.I_AB_C{:}];
        % case 'Ish(A;B|C)'
        %     required_entropies = [required_entropies, entropy_dependencies.Ish_AB_C{:}];
        % case 'Ilin(A;B|C)'
        %     required_entropies = [required_entropies, entropy_dependencies.Ilin_AB_C{:}];
    end
end
required_entropies = unique(required_entropies);
H_values = cell(1, length(required_entropies));
H_naive = cell(1, length(required_entropies));
H_shuff_all = cell(1, length(required_entropies));

opts_entropy = opts;
opts_entropy.compute_nulldist = false;
for i = 1:length(required_entropies)
    switch required_entropies{i}
        case 'H(A|C)'
            [H_values{i}, H_naive{i}, H_shuff_all{i}] = H({inputs{1}, inputs{3}}, {'H(A|B)'}, opts_entropy);
        % case 'H(B|C)'
        %     [H_values{i}, H_naive{i}, H_shuff_all{i}] = H({inputs{2}, inputs{3}}, {'H(A|B)'}, opts_entropy);
        % case 'H(A,B|C)'
        %     [H_values{i}, H_naive{i}, H_shuff_all{i}] = H({cat(1, inputs{1}, inputs{2}), inputs{3}}, {'H(A|B)'}, opts_entropy);
        % case 'Hsh(A,B|C)'
        %     [H_values{i}, H_naive{i}, H_shuff_all{i}] = H({cat(1, inputs{1}, inputs{2}), inputs{3}}, {'Hsh(A|B)'}, opts_entropy);
        case 'H(A|B,C)'
            [H_values{i}, H_naive{i}, H_shuff_all{i}] = H({inputs{1}, cat(1,inputs{2}, inputs{3})}, {'H(A|B)'}, opts_entropy);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Step 3: Compute requested Output Values                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize cell for MI values
cMI_values = cell(1, length(outputs));
cMI_naive = cell(1, length(outputs));
cMI_shuff_all = cell(1, length(outputs));

for i = 1:length(indices)
    idx = indices(i);
    nOut = nargout;
    switch possibleOutputs{idx}
        case 'I(A;B|C)'
            % I(A;B|C) = H(A|C) - H(A|B,C)
            %____________________________________________________________________________________________%
            if strcmp(opts.bias, 'shuffSub')
                H_AC_naive = H_naiv{strcmp(required_entropies, 'H(A|C)')};
                H_ABC_naive = H_naiv{strcmp(required_entropies, 'H(A|B,C)')};
                cMI_naive{i} =  H_AC_naive - H_ABC_naive;
                H_AC_shuff_all = cell2mat(H_shuff_all{strcmp(required_entropies, 'H(A|C)')});
                H_ABC_shuff_all = cell2mat(H_shuff_all{strcmp(required_entropies, 'H(A|B,C)')});
                for shuffIdx = 1:opts.shuff
                    H_AC_shuff = H_AC_shuff_all(shuffIdx);
                    H_ABC_shuff = H_ABC_shuff_all(shuffIdx);
                    cMI_shuff_all{i} =  H_AC_shuff - H_ABC_shuff;
                end
                cMI_values{i} = cMI_naive{i} - mean(cMI_shuff_all{i});
                nOut = 1;
            else
                H_AC = H_values{strcmp(required_entropies, 'H(A|C)')};
                H_ABC = H_values{strcmp(required_entropies, 'H(A|B,C)')};
                cMI_values{i} = H_AC -  H_ABC;
            end
            %____________________________________________________________________________________________%
            if nOut > 1
                H_AC_naive = H_naiv{strcmp(required_entropies, 'H(A|C)')};
                H_ABC_naive = H_naiv{strcmp(required_entropies, 'H(A|B,C)')};
                cMI_naive{i} =  H_AC_naive - H_ABC_naive;
                H_AC_shuff_all = cell2mat(H_shuff_all{strcmp(required_entropies, 'H(A|C)')});
                H_ABC_shuff_all = cell2mat(H_shuff_all{strcmp(required_entropies, 'H(A|B,C)')});
                if nOut > 2 && opts.shuff > 0
                    for shuffIdx = 1:opts.shuff
                        H_AC_shuff = H_AC_shuff_all(shuffIdx);
                        H_ABC_shuff = H_ABC_shuff_all(shuffIdx);
                        cMI_shuff_all{i} =  H_AC_shuff - H_ABC_shuff;
                    end
                end
            end
    end
end
end

