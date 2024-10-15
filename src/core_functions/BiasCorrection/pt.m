function [corrected_v, naive_v] = pt(inputs, outputs, corefunc, varargin)

% Check Outputslist
possibleOutputs = {'H(A)', 'H(A|B)', 'Hlin(A)', 'Hind(A)', 'Hind(A|B)', 'Chi(A)','Hsh(A)', 'Hsh(A|B)'};
[isMember, indices] = ismember(outputs, possibleOutputs);
if any(~isMember)
    nonMembers = outputs(~isMember);
    msg = sprintf('Invalid Outputs: %s', strjoin(nonMembers, ', '));
    error('H:invalidOutput', msg);
end

DimsA = size(inputs{1});
DimsB = size(inputs{2});
if DimsA(end) ~= DimsB(end)
    msg = sprintf('The number of trials for A (%d) and B (%d) are not consistent. Ensure both variables have the same number of trials.',DimsA(end),DimsB(end));
    error('H:InvalidInput', msg);
end
nTrials = DimsA(end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Step 2: Prepare Data (binning/reduce dimensions)                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts = varargin{1};
if ~opts.isbinned
    inputs_b = binning(inputs,opts);
    opts.isbinned = true;
    for c=1:length(inputs)
        inputs_b{c} = inputs_b{c}+1;
    end
else
    inputs_b = inputs;
end

inputs_1d = inputs_b;
if DimsA(1) > 1
    inputs_1d{1} = reduce_dim(inputs_b{1}, 1);
    if  any(strcmp(outputs,'Hlin(A)')) || any(strcmp(outputs,'Hind(A)')) || any(strcmp(outputs, 'Hind(A|B)'))
        inputs_1d{3} = inputs{1};
    end 
end
if DimsB(1) > 1
    inputs_1d{2} = reduce_dim(inputs_b{2}, 1);
end
if DimsA(2:end) ~= DimsB(2:end)
    msg = 'Inconsistent sizes of A and B';
    error('H:inconsistentSizes', msg);
end

naiveopts = opts;
naiveopts.bias = 'naive';

naive_v = H(inputs_1d, outputs, naiveopts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Step 3.B: Compute required Probability Distributions               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
entropy_distributions = struct( ...
    'H_A', {{'P(A)'}}, ...
    'H_A_B', {{'P(A|B)', 'P(B)'}}, ...
    'Hlin_A', {{'Plin(A)'}}, ...
    'Hsh_A', {{'Psh(A)'}}, ...
    'Hsh_A_B', {{'Psh(A|B)', 'P(B)'}} ...
    );

entropies_nullDist = 0;
required_distributions = {};
for ind = 1:length(indices)
    idx = indices(ind);
    switch possibleOutputs{idx}
        case 'H(A)'
            required_distributions = [required_distributions, entropy_distributions.H_A{:}];
        case 'H(A|B)'
            required_distributions = [required_distributions, entropy_distributions.H_A_B{:}];
        case 'Hlin(A)'
            required_distributions = [required_distributions, entropy_distributions.Hlin_A{:}];
        case 'Hsh(A)'
            required_distributions = [required_distributions, entropy_distributions.Hsh_A{:}];
        case 'Hsh(A|B)'
            required_distributions = [required_distributions, entropy_distributions.Hsh_A_B{:}];
    end
end
required_distributions = unique(required_distributions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Step 4.B: Compute requested Entropy Biases               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
entropy_biases = cell(1, length(outputs));
corrected_v = cell(1, length(outputs));
A = inputs{1};
redA = reduce_dim(A,1);
nA = length(unique(redA));
for i = 1:length(indices)
    idx = indices(i);
    switch possibleOutputs{idx}
        case 'H(A)'

            
            entropy_biases{i} = pt_core(redA, nA, nTrials);
        case 'H(A|B)'
            A = inputs{1};
            B = inputs{2};

            redA = reduce_dim(A,1);
            redB = reduce_dim(B,1);
            uB = unique(redB);
            nB = length(uB);
            bias_pt =0;
            for bi=1:nB
                redAs = redA(redB == uB(bi));
                nAs = length(unique(redAs));
                bias_pt =  bias_pt + pt_core(redAs, nAs, nTrials);
            end
            entropy_biases{i} = bias_pt;
        case 'Hlin(A)'
            A = inputs{1};
            bias_pt = 0;
            for row=1:size(A,1)
                Arow = A(row,:);
                nArow = length(unique(Arow));
                bias_pt =  bias_pt + pt_core(Ai, nArow, nTrials);
            end
            entropy_biases{i} = bias_pt;
        case 'Hsh(A)'
            A = inputs{1};
            redA = reduce_dim(A,1);
            nA = length(unique(redA));
            entropy_biases{i} = pt_core(redA, nA, length(redA));
        case 'Hsh(A|B)'
            A = inputs{1};
            B = inputs{2};

            redA = reduce_dim(A,1);
            redB = reduce_dim(B,1);
            uB = unique(redB);
            nB = length(uB);
            bias_pt =0;
            for bi=1:nB
                redAs = redA(redB == uB(bi));
                nAs = length(unique(redAs));
                bias_pt =  bias_pt + nAs;
            end
            entropy_biases{i} = bias_pt;
    end
    corrected_v{i} = naive_v{i} + entropy_biases{i};
end

end


function bias = pt_core(C_ptr, Rtot, N)
    % Function to compute the bias as described in Panzeri & Treves (1996)
    % Input:
    %   C_ptr - Array of values
    %   Rtot - Total number of bins
    %   N - Sample size
    % Output:
    %   bias - Computed bias value

    % Number of non-zero C values
    NnonZeroCvalues = sum(C_ptr(:) > 0);
    
    % Initialize PnonZero array
    PnonZero = zeros(1, NnonZeroCvalues);
    
    % Read non-zero probability values and compute Rnaive
    index = 1;
    for i = 1:Rtot
        if C_ptr(i) > 0
            PnonZero(index) = C_ptr(i); % / N;
            index = index + 1;
        end
    end
    
    % Rnaive is the number of non-null probability values
    Rnaive = NnonZeroCvalues;
    
    if Rnaive < Rtot
        % Initial value for Rexpected
        Rexpected = Rnaive;
        for i = 1:NnonZeroCvalues
            Rexpected = Rexpected - (1 - PnonZero(i))^N;
        end
        
        % Initial values for deltaRprevious and deltaR
        deltaRprevious = Rtot;
        deltaR = abs(Rnaive - Rexpected);
        
        R = Rnaive;
        while deltaR < deltaRprevious && R < Rtot
            R = R + 1;
            
            gamma = (R - Rnaive) * (1 - (N / (N + Rnaive))^(1 / N));
            
            Rexpected = 0;
            % Occupied bins
            for i = 1:NnonZeroCvalues
                Pbayes = (PnonZero(i) * N + 1) / (N + Rnaive) * (1 - gamma);
                Rexpected = Rexpected + 1 - (1 - Pbayes)^N;
            end
            
            % Non-occupied bins
            Pbayes = gamma / (R - Rnaive);
            Rexpected = Rexpected + (R - Rnaive) * (1 - (1 - Pbayes)^N);
            
            deltaRprevious = deltaR;
            deltaR = abs(Rnaive - Rexpected);
        end
        
        R = R - 1;
        if deltaR < deltaRprevious
            R = R + 1;
        end
    else
        R = Rtot;
    end
    
    % Estimating bias
    bias = (R - 1) / (2 * N * log(2));
end
