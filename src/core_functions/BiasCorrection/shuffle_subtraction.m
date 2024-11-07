function [corrected_v, naive_v, shuff_all, addOut] = shuffle_subtraction(inputs, outputs, corefunc, varargin)
% shuffle_subtraction - Compute bias-corrected and naive information-theoretic values
%
% Inputs:
%   - inputs: A cell array containing the input variables, where each cell corresponds to a different input source. 
%             The data can be structured as either:
%             - nDims X nTrials 
%             - nDims X nTimepoints X nTrials
%
%   - outputs: A cell array specifying the information measures to compute.
%
%   - corefunc: A function handle representing the core information-theoretic computation (e.g., @PID, @MI, @FIT, @II).
%
%   - varargin: Optional arguments as a structure (`opts`), including fields:
%              - shuff: The number of shuffling iterations for bias correction (default: 20).
%              - pid_constrained: A flag indicating whether to use a constrained PID (true/false) (default: true).
%              - chosen_atom: Specifies which PID atom ('Red', 'Syn', 'Unq1', 'Unq2') or other measure to compute (default: 'Syn').
%              - inputs_nD: Higher-dimensional input data for specific cases (if needed).
%
% Outputs:
%   - corrected_v: A cell array containing the bias-corrected values for each specified output, computed by subtracting
%                  the shuffled contributions from the naive values.
%   - naive_v: A cell array containing the naive (uncorrected) values, as computed from the original input-output 
%              relationships.
%   - shuff_all: A cell array containing the results of the shuffled computations across all trials and outputs, 
%                used for bias subtraction.
%   - addOut: Additional outputs specific to the core function used (e.g., decomposed PID atoms, FIT interaction measures).
%
% Note:
% Users can customize the core calculation function to extend the analysis for specific types of data 
% or information-theoretic measures. The function is designed to handle various dimensions of input data, 
% ensuring flexibility for different experimental setups.

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

opts = varargin{1};
if nargin > 5
    inputs_nD = varargin{2};
end
addOut = 0;
naive_opts = opts;
naive_opts.bias = 'naive';
nVars = length(inputs);
nTimepoints_all = 1;
Dims_tmp = size(inputs{1});
for var = 1:nVars
    Dims_var = size(inputs{var});
    if length(Dims_var) == 3
        nTrials = Dims_var(end);
        nTimepoints_all = [nTimepoints_all, Dims_var(2)];
    end
end
nTimepoints = max(nTimepoints_all);
corrected_v = cell(1, length(outputs));
if strcmp(func2str(corefunc), 'PID')
    naive_v = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
    shuff_all = repmat({zeros(opts.shuff, nTimepoints)},1, length(outputs));
    nSources = length(inputs)-1;
    if ~isfield(opts, 'pid_constrained')
        if nSources == 2
            opts.pid_constrained = true;
        else
            opts.pid_constrained = false;
        end
    end
    PID_naive = feval(corefunc, inputs, outputs, naive_opts);
    if opts.pid_constrained
        I1_naive  = cell2mat(MI({inputs{1}, inputs{end}}, {'I(A;B)'}, naive_opts));
        I2_naive  = cell2mat(MI({inputs{2}, inputs{end}}, {'I(A;B)'}, naive_opts));
        I12_naive = cell2mat(MI({cat(1, inputs{1}, inputs{2}), inputs{end}}, {'I(A;B)'}, naive_opts));
        I1_shuff_all = zeros(opts.shuff, nTimepoints);
        I2_shuff_all = zeros(opts.shuff, nTimepoints);
        I12_shuff_all = zeros(opts.shuff, nTimepoints);
    end
    PID_corrected = cell(1, length(outputs));
    for sIdx = 1:opts.shuff
        inputs_sh = inputs;
        for var = 1:length(inputs)
            if length(size(inputs_sh{var}))>2
                inputs_sh{var} = permute(inputs_sh{var}, [3, 1, 2]);
                shuffle_vals =  (shuffle_core(0, inputs_sh{var}, 0, [1, 0]));
                inputs_sh{var} = permute(shuffle_vals, [2, 3, 1]);
            else
                inputs_sh{var} = (shuffle_core(0, inputs_sh{var}', 0, [1, 0]))';
            end
        end
        shuff_v  = feval(corefunc, inputs_sh, outputs, naive_opts);
        for outIdx = 1:length(outputs)
            shuff_all{outIdx}(sIdx,:) = shuff_v{outIdx};
        end
        if opts.pid_constrained
            I1_shuff_all(sIdx,:) = cell2mat(MI({inputs_sh{1}, inputs_sh{end}}, {'I(A;B)'}, naive_opts));
            I2_shuff_all(sIdx,:) = cell2mat(MI({inputs_sh{2}, inputs_sh{end}}, {'I(A;B)'}, naive_opts));
            I12_shuff_all(sIdx,:) = cell2mat(MI({cat(1, inputs_sh{1}, inputs_sh{2}), inputs_sh{end}}, {'I(A;B)'}, naive_opts));
        end
    end

    if opts.pid_constrained
        I1_corrected = I1_naive-mean(I1_shuff_all,1);
        I2_corrected = I2_naive-mean(I2_shuff_all,1);
        I12_corrected = I12_naive-mean(I12_shuff_all,1);
    end
    for outIdx = 1:length(outputs)
        PID_corrected{outIdx} = PID_naive{outIdx} - mean(shuff_all{outIdx},1);
        naive_v{outIdx} = PID_naive{outIdx};
    end
    if opts.pid_constrained
        if ~isfield(opts, 'chosen_atom')
            opts.chosen_atom = 'Syn';
        end
        pos = find(strcmp(opts.chosen_atom, outputs));
        switch opts.chosen_atom
            case 'Red'
                red = PID_corrected{pos};
                syn = I12_corrected-I1_corrected-I2_corrected+red;
                unq1 = I1_corrected-red;
                unq2 = I2_corrected-red;
                %naive
                red_naive = PID_naive{pos};
                syn_naive  = I12_naive-I1_naive-I2_naive+red_naive;
                unq1_naive  = I1_naive-red_naive;
                unq2_naive  = I2_naive-red_naive;
            case 'Unq1'
                unq1 = PID_corrected{pos};
                red = I1_corrected-unq1;
                syn =  I12_corrected-I2_corrected-unq1;
                unq2 = I2_corrected-I1_corrected+unq1;
                %naive
                unq1_naive = PID_naive{pos};
                red_naive = I1_naive-unq1_naive;
                syn_naive =  I12_naive-I2_naive-unq1_naive;
                unq2_naive = I2_naive-I1_naive+unq1_naive;
            case 'Unq2'
                unq2 = PID_corrected{pos};
                red = I2_corrected-unq2;
                syn =  I12_corrected-I1_corrected-unq2;
                unq1 = I1_corrected-I2_corrected+unq2;
                %naive
                unq2_naive = PID_naive{pos};
                red_naive = I2_naive-unq2_naive;
                syn_naive =  I12_naive-I1_naive-unq2_naive;
                unq1_naive = I1_naive-I2_naive+unq2_naive;
            case 'Syn'
                syn = PID_corrected{pos};
                red = I1_corrected+I2_corrected-I12_corrected+syn;
                unq1 = I12_corrected-I2_corrected-syn;
                unq2 = I12_corrected-I1_corrected-syn;
                %naive
                syn_naive  = PID_naive{pos};
                red_naive  = I1_naive+I2_naive-I12_naive+syn_naive;
                unq1_naive  = I12_naive-I2_naive-syn_naive;
                unq2_naive  = I12_naive-I1_naive-syn_naive;
        end
        for i = 1:length(outputs)
            switch outputs{i}
                case 'Syn'
                    PID_corrected{i} = syn;
                    naive_v{i}       = syn_naive;
                case 'Red'
                    PID_corrected{i} = red;
                    naive_v{i}       = red_naive;
                case 'Unq1'
                    PID_corrected{i} = unq1;
                    naive_v{i}       = unq1_naive;
                case 'Unq2'
                    PID_corrected{i} = unq2;
                    naive_v{i}       = unq2_naive;
                case 'Unq'
                    PID_corrected{i} = unq1 + unq2;
                    naive_v{i}       = unq1_naive + unq2_naive;
                case 'Joint'
                    PID_corrected{i} = red+syn+unq1+unq2;
                    naive_v{i}       = red_naive+syn_naive+unq1_naive+unq2_naive;
                case 'Union'
                    PID_corrected{i} = (red+syn+unq1+unq2) - syn;
                    naive_v{i}       = (red_naive+syn_naive+unq1_naive+unq2_naive) - syn_naive;
                otherwise
                    PID_corrected{i} = NaN;
                    naive_v{i}       = NaN;
            end
        end
    end
    corrected_v = PID_corrected;
elseif strcmp(func2str(corefunc), 'FIT') || strcmp(func2str(corefunc), 'cFIT')
    atom1_corrected = cell(1, length(outputs));
    atom2_corrected = cell(1, length(outputs));
    shuff_all = cell(1, length(outputs));
    [naive_v,~,~,atom1_naive, atom2_naive] = feval(corefunc, inputs, outputs, naive_opts);
    addOut = cell(2,length(outputs));
    for sIdx = 1:opts.shuff
        inputs_sh = inputs;
        for var = 1:length(inputs)
            if length(size(inputs_sh{var}))>2
                inputs_sh{var} = permute(inputs_sh{var}, [3, 1, 2]);
                shuffle_vals =  (shuffle_core(0, inputs_sh{var}, 0, [1, 0]));
                inputs_sh{var} = permute(shuffle_vals, [2, 3, 1]);
            else
                inputs_sh{var} = (shuffle_core(0, inputs_sh{var}', 0, [1, 0]))';
            end
        end
        [FIT_shuff, ~, ~, atom1_shuff, atom2_shuff]  = feval(corefunc, inputs_sh, outputs, naive_opts);
        for outIdx = 1:length(outputs)
            shuff_all{outIdx} = FIT_shuff;
            addOut{1,outIdx} = [addOut{1,outIdx}, atom1_shuff{outIdx}];
            addOut{2,outIdx} = [addOut{2,outIdx}, atom2_shuff{outIdx}];
        end
    end
    for outIdx = 1:length(outputs)
        atom1_corrected{outIdx} = atom1_naive{outIdx} - mean(addOut{1,outIdx});
        atom2_corrected{outIdx} = atom2_naive{outIdx} - mean(addOut{1,outIdx});
        addOut{1,outIdx} = [addOut{1,outIdx},atom1_corrected{outIdx}];
        addOut{2,outIdx} = [addOut{2,outIdx},atom2_corrected{outIdx}];
        if  strcmp(func2str(corefunc), 'cFIT')
            atom1 = min(atom1_corrected{outIdx});
            atom2 = min(atom2_corrected{outIdx});
            corrected_v{outIdx} = atom1-atom2;
        else
            corrected_v{outIdx} = min(atom1_corrected{outIdx},atom2_corrected{outIdx});
        end
    end
elseif strcmp(func2str(corefunc), 'II')
    atom1_corrected = cell(1, length(outputs));
    atom2_corrected = cell(1, length(outputs));
    shuff_all = 0;
    atom1_shuffall = repmat({zeros(opts.shuff, nTimepoints)},1, length(outputs));
    atom2_shuffall = repmat({zeros(opts.shuff, nTimepoints)},2, length(outputs));
    [naive_v,~,~,atom1_naive, atom2_naive] = feval(corefunc, inputs, outputs, naive_opts);
    addOut = cell(2,length(outputs));
     for sIdx = 1:opts.shuff
        inputs_sh = inputs;
        for var = 1:length(inputs)
            if length(size(inputs_sh{var}))>2
                inputs_sh{var} = permute(inputs_sh{var}, [3, 1, 2]);
                shuffle_vals =  (shuffle_core(0, inputs_sh{var}, 0, [1, 0]));
                inputs_sh{var} = permute(shuffle_vals, [2, 3, 1]);
            else
                inputs_sh{var} = (shuffle_core(0, inputs_sh{var}', 0, [1, 0]))';
            end
        end
        [~, ~, ~, atom1_shuff, atom2_shuff]  = feval(corefunc, inputs_sh, outputs, naive_opts);
        for outIdx = 1:length(outputs)            
            atom1_shuffall{outIdx}(sIdx,:) = atom1_shuff{outIdx};
            atom2_shuffall{outIdx}(sIdx,:) = atom2_shuff{outIdx};
        end

     end
     for outIdx = 1:length(outputs)
         atom1_corrected{outIdx} = atom1_naive{outIdx} - mean(atom1_shuffall{outIdx},1);
         atom2_corrected{outIdx} = atom2_naive{outIdx} - mean(atom2_shuffall{outIdx},1);
         addOut{1,outIdx} = atom1_shuffall{outIdx};
         addOut{1,outIdx} = atom2_shuffall{outIdx};
         corrected_v{outIdx} = min(atom1_corrected{outIdx},atom2_corrected{outIdx});
     end
elseif strcmp(func2str(corefunc), 'MI') || strcmp(func2str(corefunc), 'TE') || strcmp(func2str(corefunc), 'cTE') || strcmp(func2str(corefunc), 'cMI')
    shuff_all = repmat({zeros(opts.shuff, nTimepoints)},1, length(outputs));
    naive_v = feval(corefunc, inputs, outputs, naive_opts);
    for sIdx = 1:opts.shuff
        inputs_sh = inputs;
        for var = 1:length(inputs)
            if length(size(inputs_sh{var}))>2
                inputs_sh{var} = permute(inputs_sh{var}, [3, 1, 2]);
                shuffle_vals =  (shuffle_core(0, inputs_sh{var}, 0, [1, 0]));
                inputs_sh{var} = permute(shuffle_vals, [2, 3, 1]);
            else
                inputs_sh{var} = (shuffle_core(0, inputs_sh{var}', 0, [1, 0]))';
            end
        end
        shuff_v  = feval(corefunc, inputs_sh, outputs, naive_opts);
        for outIdx = 1:length(outputs)
            shuff_all{outIdx}(sIdx,:) = shuff_v{outIdx};
        end
    end
     for outIdx = 1:length(outputs)      
        corrected_v{outIdx} = naive_v{outIdx} - mean(shuff_all{outIdx},1);
     end 
end

