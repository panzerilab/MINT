function [corrected_v, plugin_v, shuff_all, addOut] = shuffle_subtraction(inputs, outputs, corefunc, varargin)
% shuffle_subtraction - Compute bias-corrected and plugin information-theoretic values
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
%                  the shuffled contributions from the plugin values.
%   - plugin_v: A cell array containing the plugin (uncorrected) values, as computed from the original input-output 
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
plugin_opts = opts;
plugin_opts.bias = 'plugin';
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
    plugin_v = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
    shuff_all = repmat({zeros(opts.shuff, nTimepoints)},1, length(outputs));
    nSources = length(inputs)-1;
    if ~isfield(opts, 'pid_constrained')
        if nSources == 2
            opts.pid_constrained = true;
        else
            opts.pid_constrained = false;
        end
    end
    PID_plugin = feval(corefunc, inputs, outputs, plugin_opts);
    if opts.pid_constrained
        I1_plugin  = cell2mat(MI({inputs{1}, inputs{end}}, {'I(A;B)'}, plugin_opts));
        I2_plugin  = cell2mat(MI({inputs{2}, inputs{end}}, {'I(A;B)'}, plugin_opts));
        I12_plugin = cell2mat(MI({cat(1, inputs{1}, inputs{2}), inputs{end}}, {'I(A;B)'}, plugin_opts));
        I1_shuff_all = zeros(opts.shuff, nTimepoints);
        I2_shuff_all = zeros(opts.shuff, nTimepoints);
        I12_shuff_all = zeros(opts.shuff, nTimepoints);
    end
    PID_corrected = cell(1, length(outputs));
    if opts.parallel
        shuff_all_par = zeros(length(outputs),opts.shuff, nTimepoints);
        noutputs = length(outputs);
        parfor sIdx = 1:opts.shuff
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
            shuff_v  = feval(corefunc, inputs_sh, outputs, plugin_opts);
            for outIdx = 1:noutputs
                shuff_all_par(outIdx,sIdx,:) = shuff_v{outIdx};
            end
            if opts.pid_constrained
                I1_shuff_all(sIdx,:) = cell2mat(MI({inputs_sh{1}, inputs_sh{end}}, {'I(A;B)'}, plugin_opts));
                I2_shuff_all(sIdx,:) = cell2mat(MI({inputs_sh{2}, inputs_sh{end}}, {'I(A;B)'}, plugin_opts));
                I12_shuff_all(sIdx,:) = cell2mat(MI({cat(1, inputs_sh{1}, inputs_sh{2}), inputs_sh{end}}, {'I(A;B)'}, plugin_opts));
            end
        end
        for spar =1:noutputs
            shuff_all{spar}(:,:) = reshape(shuff_all_par(spar,:,:), [opts.shuff, nTimepoints]);
        end

    else
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
            shuff_v  = feval(corefunc, inputs_sh, outputs, plugin_opts);
            for outIdx = 1:length(outputs)
                shuff_all{outIdx}(sIdx,:) = shuff_v{outIdx};
            end
            if opts.pid_constrained
                I1_shuff_all(sIdx,:) = cell2mat(MI({inputs_sh{1}, inputs_sh{end}}, {'I(A;B)'}, plugin_opts));
                I2_shuff_all(sIdx,:) = cell2mat(MI({inputs_sh{2}, inputs_sh{end}}, {'I(A;B)'}, plugin_opts));
                I12_shuff_all(sIdx,:) = cell2mat(MI({cat(1, inputs_sh{1}, inputs_sh{2}), inputs_sh{end}}, {'I(A;B)'}, plugin_opts));
            end
        end
    end

    if opts.pid_constrained
        I1_corrected = I1_plugin-mean(I1_shuff_all,1);
        I2_corrected = I2_plugin-mean(I2_shuff_all,1);
        I12_corrected = I12_plugin-mean(I12_shuff_all,1);
    end
    for outIdx = 1:length(outputs)
        PID_corrected{outIdx} = PID_plugin{outIdx} - mean(shuff_all{outIdx},1);
        plugin_v{outIdx} = PID_plugin{outIdx};
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
                %plugin
                red_plugin = PID_plugin{pos};
                syn_plugin  = I12_plugin-I1_plugin-I2_plugin+red_plugin;
                unq1_plugin  = I1_plugin-red_plugin;
                unq2_plugin  = I2_plugin-red_plugin;
            case 'Unq1'
                unq1 = PID_corrected{pos};
                red = I1_corrected-unq1;
                syn =  I12_corrected-I2_corrected-unq1;
                unq2 = I2_corrected-I1_corrected+unq1;
                %plugin
                unq1_plugin = PID_plugin{pos};
                red_plugin = I1_plugin-unq1_plugin;
                syn_plugin =  I12_plugin-I2_plugin-unq1_plugin;
                unq2_plugin = I2_plugin-I1_plugin+unq1_plugin;
            case 'Unq2'
                unq2 = PID_corrected{pos};
                red = I2_corrected-unq2;
                syn =  I12_corrected-I1_corrected-unq2;
                unq1 = I1_corrected-I2_corrected+unq2;
                %plugin
                unq2_plugin = PID_plugin{pos};
                red_plugin = I2_plugin-unq2_plugin;
                syn_plugin =  I12_plugin-I1_plugin-unq2_plugin;
                unq1_plugin = I1_plugin-I2_plugin+unq2_plugin;
            case 'Syn'
                syn = PID_corrected{pos};
                red = I1_corrected+I2_corrected-I12_corrected+syn;
                unq1 = I12_corrected-I2_corrected-syn;
                unq2 = I12_corrected-I1_corrected-syn;
                %plugin
                syn_plugin  = PID_plugin{pos};
                red_plugin  = I1_plugin+I2_plugin-I12_plugin+syn_plugin;
                unq1_plugin  = I12_plugin-I2_plugin-syn_plugin;
                unq2_plugin  = I12_plugin-I1_plugin-syn_plugin;
        end
        for i = 1:length(outputs)
            switch outputs{i}
                case 'Syn'
                    PID_corrected{i} = syn;
                    plugin_v{i}       = syn_plugin;
                case 'Red'
                    PID_corrected{i} = red;
                    plugin_v{i}       = red_plugin;
                case 'Unq1'
                    PID_corrected{i} = unq1;
                    plugin_v{i}       = unq1_plugin;
                case 'Unq2'
                    PID_corrected{i} = unq2;
                    plugin_v{i}       = unq2_plugin;
                case 'Unq'
                    PID_corrected{i} = unq1 + unq2;
                    plugin_v{i}       = unq1_plugin + unq2_plugin;
                case 'Joint'
                    PID_corrected{i} = red+syn+unq1+unq2;
                    plugin_v{i}       = red_plugin+syn_plugin+unq1_plugin+unq2_plugin;
                case 'Union'
                    PID_corrected{i} = (red+syn+unq1+unq2) - syn;
                    plugin_v{i}       = (red_plugin+syn_plugin+unq1_plugin+unq2_plugin) - syn_plugin;
                otherwise
                    PID_corrected{i} = NaN;
                    plugin_v{i}       = NaN;
            end
        end
    end
    corrected_v = PID_corrected;
elseif strcmp(func2str(corefunc), 'FIT') || strcmp(func2str(corefunc), 'cFIT')
    atom1_corrected = cell(1, length(outputs));
    atom2_corrected = cell(1, length(outputs));
    shuff_all = cell(1, length(outputs));
    [plugin_v,~,~,atom1_plugin, atom2_plugin] = feval(corefunc, inputs, outputs, plugin_opts);
    addOut = cell(2,length(outputs));
    if opts.parallel
        shuff_all_par = zeros(length(outputs),opts.shuff, nTimepoints);
        addOut_par    = zeros(2,length(outputs),opts.shuff, nTimepoints);
        noutputs = length(outputs);
        parfor sIdx = 1:opts.shuff
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
            [FIT_shuff, ~, ~, atom1_shuff, atom2_shuff]  = feval(corefunc, inputs_sh, outputs, plugin_opts);
            for outIdx = 1:noutputs
                shuff_all_par(outIdx,sIdx,:) = FIT_shuff;
                addOut(:,outIdx,sIdx,:) =  [atom1_shuff{outIdx}; atom2_shuff{outIdx}];
            end
        end   
        for spar =1:noutputs
            shuff_all{spar}(:,:) = reshape(shuff_all_par(spar,:,:),opts.shuff, nTimepoints);
            addOut{spar}(:,:)    = reshape(addOut_par(spar,:,:),opts.shuff, nTimepoints);
        end
    else
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
            [FIT_shuff, ~, ~, atom1_shuff, atom2_shuff]  = feval(corefunc, inputs_sh, outputs, plugin_opts);
            for outIdx = 1:length(outputs)
                shuff_all{outIdx} = FIT_shuff;
                addOut{1,outIdx} = [addOut{1,outIdx}, atom1_shuff{outIdx}];
                addOut{2,outIdx} = [addOut{2,outIdx}, atom2_shuff{outIdx}];
            end
        end
    end
    for outIdx = 1:length(outputs)
        atom1_corrected{outIdx} = atom1_plugin{outIdx} - mean(addOut{1,outIdx});
        atom2_corrected{outIdx} = atom2_plugin{outIdx} - mean(addOut{1,outIdx});
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
    [plugin_v,~,~,atom1_plugin, atom2_plugin] = feval(corefunc, inputs, outputs, plugin_opts);
    addOut = cell(2,length(outputs));
    if opts.parallel
        noutputs = length(outputs);
        atom1_shuffall_par = zeros(noutputs, opts.shuff, nTimepoints);
        atom2_shuffall_par = zeros(noutputs, opts.shuff, nTimepoints);
        parfor sIdx = 1:opts.shuff
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
            [~, ~, ~, atom1_shuff, atom2_shuff]  = feval(corefunc, inputs_sh, outputs, plugin_opts);
            for outIdx = 1:noutputs            
                atom1_shuffall_par(outIdx,sIdx,:) = atom1_shuff{outIdx};
                atom2_shuffall_par(outIdx,sIdx,:) = atom2_shuff{outIdx};
            end
        end
        for spar=1:noutputs
            atom1_shuffall{noutputs} = reshape(atom1_shuffall_par(noutputs,:,:),opts.shuff, nTimepoints);
            atom2_shuffall{noutputs} = reshape(atom2_shuffall_par(noutputs,:,:),opts.shuff, nTimepoints);
        end

    else
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
            [~, ~, ~, atom1_shuff, atom2_shuff]  = feval(corefunc, inputs_sh, outputs, plugin_opts);
            for outIdx = 1:length(outputs)            
                atom1_shuffall{outIdx}(sIdx,:) = atom1_shuff{outIdx};
                atom2_shuffall{outIdx}(sIdx,:) = atom2_shuff{outIdx};
            end
        end
    end
    for outIdx = 1:length(outputs)
         atom1_corrected{outIdx} = atom1_plugin{outIdx} - mean(atom1_shuffall{outIdx},1);
         atom2_corrected{outIdx} = atom2_plugin{outIdx} - mean(atom2_shuffall{outIdx},1);
         addOut{1,outIdx} = atom1_shuffall{outIdx};
         addOut{1,outIdx} = atom2_shuffall{outIdx};
         corrected_v{outIdx} = min(atom1_corrected{outIdx},atom2_corrected{outIdx});
     end
elseif strcmp(func2str(corefunc), 'MI') || strcmp(func2str(corefunc), 'TE') || strcmp(func2str(corefunc), 'cTE') || strcmp(func2str(corefunc), 'cMI')
    shuff_all = repmat({zeros(opts.shuff, nTimepoints)},1, length(outputs));
    plugin_v = feval(corefunc, inputs, outputs, plugin_opts);
    if opts.parallel
        noutputs = length(outputs);
        shuffall_par = zeros(noutputs, opts.shuff, nTimepoints);
        parfor sIdx = 1:opts.shuff
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
            shuff_v  = feval(corefunc, inputs_sh, outputs, plugin_opts);
            for outIdx = 1:noutputs
                shuffall_par(outIdx,sIdx,:) = shuff_v{outIdx};
            end
        end
        for spar=1:noutputs
            shuff_all{noutputs} = reshape(shuffall_par(noutputs,:,:),opts.shuff, nTimepoints);
        end
    else
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
            shuff_v  = feval(corefunc, inputs_sh, outputs, plugin_opts);
            for outIdx = 1:length(outputs)
                shuff_all{outIdx}(sIdx,:) = shuff_v{outIdx};
            end
        end
    end
    for outIdx = 1:length(outputs)      
        corrected_v{outIdx} = plugin_v{outIdx} - mean(shuff_all{outIdx},1);
    end 
end

