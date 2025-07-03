function [corrected_v, plugin_v] = extrapolation(inputs, outputs, corr, corefunc, varargin)
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
if ~isfield(opts, 'timeseries')
    opts.timeseries = false;
end
if ~isfield(opts, 'shuff')
    opts.shuff = 0;
end
if opts.timeseries
    nTimepoints = size(inputs{1},2);
else
    nTimepoints = 1;
end

corrected_v = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
xtrp = opts.xtrp;

if strcmp(func2str(corefunc), 'PID')
    plugin_opts = opts;
    plugin_opts.bias = 'plugin';
    nSources = length(inputs)-1;
    if ~isfield(opts, 'pid_constrained')
        if nSources == 2
            opts.pid_constrained = true;
        else
            opts.pid_constrained = false;
        end
    end
    plugin_v = feval(corefunc, inputs, outputs, plugin_opts);
    corrected_v = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
    if opts.pid_constrained
        I1_plugin  = cell2mat(MI({inputs{1}, inputs{end}}, {'I(A;B)'}, plugin_opts));
        I2_plugin  = cell2mat(MI({inputs{2}, inputs{end}}, {'I(A;B)'}, plugin_opts));
        I12_plugin = cell2mat(MI({cat(1, inputs{1}, inputs{2}), inputs{end}}, {'I(A;B)'}, plugin_opts));
        I1_corrected = zeros(1, nTimepoints);
        I2_corrected = zeros(1, nTimepoints);
        I12_corrected = zeros(1, nTimepoints);
    end
    nTrials = size(inputs{1},length(size(inputs{1})));

    if opts.parallel
        noutputs = length(outputs);
        corrected_v_par = zeros(opts.xtrp,nTimepoints);
        I1_corrected_par = zeros(opts.xtrp,nTimepoints);
        I2_corrected_par = zeros(opts.xtrp,nTimepoints);
        I12_corrected_par = zeros(opts.xtrp,nTimepoints);
        parfor i=1:opts.xtrp
            randidx = randperm(nTrials, nTrials);
            npartition = [1 2 4];
            PID2 =  repmat({zeros(1,nTimepoints)}, 1, length(outputs));
            PID4 =  repmat({zeros(1,nTimepoints)}, 1, length(outputs));
            I1_2  = 0;
            I2_2  = 0;
            I12_2 = 0;
            I1_4  = 0;
            I2_4  = 0;
            I12_4 = 0;
            inputs_s = inputs;
            idx = repmat({':'}, 1, ndims(inputs{1}));
            idx{end} = randidx;
            for var = 1:length(inputs)
                inputs_s{var} = inputs{var}(idx{:});
            end
            for np=2:length(npartition)
                for pidx = 1:npartition(np)
                    inputs_p = partition(inputs_s, npartition(np), pidx,1);
                    PID_p = feval(corefunc, inputs_p, outputs, plugin_opts);
                    for outIdx = 1:noutputs
                        if npartition(np)==2
                            PID2{outIdx} = PID2{outIdx} + PID_p{outIdx}/2;
                        elseif npartition(np)==4
                            PID4{outIdx} =  PID4{outIdx} + PID_p{outIdx}/4;
                        end
                    end
                    if opts.pid_constrained
                        I1_p = MI({inputs_p{1}, inputs_p{end}}, {'I(A;B)'}, plugin_opts);
                        I2_p = MI({inputs_p{2}, inputs_p{end}}, {'I(A;B)'}, plugin_opts);
                        I12_p = MI({cat(1, inputs_p{1}, inputs_p{2}), inputs_p{end}}, {'I(A;B)'}, plugin_opts);
                        if npartition(np)==2
                            I1_2 = I1_2 + I1_p{1}/2;
                            I2_2 = I2_2 + I2_p{1}/2;
                            I12_2 = I12_2 + I12_p{1}/2;
                        elseif npartition(np)==4
                            I1_4 = I1_4 + I1_p{1}/4;
                            I2_4 = I2_4 + I2_p{1}/4;
                            I12_4 = I12_4 + I12_p{1}/4;
                        end
                    end
                end
            end
            x_extrap = npartition./nTrials;
            if strcmp(opts.bias,'qe')
                for t = 1:nTimepoints
                    for outIdx = 1:noutputs
                        y = [plugin_v{outIdx}(t), PID2{outIdx}(t), PID4{outIdx}(t)];
                        p = polyfit(x_extrap, y, 2);
                        corrected_v_par(outIdx, i, t) = p(3);
                    end
                end
            elseif strcmp(opts.bias,'le')
                for t = 1:nTimepoints
                    for outIdx = 1:noutputs
                        y = [plugin_v{outIdx}(t), PID2{outIdx}(t)];
                        p = polyfit(x_extrap, y, 1);
                        corrected_v_par(outIdx, i, t) = p(2);
                    end
                end
            end
            if opts.pid_constrained
                if strcmp(opts.bias,'qe')
                    for t = 1:nTimepoints
                        y = [I1_plugin(t), I1_2(t), I1_4(t)];
                        p = polyfit(x_extrap, y, 2);
                        I1_corrected_par(i,t) = p(3);
                        y = [I2_plugin(t), I2_2(t), I2_4(t)];
                        p = polyfit(x_extrap, y, 2);
                        I2_corrected_par(i,t) =  p(3);
                        y = [I12_plugin(t), I12_2(t), I12_4(t)];
                        p = polyfit(x_extrap, y, 2);
                        I12_corrected_par(i,t) = p(3);
                    end
                elseif strcmp(opts.bias,'le')
                    for t = 1:nTimepoints
                        y = [I1_plugin(t), I1_2(t)];
                        p = polyfit(x_extrap, y, 1);
                        I1_corrected_par(i,t) = p(2);
                        y = [I2_plugin(t), I2_2(t)];
                        p = polyfit(x_extrap, y, 1);
                        I2_corrected_par(i,t) = p(2);
                        y = [I12_plugin(t), I12_2(t)];
                        p = polyfit(x_extrap, y, 1);
                        I12_corrected_par(i,t) = p(2);
                    end
                end
            end
        end
        for spar=1:noutputs
            corrected_v{spar}(1,:) = mean(corrected_v_par(spar,:,:),2);
        end
        I1_corrected = mean(I1_corrected_par,1);
        I2_corrected = mean(I2_corrected_par,1);
        I12_corrected = mean(I12_corrected_par,1);

    else
        for i=1:opts.xtrp
            randidx = randperm(nTrials, nTrials);
            npartition = [1 2 4];
            PID2 =  repmat({zeros(1,nTimepoints)}, 1, length(outputs));
            PID4 =  repmat({zeros(1,nTimepoints)}, 1, length(outputs));
            I1_2  = 0;
            I2_2  = 0;
            I12_2 = 0;
            I1_4  = 0;
            I2_4  = 0;
            I12_4 = 0;
            inputs_s = inputs;
            idx = repmat({':'}, 1, ndims(inputs{1}));
            idx{end} = randidx;
            for var = 1:length(inputs)
                inputs_s{var} = inputs{var}(idx{:});
            end
            for np=2:length(npartition)
                for pidx = 1:npartition(np)
                    inputs_p = partition(inputs_s, npartition(np), pidx,1);
                    PID_p = feval(corefunc, inputs_p, outputs, plugin_opts);
                    for outIdx = 1:length(outputs)
                        if npartition(np)==2
                            PID2{outIdx} = PID2{outIdx} + PID_p{outIdx}/2;
                        elseif npartition(np)==4
                            PID4{outIdx} =  PID4{outIdx} + PID_p{outIdx}/4;
                        end
                    end
                    if opts.pid_constrained
                        I1_p = MI({inputs_p{1}, inputs_p{end}}, {'I(A;B)'}, plugin_opts);
                        I2_p = MI({inputs_p{2}, inputs_p{end}}, {'I(A;B)'}, plugin_opts);
                        I12_p = MI({cat(1, inputs_p{1}, inputs_p{2}), inputs_p{end}}, {'I(A;B)'}, plugin_opts);
                        if npartition(np)==2
                            I1_2 = I1_2 + I1_p{1}/2;
                            I2_2 = I2_2 + I2_p{1}/2;
                            I12_2 = I12_2 + I12_p{1}/2;
                        elseif npartition(np)==4
                            I1_4 = I1_4 + I1_p{1}/4;
                            I2_4 = I2_4 + I2_p{1}/4;
                            I12_4 = I12_4 + I12_p{1}/4;
                        end
                    end
                end
            end
            x_extrap = npartition./nTrials;
            if strcmp(opts.bias,'qe')
                for t = 1:nTimepoints
                    for outIdx = 1:length(outputs)
                        y = [plugin_v{outIdx}(t), PID2{outIdx}(t), PID4{outIdx}(t)];
                        p = polyfit(x_extrap, y, 2);
                        corrected_v{outIdx}(1,t)  =  corrected_v{outIdx}(t) + p(3)/xtrp;
                    end
                end
            elseif strcmp(opts.bias,'le')
                for t = 1:nTimepoints
                    for outIdx = 1:length(outputs)
                        y = [plugin_v{outIdx}(t), PID2{outIdx}(t)];
                        p = polyfit(x_extrap, y, 1);
                        corrected_v{outIdx}(1,t) = corrected_v{outIdx}(t) + p(2)/xtrp;
                    end
                end
            end
            if opts.pid_constrained
                if strcmp(opts.bias,'qe')
                    for t = 1:nTimepoints
                        y = [I1_plugin(t), I1_2(t), I1_4(t)];
                        p = polyfit(x_extrap, y, 2);
                        I1_corrected(t) = I1_corrected(t) + p(3)/xtrp;
                        y = [I2_plugin(t), I2_2(t), I2_4(t)];
                        p = polyfit(x_extrap, y, 2);
                        I2_corrected(t) = I2_corrected(t) + p(3)/xtrp;
                        y = [I12_plugin(t), I12_2(t), I12_4(t)];
                        p = polyfit(x_extrap, y, 2);
                        I12_corrected(t) = I12_corrected(t) + p(3)/xtrp;
                    end
                elseif strcmp(opts.bias,'le')
                    for t = 1:nTimepoints
                        y = [I1_plugin(t), I1_2(t)];
                        p = polyfit(x_extrap, y, 1);
                        I1_corrected(t) = I1_corrected(t) +  p(2)/xtrp;
                        y = [I2_plugin(t), I2_2(t)];
                        p = polyfit(x_extrap, y, 1);
                        I2_corrected(t) = I2_corrected(t) +  p(2)/xtrp;
                        y = [I12_plugin(t), I12_2(t)];
                        p = polyfit(x_extrap, y, 1);
                        I12_corrected(t) = I12_corrected(t) +  p(2)/xtrp;
                    end
                end
            end
        end

    end
    if opts.pid_constrained
        if ~isfield(opts, 'chosen_atom')
            opts.chosen_atom = 'Syn';
        end
        pos = find(strcmp(opts.chosen_atom, outputs));
        switch opts.chosen_atom
            case'Red'
                red = corrected_v{pos};
                syn = I12_corrected-I1_corrected-I2_corrected+red;
                unq1 = I1_corrected-red;
                unq2 = I2_corrected-red;
                %plugin
                red_plugin = plugin_v{pos};
                syn_plugin  = I12_plugin-I1_plugin-I2_plugin+red_plugin;
                unq1_plugin  = I1_plugin-red_plugin;
                unq2_plugin  = I2_plugin-red_plugin;
            case 'Unq1'
                unq1 = corrected_v{pos};
                red = I1_corrected-unq1;
                syn =  I12_corrected-I2_corrected-unq1;
                unq2 = I2_corrected-I1_corrected+unq1;
                %plugin
                unq1_plugin = plugin_v{pos};
                red_plugin = I1_plugin-unq1_plugin;
                syn_plugin =  I12_plugin-I2_plugin-unq1_plugin;
                unq2_plugin = I2_plugin-I1_plugin+unq1_plugin;
            case 'Unq2'
                unq2 = corrected_v{pos};
                red = I2_corrected-unq2;
                syn =  I12_corrected-I1_corrected-unq2;
                unq1 = I1_corrected-I2_corrected+unq2;
                %plugin
                unq2_plugin = plugin_v{pos};
                red_plugin = I2_plugin-unq2_plugin;
                syn_plugin =  I12_plugin-I1_plugin-unq2_plugin;
                unq1_plugin = I1_plugin-I2_plugin+unq2_plugin;
            case 'Syn'
                syn = corrected_v{pos};
                red = I1_corrected+I2_corrected-I12_corrected+syn;
                unq1 = I12_corrected-I2_corrected-syn;
                unq2 = I12_corrected-I1_corrected-syn;
                %plugin
                syn_plugin  = plugin_v{pos};
                red_plugin  = I1_plugin+I2_plugin-I12_plugin+syn_plugin;
                unq1_plugin  = I12_plugin-I2_plugin-syn_plugin;
                unq2_plugin  = I12_plugin-I1_plugin-syn_plugin;
        end
        for i = 1:length(outputs)
            switch outputs{i}
                case 'Syn'
                    corrected_v{i} = syn;
                    plugin_v{i}       = syn_plugin;
                case 'Red'
                    corrected_v{i} = red;
                    plugin_v{i}       = red_plugin;
                case 'Unq1'
                    corrected_v{i} = unq1;
                    plugin_v{i}       = unq1_plugin;
                case 'Unq2'
                    corrected_v{i} = unq2;
                    plugin_v{i}       = unq2_plugin;
                case 'Unq'
                    corrected_v{i} = unq1 + unq2;
                    plugin_v{i}       = unq1_plugin + unq2_plugin;
                case 'Joint'
                    corrected_v{i} = red+syn+unq1+unq2;
                    plugin_v{i}       = red_plugin+syn_plugin+unq1_plugin+unq2_plugin;
                case 'Union'
                    corrected_v{i} = (red+syn+unq1+unq2) - syn;
                    plugin_v{i}       = (red_plugin+syn_plugin+unq1_plugin+unq2_plugin) - syn_plugin;
                otherwise
                    corrected_v{i} = NaN;
                    plugin_v{i}       = NaN;
            end
        end
    end
elseif strcmp(func2str(corefunc), 'FIT') || strcmp(func2str(corefunc), 'cFIT')
    plugin_opts = opts;
    plugin_opts.bias = 'plugin';
    [~,~,~,atom1_plugin, atom2_plugin] = feval(corefunc, inputs, outputs, plugin_opts);
    numAtoms = length(atom1_plugin);
    atom1_corr = repmat({zeros(1,numAtoms)}, 1, length(outputs));
    atom2_corr = repmat({zeros(1,numAtoms)}, 1, length(outputs));
    plugin_opts = opts;
    plugin_opts.bias = 'plugin';
    plugin_opts.recall = true;
    nTrials = size(inputs{1},length(size(inputs{1})));

    if opts.parallel
    else
    end
    for i=1:opts.xtrp
        randidx = randperm(nTrials, nTrials);
        npartition = [1 2 4];
        atom1_2 =  repmat({zeros(1,numAtoms)}, 1, length(outputs));
        atom1_4 =  repmat({zeros(1,numAtoms)}, 1, length(outputs));
        atom2_2 =  repmat({zeros(1,numAtoms)}, 1, length(outputs));
        atom2_4 =  repmat({zeros(1,numAtoms)}, 1, length(outputs));
        inputs_s = inputs;
        idx = repmat({':'}, 1, ndims(inputs{1}));
        idx{end} = randidx;
        for var = 1:length(inputs)
            inputs_s{var} = inputs{var}(idx{:});
        end
        for np=2:length(npartition)
            for pidx = 1:npartition(np)
                inputs_p = partition(inputs_s, npartition(np), pidx,1);
                [~,~,~,atom1_tmp, atom2_tmp] = feval(corefunc, inputs_p,outputs, plugin_opts);
                for idx = 1:numAtoms
                    for outIdx = 1:length(outputs)
                        if npartition(np)==2
                            atom1_2{outIdx}(idx) = atom1_2{outIdx}(idx) + atom1_tmp{outIdx}(idx)/2;
                            atom2_2{outIdx}(idx) = atom2_2{outIdx}(idx) + atom2_tmp{outIdx}(idx)/2;
                        elseif npartition(np)==4
                            atom1_4{outIdx}(idx)  = atom1_4{outIdx}(idx)  + atom1_tmp{outIdx}(idx)/4;
                            atom2_4{outIdx}(idx)  = atom2_4{outIdx}(idx)  + atom2_tmp{outIdx}(idx)/4;
                        end
                    end
                end
            end
        end
        x_extrap = npartition./nTrials;
        for outIdx = 1:length(outputs)
            for idx = 1:numAtoms
                if strcmp(opts.bias,'qe')
                    atom1_corr_tmp = polyfit(x_extrap, [atom1_plugin{outIdx}(idx), atom1_2{outIdx}(idx), atom1_4{outIdx}(idx)], 2);
                    atom2_corr_tmp = polyfit(x_extrap, [atom2_plugin{outIdx}(idx), atom2_2{outIdx}(idx), atom2_4{outIdx}(idx)], 2);
                    atom1_corr{outIdx}(idx) = atom1_corr{outIdx}(idx) + (atom1_corr_tmp(3)/xtrp);
                    atom2_corr{outIdx}(idx) = atom2_corr{outIdx}(idx) + (atom2_corr_tmp(3)/xtrp);
                elseif strcmp(opts.bias,'le')
                    atom1_corr_tmp = polyfit(x_extrap(1:2), [atom1_plugin{outIdx}(idx), atom1_2{outIdx}(idx)], 1);
                    atom2_corr_tmp = polyfit(x_extrap(1:2), [atom2_plugin{outIdx}(idx), atom2_2{outIdx}(idx)], 1);
                    atom1_corr{outIdx}(idx) = atom1_corr{outIdx}(idx) + (atom1_corr_tmp(2)/xtrp);
                    atom2_corr{outIdx}(idx) = atom2_corr{outIdx}(idx) + (atom2_corr_tmp(2)/xtrp);
                end
            end
        end
    end
    for outIdx = 1:length(outputs)
        if strcmp(func2str(corefunc), 'cFIT')
            atom1 = min(atom1_corr{outIdx});
            atom2 = min(atom2_corr{outIdx});
            corrected_v{outIdx} = atom1-atom2;
            atom1_plugin = min(atom1_plugin{outIdx});
            atom2_plugin = min(atom2_plugin{outIdx});
            plugin_v{outIdx} = atom1_plugin-atom2_plugin;
        else
            corrected_v{outIdx} = min(atom1_corr{outIdx},atom2_corr{outIdx});
            plugin_v{outIdx} = min(atom1_plugin{outIdx}, atom2_plugin{outIdx});
        end
    end
elseif strcmp(func2str(corefunc), 'II')
    atom1_corr = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
    atom2_corr = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
    plugin_opts = opts;
    plugin_opts.bias = 'plugin';
    plugin_opts.recall = true;
    [plugin_v,~,~,atom1_plugin, atom2_plugin] = feval(corefunc, inputs, outputs, plugin_opts);
    nTrials = size(inputs{1},length(size(inputs{1})));
    if opts.pid_constrained
        if ~isfield(opts, 'chosen_atom')
            opts.chosen_atom = 'Syn';
        end

        for i=1:opts.xtrp
            randidx = randperm(nTrials, nTrials);
            npartition = [1 2 4];
            PID2S =  repmat({zeros(1,nTimepoints)}, 1, 4);
            PID4S =  repmat({zeros(1,nTimepoints)}, 1, 4);
            PID2C =  repmat({zeros(1,nTimepoints)}, 1, 4);
            PID4C =  repmat({zeros(1,nTimepoints)}, 1, 4);
            I1S_2  = 0;
            I2S_2  = 0;
            I12S_2 = 0;
            I1S_4  = 0;
            I2S_4  = 0;
            I12S_4 = 0;

            I1C_2  = 0;
            I2C_2  = 0;
            I12C_2 = 0;
            I1C_4  = 0;
            I2C_4  = 0;
            I12C_4 = 0;
            inputs_s = inputs;
            idx = repmat({':'}, 1, ndims(inputs{1}));
            idx{end} = randidx;
            for var = 1:length(inputs)
                inputs_s{var} = inputs{var}(idx{:});
            end
            for np=2:length(npartition)
                for pidx = 1:npartition(np)
                    inputs_p = partition(inputs_s, npartition(np), pidx,1);
                    PID_pC = feval(corefunc, inputs_p, {'PID_atoms'}, plugin_opts);
                    PID_pS = feval(corefunc, {inputs_p{3}, inputs_p{2}, inputs_p{1}}, {'PID_atoms'}, plugin_opts);
                    for outIdx = 1:length(outputs_pid)
                        if npartition(np)==2
                            PID2S{outIdx} = PID2S{outIdx} + PID_pS{outIdx}/2;
                            PID2C{outIdx} = PID2C{outIdx} + PID_pC{outIdx}/2;
                        elseif npartition(np)==4
                            PID4S{outIdx} = PID4S{outIdx} + PID_pS{outIdx}/4;
                            PID4C{outIdx} = PID4C{outIdx} + PID_pC{outIdx}/4;
                        end
                    end
                    I1_pC = MI({inputs_p{1}, inputs_p{end}}, {'I(A;B)'}, plugin_opts);
                    I2_pC = MI({inputs_p{2}, inputs_p{end}}, {'I(A;B)'}, plugin_opts);
                    I12_pC = MI({cat(1, inputs_p{1}, inputs_p{2}), inputs_p{end}}, {'I(A;B)'}, plugin_opts);

                    I1_pS = I1_pC ; %MI({inputs_p{1}, inputs_p{end}}, {'I(A;B)'}, plugin_opts);
                    I2_pS = MI({inputs_p{2}, inputs_p{1}}, {'I(A;B)'}, plugin_opts);
                    I12_pS = MI({cat(1, inputs_p{end}, inputs_p{2}), inputs_p{1}}, {'I(A;B)'}, plugin_opts);
                    if npartition(np)==2
                        I1S_2 = I1S_2   + I1_pS{1}/2;
                        I2S_2 = I2S_2   + I2_pS{1}/2;
                        I12S_2 = I12S_2 + I12_pS{1}/2;

                        I1C_2 = I1C_2   + I1_pC{1}/2;
                        I2C_2 = I2C_2   + I2_pC{1}/2;
                        I12C_2 = I12C_2 + I12_pC{1}/2;
                    elseif npartition(np)==4
                        I1S_4  = I1S_4  + I1_pS{1} /4;
                        I2S_4  = I2S_4  + I2_pS{1} /4;
                        I12S_4 = I12S_4 + I12_pS{1}/4;

                        I1C_4  = I1C_4  + I1_pC{1} /4;
                        I2C_4  = I2C_4  + I2_pC{1} /4;
                        I12C_4 = I12C_4 + I12_pC{1}/4;

                    end
                end
            end
            x_extrap = npartition./nTrials;
            if strcmp(opts.bias,'qe')
                for t = 1:nTimepoints
                    for outIdx = 1:length(outputs)
                        y = [plugin_v{outIdx}(t), PID2{outIdx}(t), PID4{outIdx}(t)];
                        p = polyfit(x_extrap, y, 2);
                        corrected_v{outIdx}(1,t)  =  corrected_v{outIdx}(t) + p(3)/xtrp;
                    end
                end
            elseif strcmp(opts.bias,'le')
                for t = 1:nTimepoints
                    for outIdx = 1:length(outputs)
                        y = [plugin_v{outIdx}(t), PID2{outIdx}(t)];
                        p = polyfit(x_extrap, y, 1);
                        corrected_v{outIdx}(1,t) = corrected_v{outIdx}(t) + p(2)/xtrp;
                    end
                end
            end
            if opts.pid_constrained
                if strcmp(opts.bias,'qe')
                    for t = 1:nTimepoints
                        y = [I1_plugin(t), I1_2(t), I1_4(t)];
                        p = polyfit(x_extrap, y, 2);
                        I1_corrected(t) = I1_corrected(t) + p(3)/xtrp;
                        y = [I2_plugin(t), I2_2(t), I2_4(t)];
                        p = polyfit(x_extrap, y, 2);
                        I2_corrected(t) = I2_corrected(t) + p(3)/xtrp;
                        y = [I12_plugin(t), I12_2(t), I12_4(t)];
                        p = polyfit(x_extrap, y, 2);
                        I12_corrected(t) = I12_corrected(t) + p(3)/xtrp;
                    end
                elseif strcmp(opts.bias,'le')
                    for t = 1:nTimepoints
                        y = [I1_plugin(t), I1_2(t)];
                        p = polyfit(x_extrap, y, 1);
                        I1_corrected(t) = I1_corrected(t) +  p(2)/xtrp;
                        y = [I2_plugin(t), I2_2(t)];
                        p = polyfit(x_extrap, y, 1);
                        I2_corrected(t) = I2_corrected(t) +  p(2)/xtrp;
                        y = [I12_plugin(t), I12_2(t)];
                        p = polyfit(x_extrap, y, 1);
                        I12_corrected(t) = I12_corrected(t) +  p(2)/xtrp;
                    end
                end
            end
        end

    else

        for i=1:opts.xtrp
            randidx = randperm(nTrials, nTrials);
            npartition = [1 2 4];
            atom1_2 = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
            atom1_4 = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
            atom2_2 = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
            atom2_4 = repmat({zeros(1,nTimepoints)}, 1, length(outputs));

            inputs_s = inputs;
            idx = repmat({':'}, 1, ndims(inputs{1}));
            idx{end} = randidx;
            for var = 1:length(inputs)
                inputs_s{var} = inputs{var}(idx{:});
            end
            for np=2:length(npartition)
                for pidx = 1:npartition(np)
                    inputs_p = partition(inputs_s, npartition(np), pidx,1);
                    [~,~,~,atom1_tmp, atom2_tmp] = feval(corefunc, inputs_p,outputs, plugin_opts);
                    for outIdx = 1:length(outputs)
                        if npartition(np)==2
                            atom1_2{outIdx} = atom1_2{outIdx} + atom1_tmp{outIdx}/2;
                            atom2_2{outIdx} = atom2_2{outIdx} + atom2_tmp{outIdx}/2;
                        elseif npartition(np)==4
                            atom1_4{outIdx} = atom1_4{outIdx} + atom1_tmp{outIdx}/4;
                            atom2_4{outIdx} = atom2_4{outIdx} + atom2_tmp{outIdx}/4;
                        end
                    end
                end
            end
            x_extrap = npartition./nTrials;

            for outIdx = 1:length(outputs)
                if strcmp(opts.bias,'qe')
                    for t = 1:nTimepoints
                        y = [atom1_plugin{outIdx}(t), atom1_2{outIdx}(t), atom1_4{outIdx}(t)];
                        p = polyfit(x_extrap, y, 2);
                        atom1_corr{outIdx}(1,t)  =  atom1_corr{outIdx}(t) + p(3)/xtrp;

                        y = [atom2_plugin{outIdx}(t), atom2_2{outIdx}(t), atom2_4{outIdx}(t)];
                        p = polyfit(x_extrap, y, 2);
                        atom2_corr{outIdx}(1,t)  =  atom2_corr{outIdx}(t) + p(3)/xtrp;
                    end
                elseif strcmp(opts.bias,'le')
                    for t = 1:nTimepoints
                        y = [atom1_plugin{outIdx}(t), atom1_2{outIdx}(t)];
                        p = polyfit(x_extrap, y, 1);
                        atom1_corr{outIdx}(1,t)  =  atom1_corr{outIdx}(t) + p(2)/xtrp;

                        y = [atom2_plugin{outIdx}(t), atom2_2{outIdx}(t)];
                        p = polyfit(x_extrap, y, 1);
                        atom2_corr{outIdx}(1,t)  =  atom2_corr{outIdx}(t) + p(2)/xtrp;
                    end
                end
            end
        end
        for outIdx = 1:length(outputs)
            corrected_v{outIdx} = min(atom1_corr{outIdx},atom2_corr{outIdx});
        end
    end
elseif strcmp(func2str(corefunc), 'H') || strcmp(func2str(corefunc), 'MI') || strcmp(func2str(corefunc), 'TE') || strcmp(func2str(corefunc), 'cTE') || strcmp(func2str(corefunc), 'cMI')
    plugin_opts = opts;
    plugin_opts.bias = 'plugin';
    plugin_opts.computeNulldist = false;
    plugin_v = feval(corefunc, inputs, outputs, plugin_opts);
    nTrials = size(inputs{1},length(size(inputs{1})));
    for i=1:opts.xtrp
        randidx = randperm(nTrials, nTrials);
        npartition = [1 2 4];
        part2 = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
        part4 = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
        inputs_s = inputs;
        for var = 1:length(inputs)
            idx = repmat({':'}, 1, ndims(inputs{var}));
            idx{end} = randidx;
            inputs_s{var} = inputs{var}(idx{:});
        end
        for np=2:length(npartition)
            for pidx = 1:npartition(np)
                inputs_p = partition(inputs_s, npartition(np), pidx, 0);
                value_tmp = feval(corefunc, inputs_p, outputs, plugin_opts);
                for outIdx = 1:length(outputs)
                    for t = 1:nTimepoints
                        if npartition(np)==2
                            part2{outIdx}(1,t) = part2{outIdx}(1,t) + value_tmp{outIdx}(1,t)/2;
                        elseif npartition(np)==4
                            part4{outIdx}(1,t) = part4{outIdx}(1,t) + value_tmp{outIdx}(1,t)/4;
                        end
                    end
                end
            end
        end
        x_extrap = npartition./nTrials;
        for outIdx = 1:length(outputs)
            for t = 1:nTimepoints
                if strcmp(opts.bias,'qe') ||strcmp(opts.bias,'qe_shuffSub')
                    p = polyfit(x_extrap, [plugin_v{outIdx}(1,t), part2{outIdx}(1,t), part4{outIdx}(1,t)], 2);
                    corrected_v{outIdx}(1,t) = corrected_v{outIdx}(1,t) +(p(3)/xtrp);
                elseif strcmp(opts.bias,'le')||strcmp(opts.bias,'le_shuffSub')
                    p = polyfit(x_extrap(1:2), [plugin_v{outIdx}(1,t), part2{outIdx}(1,t)], 1);
                    corrected_v{outIdx}(1,t) = corrected_v{outIdx}(1,t) +(p(2)/xtrp);
                end
            end
        end
    end
end






