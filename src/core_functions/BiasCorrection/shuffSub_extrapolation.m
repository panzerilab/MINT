function [corrected_v, naive_v] = shuffSub_extrapolation(inputs, outputs, corr, corefunc, varargin)
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
    opts.shuff = 1;
end
if opts.timeseries
    nTimepoints = size(inputs{1},2);
else
    nTimepoints = 1;
end
xtrp = opts.xtrp;
nTrials = size(inputs{1},length(size(inputs{1})));

naive_opts = opts;
naive_opts.bias = 'naive';
naive_opts.computeNulldist = false;


if strcmp(func2str(corefunc), 'MI') || strcmp(func2str(corefunc), 'TE')  || strcmp(func2str(corefunc), 'cTE')  || strcmp(func2str(corefunc), 'cMI')
    corrected_v = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
    correctedSh = repmat({zeros(opts.shuff,nTimepoints)}, 1, length(outputs));
    correctedQE = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
    naive_v = feval(corefunc, inputs, outputs, naive_opts);
    naive_shuff = repmat({zeros(opts.shuff, nTimepoints)},1, length(outputs));
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
            naive_shuff{outIdx}(sIdx,:) = shuff_v{outIdx};
        end
    end
    for i=1:opts.xtrp
        randidx = randperm(nTrials, nTrials);
        npartition = [1 2 4];
        part2 = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
        part4 = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
        shuff2 = repmat({zeros(opts.shuff,nTimepoints)}, 1, length(outputs));
        shuff4 = repmat({zeros(opts.shuff,nTimepoints)}, 1, length(outputs));
        inputs_s = inputs;
        idx = repmat({':'}, 1, ndims(inputs{1}));
        idx{end} = randidx;
        for var = 1:length(inputs)
            inputs_s{var} = inputs{var}(idx{:});
        end
        for np=2:length(npartition)
            for pidx = 1:npartition(np)
                inputs_p = partition(inputs_s, npartition(np), pidx, 0);
                value_tmp = feval(corefunc, inputs_p, outputs, naive_opts);
                for outIdx = 1:length(outputs)
                    for t = 1:nTimepoints
                        if npartition(np)==2
                            part2{outIdx}(1,t) = part2{outIdx}(1,t) + value_tmp{outIdx}(1,t)/2;
                        elseif npartition(np)==4
                            part4{outIdx}(1,t) = part4{outIdx}(1,t) + value_tmp{outIdx}(1,t)/4;
                        end
                    end
                end
                % shuffle
                for shuffIdx = 1:opts.shuff
                    inputs_sh = inputs_p;
                    for var = 1:length(inputs_p)
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
                        for t = 1:nTimepoints
                            if npartition(np)==2
                                shuff2{outIdx}(shuffIdx,t) = shuff2{outIdx}(shuffIdx,t) + shuff_v{outIdx}(1,t);
                            elseif npartition(np)==4
                                shuff4{outIdx}(shuffIdx,t) = shuff4{outIdx}(shuffIdx,t) + shuff_v{outIdx}(1,t);
                            end
                        end
                    end
                end
            end
        end
        x_extrap = npartition./nTrials;
        for outIdx = 1:length(outputs)
            for t = 1:nTimepoints
                if strcmp(corr,'qe_shuffSub')
                    p = polyfit(x_extrap, [naive_v{outIdx}(1,t), part2{outIdx}(1,t), part4{outIdx}(1,t)], 2);
                    correctedQE{outIdx}(1,t) = correctedQE{outIdx}(1,t) +(p(3)/xtrp);
                    for shuffIdx = 1:opts.shuff
                        y = [naive_shuff{outIdx}(shuffIdx,t),  shuff2{outIdx}(shuffIdx,t),  shuff4{outIdx}(shuffIdx,t)];
                        p = polyfit(x_extrap, y, 2);
                        correctedSh{outIdx}(shuffIdx,t)  =  correctedSh{outIdx}(shuffIdx,t) + (p(3)/xtrp);
                    end
                elseif strcmp(corr,'le_shuffSub')
                    p = polyfit(x_extrap(1:2), [naive_v{outIdx}(1,t), part2{outIdx}(1,t)], 1);
                    correctedQE{outIdx}(1,t) = correctedQE{outIdx}(1,t) +(p(2)/xtrp);
                    for shuffIdx = 1:opts.shuff
                        y = [naive_shuff{outIdx}(shuffIdx,t),  shuff2{outIdx}(shuffIdx,t)];
                        p = polyfit(x_extrap(1:2), y, 1);
                        correctedSh{outIdx}(shuffIdx,t)  =  correctedSh{outIdx}(t) + (p(2)/xtrp);
                    end
                end
            end
        end
    end
    for t = 1:nTimepoints
        for outIdx = 1:length(outputs)
            corrected_v{outIdx}(1,t) = correctedQE{outIdx}(1,t) - mean(correctedSh{outIdx}(:,t));
        end
    end
elseif strcmp(func2str(corefunc), 'PID')
    corrected_v = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
    correctedSh = repmat({zeros(opts.shuff,nTimepoints)}, 1, length(outputs));
    correctedQE = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
    nSources = length(inputs)-1;
    if ~isfield(opts, 'pid_constrained')
        if nSources == 2
            opts.pid_constrained = true;
        else
            opts.pid_constrained = false;
        end
    end
    % not shuffled
    naive_v = feval(corefunc, inputs, outputs, naive_opts);
    if opts.pid_constrained
        I_naive.I1  = cell2mat(MI({inputs{1}, inputs{end}}, {'I(A;B)'}, naive_opts));
        I_naive.I2  = cell2mat(MI({inputs{2}, inputs{end}}, {'I(A;B)'}, naive_opts));
        I_naive.I12 = cell2mat(MI({cat(1, inputs{1}, inputs{2}), inputs{end}}, {'I(A;B)'}, naive_opts));
        corrected_I1 = zeros(1, nTimepoints);
        corrected_I2 = zeros(1, nTimepoints);
        corrected_I12 = zeros(1, nTimepoints);
        I_corrected.I1QE = zeros(1, nTimepoints);
        I_corrected.I2QE = zeros(1, nTimepoints);
        I_corrected.I12QE = zeros(1, nTimepoints);
        I_corrected.I1Sh = zeros(opts.shuff,nTimepoints);
        I_corrected.I2Sh = zeros(opts.shuff,nTimepoints);
        I_corrected.I12Sh = zeros(opts.shuff,nTimepoints);
    end
    % shuffled
    naive_shuff = repmat({zeros(opts.shuff, nTimepoints)},1, length(outputs));
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
            naive_shuff{outIdx}(sIdx,:) = shuff_v{outIdx};
        end
        if opts.pid_constrained
            I_naive.I1sh  = cell2mat(MI({inputs_sh{1}, inputs_sh{end}}, {'I(A;B)'}, naive_opts));
            I_naive.I2sh  = cell2mat(MI({inputs_sh{2}, inputs_sh{end}}, {'I(A;B)'}, naive_opts));
            I_naive.I12sh = cell2mat(MI({cat(1, inputs_sh{1}, inputs_sh{2}), inputs_sh{end}}, {'I(A;B)'}, naive_opts));
        end
    end
    for i=1:opts.xtrp
        randidx = randperm(nTrials, nTrials);
        npartition = [1 2 4];
        part2 = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
        part4 = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
        shuff2 = repmat({zeros(opts.shuff,nTimepoints)}, 1, length(outputs));
        shuff4 = repmat({zeros(opts.shuff,nTimepoints)}, 1, length(outputs));
        inputs_s = inputs;
        idx = repmat({':'}, 1, ndims(inputs{1}));
        idx{end} = randidx;
        for var = 1:length(inputs)
            inputs_s{var} = inputs{var}(idx{:});
        end
        I_part.I1_2_sh = zeros(opts.shuff, nTimepoints);
        I_part.I2_2_sh = zeros(opts.shuff, nTimepoints);
        I_part.I12_2_sh = zeros(opts.shuff, nTimepoints);
        I_part.I1_4_sh = zeros(opts.shuff, nTimepoints);
        I_part.I2_4_sh = zeros(opts.shuff, nTimepoints);
        I_part.I12_4_sh = zeros(opts.shuff, nTimepoints);
        for np=2:length(npartition)
            for pidx = 1:npartition(np)
                inputs_p = partition(inputs_s, npartition(np), pidx, 0);
                value_tmp = feval(corefunc, inputs_p, outputs, naive_opts);
                for outIdx = 1:length(outputs)
                    for t = 1:nTimepoints
                        if npartition(np)==2
                            part2{outIdx}(1,t) = part2{outIdx}(1,t) + value_tmp{outIdx}(1,t)/2;
                        elseif npartition(np)==4
                            part4{outIdx}(1,t) = part4{outIdx}(1,t) + value_tmp{outIdx}(1,t)/4;
                        end
                    end
                end
                if opts.pid_constrained
                    if npartition(np)==2
                        I_part.I1_2  = cell2mat(MI({inputs_p{1}, inputs_p{end}}, {'I(A;B)'}, naive_opts));
                        I_part.I2_2  = cell2mat(MI({inputs_p{2}, inputs_p{end}}, {'I(A;B)'}, naive_opts));
                        I_part.I12_2 = cell2mat(MI({cat(1, inputs_p{1}, inputs_p{2}), inputs_p{end}}, {'I(A;B)'}, naive_opts));
                    elseif npartition(np)==4
                        I_part.I1_4  = cell2mat(MI({inputs_p{1}, inputs_p{end}}, {'I(A;B)'}, naive_opts));
                        I_part.I2_4  = cell2mat(MI({inputs_p{2}, inputs_p{end}}, {'I(A;B)'}, naive_opts));
                        I_part.I12_4 = cell2mat(MI({cat(1, inputs_p{1}, inputs_p{2}), inputs_p{end}}, {'I(A;B)'}, naive_opts));
                    end
                end
                % shuffle
                for shuffIdx = 1:opts.shuff
                    inputs_sh = inputs_p;
                    for var = 1:length(inputs_p)
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
                        if npartition(np)==2
                            shuff2{outIdx}(shuffIdx,:)   =  shuff2{outIdx}(shuffIdx,:) + shuff_v{outIdx}/2;
                        elseif npartition(np)==4
                            shuff4{outIdx}(shuffIdx,:)   =  shuff4{outIdx}(shuffIdx,:) + shuff_v{outIdx}/4;
                        end

                    end
                    if opts.pid_constrained
                        if npartition(np)==2
                            I_part.I1_2_sh(shuffIdx,:)  = I_part.I1_2_sh(shuffIdx,:)  + cell2mat(MI({inputs_sh{1}, inputs_sh{end}}, {'I(A;B)'}, naive_opts))/2;
                            I_part.I2_2_sh(shuffIdx,:)  = I_part.I2_2_sh(shuffIdx,:)  + cell2mat(MI({inputs_sh{2}, inputs_sh{end}}, {'I(A;B)'}, naive_opts))/2;
                            I_part.I12_2_sh(shuffIdx,:) = I_part.I12_2_sh(shuffIdx,:) + cell2mat(MI({cat(1, inputs_sh{1}, inputs_sh{2}), inputs_sh{end}}, {'I(A;B)'}, naive_opts))/2;
                        elseif npartition(np)==4
                            I_part.I1_4_sh(shuffIdx,:)  = I_part.I1_4_sh(shuffIdx,:)  + cell2mat(MI({inputs_sh{1}, inputs_sh{end}}, {'I(A;B)'}, naive_opts))/4;
                            I_part.I2_4_sh(shuffIdx,:)  = I_part.I2_4_sh(shuffIdx,:)  + cell2mat(MI({inputs_sh{2}, inputs_sh{end}}, {'I(A;B)'}, naive_opts))/4;
                            I_part.I12_4_sh(shuffIdx,:) = I_part.I12_4_sh(shuffIdx,:) + cell2mat(MI({cat(1, inputs_sh{1}, inputs_sh{2}), inputs_sh{end}}, {'I(A;B)'}, naive_opts))/4;
                        end
                    end
                end
            end
        end
        x_extrap = npartition./nTrials;
        for outIdx = 1:length(outputs)
            for t = 1:nTimepoints
                if strcmp(corr,'qe_shuffSub')
                    p = polyfit(x_extrap, [naive_v{outIdx}(1,t), part2{outIdx}(1,t), part4{outIdx}(1,t)], 2);
                    correctedQE{outIdx}(1,t) = correctedQE{outIdx}(1,t) +(p(3)/xtrp);
                    for shuffIdx = 1:opts.shuff
                        y = [naive_shuff{outIdx}(shuffIdx,t),  shuff2{outIdx}(shuffIdx,t),  shuff4{outIdx}(shuffIdx,t)];
                        p = polyfit(x_extrap, y, 2);
                        correctedSh{outIdx}(shuffIdx,t)  =  correctedSh{outIdx}(shuffIdx,t) + p(3)/xtrp;
                    end
                elseif strcmp(corr,'le_shuffSub')
                    p = polyfit(x_extrap(1:2), [naive_v{outIdx}(1,t), part2{outIdx}(1,t)], 1);
                    correctedQE{outIdx}(1,t) = correctedQE{outIdx}(1,t) +(p(2)/xtrp);
                    for shuffIdx = 1:opts.shuff
                        y = [naive_shuff{outIdx}(shuffIdx,t),  shuff2{outIdx}(shuffIdx,t)];
                        p = polyfit(x_extrap(1:2), y, 1);
                        correctedSh{outIdx}(shuffIdx,t)  =  correctedSh{outIdx}(t) + p(2)/xtrp;
                    end
                end
            end
        end
        if opts.pid_constrained
            if strcmp(opts.bias,'qe_shuffSub')
                for t = 1:nTimepoints
                    y = [I_naive.I1(t), I_part.I1_2(t), I_part.I1_4(t)];
                    p = polyfit(x_extrap, y, 2);
                    I_corrected.I1QE(t) = I_corrected.I1QE(t) + p(3)/xtrp;
                    y = [I_naive.I2(t), I_part.I2_2(t), I_part.I2_4(t)];
                    p = polyfit(x_extrap, y, 2);
                    I_corrected.I2QE(t) = I_corrected.I2QE(t) + p(3)/xtrp;
                    y = [I_naive.I12(t), I_part.I12_2(t), I_part.I12_4(t)];
                    p = polyfit(x_extrap, y, 2);
                    I_corrected.I12QE(t) = I_corrected.I12QE(t) + p(3)/xtrp;
                    for shuffIdx = 1:opts.shuff
                        y = [I_naive.I1sh(t), I_part.I1_2_sh(t), I_part.I1_4_sh(t)];
                        p = polyfit(x_extrap, y, 2);
                        I_corrected.I1Sh(shuffIdx,t) =  I_corrected.I1Sh(shuffIdx,t) + p(3)/xtrp;
                        y = [I_naive.I2sh(t), I_part.I2_2_sh(t), I_part.I2_4_sh(t)];
                        p = polyfit(x_extrap, y, 2);
                        I_corrected.I2Sh(shuffIdx,t) = I_corrected.I2Sh(shuffIdx,t) + p(3)/xtrp;
                        y = [I_naive.I12sh(t), I_part.I12_2_sh(t), I_part.I12_4_sh(t)];
                        p = polyfit(x_extrap, y, 2);
                        I_corrected.I12Sh(shuffIdx,t) = I_corrected.I12Sh(shuffIdx,t) + p(3)/xtrp;
                    end
                end
            elseif strcmp(opts.bias,'le_shuffSub')
                for t = 1:nTimepoints
                    y = [I_naive.I1(t), I_part.I1_2(t)];
                    p = polyfit(x_extrap(1:2), y, 1);
                    I_corrected.I1QE(t) = I_corrected.I1QE(t) + p(2)/xtrp;
                    y = [I_naive.I2(t), I_part.I2_2(t)];
                    p = polyfit(x_extrap(1:2), y, 1);
                    I_corrected.I2QE(t) = I_corrected.I2QE(t) + p(2)/xtrp;
                    y = [I_naive.I12(t), I_part.I12_2(t)];
                    p = polyfit(x_extrap(1:2), y, 1);
                    I_corrected.I12QE(t) = I_corrected.I12QE(t) + p(2)/xtrp;
                    for shuffIdx = 1:opts.shuff
                        y = [I_naive.I1sh(t), I_part.I1_2_sh(t)];
                        p = polyfit(x_extrap(1:2), y, 1);
                        I_corrected.I1Sh(shuffIdx,t) =  I_corrected.I1Sh(shuffIdx,t) + p(2)/xtrp;
                        y = [I_naive.I2sh(t), I_part.I2_2_sh(t)];
                        p = polyfit(x_extrap(1:2), y, 1);
                        I_corrected.I2Sh(shuffIdx,t) = I_corrected.I2Sh(shuffIdx,t) + p(2)/xtrp;
                        y = [I_naive.I12sh(t), I_part.I12_2_sh(t)];
                        p = polyfit(x_extrap(1:2), y, 1);
                        I_corrected.I12Sh(shuffIdx,t) = I_corrected.I12Sh(shuffIdx,t) + p(2)/xtrp;
                    end
                end
            end
        end
    end
    for t = 1:nTimepoints
        for outIdx = 1:length(outputs)
            corrected_v{outIdx}(1,t) = correctedQE{outIdx}(1,t) - mean(correctedSh{outIdx}(:,t));
        end
        if opts.pid_constrained
            corrected_I1(1,t) = I_corrected.I1QE(t) - mean(I_corrected.I1Sh(:,t));
            corrected_I2(1,t) = I_corrected.I2QE(t) - mean(I_corrected.I2Sh(:,t));
            corrected_I12(1,t) = I_corrected.I12QE(t) - mean(I_corrected.I12Sh(:,t));
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
                syn = corrected_I12-corrected_I1-corrected_I2+red;
                unq1 = corrected_I1-red;
                unq2 = corrected_I2-red;
                %naive
                red_naive = naive_v{pos};
                syn_naive  = I_naive.I12-I_naive.I1-I_naive.I2+red_naive;
                unq1_naive  = I_naive.I1-red_naive;
                unq2_naive  = I_naive.I2-red_naive;
            case 'Unq1'
                unq1 = corrected_v{pos};
                red = corrected_I1-unq1;
                syn =  corrected_I12-corrected_I2-unq1;
                unq2 = corrected_I2-corrected_I1+unq1;
                %naive
                unq1_naive = naive_v{pos};
                red_naive = I_naive.I1-unq1_naive;
                syn_naive =  I_naive.I12-I_naive.I2-unq1_naive;
                unq2_naive = I_naive.I2-I_naive.I1+unq1_naive;
            case 'Unq2'
                unq2 = corrected_v{pos};
                red = corrected_I2-unq2;
                syn =  corrected_I12-corrected_I1-unq2;
                unq1 = corrected_I1-corrected_I2+unq2;
                %naive
                unq2_naive = naive_v{pos};
                red_naive = I_naive.I2-unq2_naive;
                syn_naive =  I_naive.I12-I_naive.I1-unq2_naive;
                unq1_naive = I_naive.I1-I_naive.I2+unq2_naive;
            case 'Syn'
                syn = corrected_v{pos};
                red = corrected_I1+corrected_I2-corrected_I12+syn;
                unq1 = corrected_I12-corrected_I2-syn;
                unq2 = corrected_I12-corrected_I1-syn;
                %naive
                syn_naive  = naive_v{pos};
                red_naive  = I_naive.I1+I_naive.I2-I_naive.I12+syn_naive;
                unq1_naive  = I_naive.I12-I_naive.I2-syn_naive;
                unq2_naive  = I_naive.I12-I_naive.I1-syn_naive;
        end
        for i = 1:length(outputs)
            switch outputs{i}
                case 'Syn'
                    corrected_v{i} = syn;
                    naive_v{i}       = syn_naive;
                case 'Red'
                    corrected_v{i} = red;
                    naive_v{i}       = red_naive;
                case 'Unq1'
                    corrected_v{i} = unq1;
                    naive_v{i}       = unq1_naive;
                case 'Unq2'
                    corrected_v{i} = unq2;
                    naive_v{i}       = unq2_naive;
                case 'Unq'
                    corrected_v{i} = unq1 + unq2;
                    naive_v{i}       = unq1_naive + unq2_naive;
                case 'Joint'
                    corrected_v{i} = red+syn+unq1+unq2;
                    naive_v{i}       = red_naive+syn_naive+unq1_naive+unq2_naive;
                case 'Union'
                    corrected_v{i} = (red+syn+unq1+unq2) - syn;
                    naive_v{i}       = (red_naive+syn_naive+unq1_naive+unq2_naive) - syn_naive;
                otherwise
                    corrected_v{i} = NaN;
                    naive_v{i}       = NaN;
            end
        end
    end
elseif strcmp(func2str(corefunc), 'cFIT') || strcmp(func2str(corefunc), 'FIT')
    [~,~,~,atom1_naive, atom2_naive] = feval(corefunc, inputs, outputs, naive_opts);
    numAtoms = length(atom1_naive);
    naiveSh.atom1 = repmat({zeros(numAtoms,opts.shuff)},1, length(outputs));
    naiveSh.atom2 = repmat({zeros(numAtoms,opts.shuff)},1, length(outputs));
    correctedSh.atom1 = repmat({zeros(numAtoms,opts.shuff)}, 1, length(outputs));
    correctedQE.atom1 = repmat({zeros(1,numAtoms)}, 1, length(outputs));
    correctedSh.atom2 = repmat({zeros(numAtoms,opts.shuff)}, 1, length(outputs));
    correctedQE.atom2 = repmat({zeros(1,numAtoms)}, 1, length(outputs));
    corrected_atom.atom1 = repmat({zeros(1,numAtoms)}, 1, length(outputs));
    corrected_atom.atom2 = repmat({zeros(1,numAtoms)}, 1, length(outputs));
    corrected_v = repmat({zeros(1)}, 1, length(outputs));
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
            for idx = 1:length(atom1_naive)
                naiveSh.atom1{outIdx}(idx,sIdx) = atom1_shuff{outIdx}(idx);
                naiveSh.atom2{outIdx}(idx,sIdx) = atom2_shuff{outIdx}(idx);
            end 
        end
    end
    for i=1:opts.xtrp
        randidx = randperm(nTrials, nTrials);
        npartition = [1 2 4];
        part2.atom1 = repmat({zeros(1,numAtoms)}, 1, length(outputs));
        part4.atom1 = repmat({zeros(1,numAtoms)}, 1, length(outputs));
        part2.atom2 = repmat({zeros(1,numAtoms)}, 1, length(outputs));
        part4.atom2 = repmat({zeros(1,numAtoms)}, 1, length(outputs));
        shuff2.atom1 = repmat({zeros(numAtoms,opts.shuff)}, 1, length(outputs));
        shuff4.atom1 = repmat({zeros(numAtoms,opts.shuff)}, 1, length(outputs));
        shuff2.atom2 = repmat({zeros(numAtoms,opts.shuff)}, 1, length(outputs));
        shuff4.atom2 = repmat({zeros(numAtoms,opts.shuff)}, 1, length(outputs));
        inputs_s = inputs;
        idx = repmat({':'}, 1, ndims(inputs{1}));
        idx{end} = randidx;
        for var = 1:length(inputs)
            inputs_s{var} = inputs{var}(idx{:});
        end
        for np=2:length(npartition)
            for pidx = 1:npartition(np)
                inputs_p = partition(inputs_s, npartition(np), pidx, 0);
                [~,~,~,atom1_naive, atom2_naive] = feval(corefunc, inputs_p, outputs, naive_opts);
                for idx = 1:length(atom1_naive)
                    for outIdx = 1:length(outputs)
                        if npartition(np)==2
                            part2.atom1{outIdx}(idx) = part2.atom1{outIdx}(idx) + atom1_naive{outIdx}(idx)/2;
                            part2.atom2{outIdx}(idx) = part2.atom2{outIdx}(idx) + atom2_naive{outIdx}(idx)/2;
                        elseif npartition(np)==4
                            part4.atom1{outIdx}(idx) = part4.atom1{outIdx}(idx) + atom1_naive{outIdx}(idx)/4;
                            part4.atom2{outIdx}(idx) = part4.atom2{outIdx}(idx) + atom2_naive{outIdx}(idx)/4;
                        end
                    end
                end
                % shuffle
                for shuffIdx = 1:opts.shuff
                    inputs_sh = inputs_p;
                    for var = 1:length(inputs_p)
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
                        for idx = 1:length(atom1_naive)
                            if npartition(np)==2
                                shuff2.atom1{outIdx}(idx,shuffIdx) = shuff2.atom1{outIdx}(idx,shuffIdx) + atom1_shuff{outIdx}(idx)/2;
                                shuff2.atom2{outIdx}(idx,shuffIdx) = shuff2.atom2{outIdx}(idx,shuffIdx) + atom2_shuff{outIdx}(idx)/2;
                            elseif npartition(np)==4
                                shuff4.atom1{outIdx}(idx,shuffIdx) = shuff4.atom1{outIdx}(idx,shuffIdx) + atom1_shuff{outIdx}(idx)/4;
                                shuff4.atom2{outIdx}(idx,shuffIdx) = shuff4.atom2{outIdx}(idx,shuffIdx) + atom2_shuff{outIdx}(idx)/4;
                            end
                        end
                    end
                end
            end
        end
        x_extrap = npartition./nTrials;
        for outIdx = 1:length(outputs)
            if strcmp(corr,'qe_shuffSub')
                for idx = 1:length(atom1_naive)
                    p = polyfit(x_extrap, [atom1_naive{outIdx}(idx), part2.atom1{outIdx}(idx), part4.atom1{outIdx}(idx)], 2);
                    correctedQE.atom1{outIdx}(idx) = correctedQE.atom1{outIdx}(idx) +(p(3)/xtrp);
                    p = polyfit(x_extrap, [atom2_naive{outIdx}(idx), part2.atom2{outIdx}(idx), part4.atom2{outIdx}(idx)], 2);
                    correctedQE.atom2{outIdx}(idx) = correctedQE.atom2{outIdx}(idx) +(p(3)/xtrp);
                    for shuffIdx = 1:opts.shuff
                        y = [naiveSh.atom1{outIdx}(idx,shuffIdx),  shuff2.atom1{outIdx}(idx,shuffIdx),  shuff4.atom2{outIdx}(idx,shuffIdx)];
                        p = polyfit(x_extrap, y, 2);
                        correctedSh.atom1{outIdx}(idx,shuffIdx)  =  correctedSh.atom1{outIdx}(idx,shuffIdx) + (p(3)/xtrp);
                        y = [naiveSh.atom2{outIdx}(idx,shuffIdx),  shuff2.atom2{outIdx}(idx,shuffIdx),  shuff4.atom2{outIdx}(idx,shuffIdx)];
                        p = polyfit(x_extrap, y, 2);
                        correctedSh.atom2{outIdx}(idx,shuffIdx)  =  correctedSh.atom2{outIdx}(idx,shuffIdx) + (p(3)/xtrp);
                    end
                end
            elseif strcmp(corr,'le_shuffSub')
                for idx = 1:length(atom1_naive)
                    p = polyfit(x_extrap(1:2), [atom1_naive{outIdx}(idx), part2.atom1{outIdx}(idx)], 1);
                    correctedQE.atom1{outIdx}(idx) = correctedQE.atom1{outIdx}(idx) +(p(2)/xtrp);
                    p = polyfit(x_extrap(1:2), [atom2_naive{outIdx}(idx), part2.atom2{outIdx}(idx)], 1);
                    correctedQE.atom2{outIdx}(idx) = correctedQE.atom2{outIdx}(idx) +(p(2)/xtrp);
                    for shuffIdx = 1:opts.shuff
                        y = [naiveSh.atom1{outIdx}(idx,shuffIdx), shuff2.atom1{outIdx}(idx,shuffIdx)];
                        p = polyfit(x_extrap(1:2), y, 1);
                        correctedSh.atom1{outIdx}(idx,shuffIdx)  =  correctedSh.atom1{outIdx}(idx,shuffIdx) + (p(2)/xtrp);
                        y = [naiveSh.atom2{outIdx}(idx,shuffIdx), shuff2.atom2{outIdx}(idx,shuffIdx)];
                        p = polyfit(x_extrap(1:2), y, 1);
                        correctedSh.atom2{outIdx}(idx,shuffIdx)  =  correctedSh.atom2{outIdx}(idx,shuffIdx) + (p(2)/xtrp);
                    end
                end
            end
        end
        if strcmp(func2str(corefunc), 'FIT')
            for outIdx = 1:length(outputs)
                corrected_atom.atom1{outIdx} = correctedQE.atom1{outIdx} - mean(correctedSh.atom1{outIdx});
                corrected_atom.atom2{outIdx} = correctedQE.atom2{outIdx} - mean(correctedSh.atom2{outIdx});
                corrected_v{outIdx} = min(corrected_atom.atom1{outIdx},corrected_atom.atom2{outIdx});
                naive_v_all{outIdx} = min(atom1_naive{outIdx},atom2_naive{outIdx});
            end
        elseif strcmp(func2str(corefunc), 'cFIT') 
            for outIdx = 1:length(outputs)
                corrected_atom.atom1{outIdx} = correctedQE.atom1{outIdx} - mean(correctedSh.atom1{outIdx});
                corrected_atom.atom2{outIdx} = correctedQE.atom2{outIdx} - mean(correctedSh.atom2{outIdx});
                atom1 = min(corrected_atom.atom1{outIdx});
                atom2 = min(corrected_atom.atom2{outIdx});
                atom1_n = min(atom1_naive{outIdx});
                atom2_n = min(atom2_naive{outIdx});
                corrected_v{outIdx} = atom1-atom2; 
                naive_v_all{outIdx} = atom1_n-atom2_n; 
            end
        end 
    end
    naive_v = naive_v_all;
elseif strcmp(func2str(corefunc), 'II')
    [~,~,~,atom1_naive, atom2_naive] = feval(corefunc, inputs, outputs, naive_opts);
    naiveSh.atom1 = repmat({zeros(opts.shuff,nTimepoints)},1, length(outputs));
    naiveSh.atom2 = repmat({zeros(opts.shuff,nTimepoints)},1, length(outputs));
    correctedSh.atom1 = repmat({zeros(opts.shuff,nTimepoints)}, 1, length(outputs));
    correctedQE.atom1 = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
    correctedSh.atom2 = repmat({zeros(opts.shuff,nTimepoints)}, 1, length(outputs));
    correctedQE.atom2 = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
    corrected_atom.atom1 = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
    corrected_atom.atom2 = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
    corrected_v = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
    
    for i=1:opts.xtrp
        randidx = randperm(nTrials, nTrials);
        npartition = [1 2 4];
        part2.atom1 = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
        part4.atom1 = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
        part2.atom2 = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
        part4.atom2 = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
        shuff2.atom1 = repmat({zeros(opts.shuff,nTimepoints)}, 1, length(outputs));
        shuff4.atom1 = repmat({zeros(opts.shuff,nTimepoints)}, 1, length(outputs));
        shuff2.atom2 = repmat({zeros(opts.shuff,nTimepoints)}, 1, length(outputs));
        shuff4.atom2 = repmat({zeros(opts.shuff,nTimepoints)}, 1, length(outputs));

        inputs_s = inputs;
        idx = repmat({':'}, 1, ndims(inputs{1}));
        idx{end} = randidx;
        for var = 1:length(inputs)
            inputs_s{var} = inputs{var}(idx{:});
        end
        for np=2:length(npartition)
            for pidx = 1:npartition(np)
                inputs_p = partition(inputs_s, npartition(np), pidx,0);
                [~,~,~,atom1_tmp, atom2_tmp] = feval(corefunc, inputs_p,outputs, naive_opts);
                for outIdx = 1:length(outputs)
                    if npartition(np)==2
                        part2.atom1{outIdx} = part2.atom1{outIdx} + atom1_naive{outIdx}/2;
                        part2.atom2{outIdx} = part2.atom2{outIdx} + atom2_naive{outIdx}/2;
                    elseif npartition(np)==4
                        part4.atom1{outIdx} = part4.atom1{outIdx} + atom1_naive{outIdx}/4;
                        part4.atom2{outIdx} = part4.atom2{outIdx} + atom2_naive{outIdx}/4;
                    end
                end
                % shuffle
                for shuffIdx = 1:opts.shuff
                    inputs_sh = inputs_p;
                    for var = 1:length(inputs_p)
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
                        if npartition(np)==2
                            shuff2.atom1{outIdx}(shuffIdx,:) = shuff2.atom1{outIdx}(shuffIdx,:) + atom1_shuff{outIdx}/2;
                            shuff2.atom2{outIdx}(shuffIdx,:) = shuff2.atom2{outIdx}(shuffIdx,:) + atom2_shuff{outIdx}/2;
                        elseif npartition(np)==4
                            shuff4.atom1{outIdx}(shuffIdx,:) = shuff4.atom1{outIdx}(shuffIdx,:) + atom1_shuff{outIdx}/4;
                            shuff4.atom2{outIdx}(shuffIdx,:) = shuff4.atom2{outIdx}(shuffIdx,:) + atom2_shuff{outIdx}/4;
                        end

                    end
                end
            end
        end
        x_extrap = npartition./nTrials;
         for outIdx = 1:length(outputs)
            if strcmp(corr,'qe_shuffSub')
                for t = 1:nTimepoints
                    p = polyfit(x_extrap, [atom1_naive{outIdx}(t), part2.atom1{outIdx}(t), part4.atom1{outIdx}(t)], 2);
                    correctedQE.atom1{outIdx}(t) = correctedQE.atom1{outIdx}(t) +(p(3)/xtrp);
                    p = polyfit(x_extrap, [atom2_naive{outIdx}(t), part2.atom2{outIdx}(t), part4.atom2{outIdx}(t)], 2);
                    correctedQE.atom2{outIdx}(t) = correctedQE.atom2{outIdx}(t) +(p(3)/xtrp);
                    for shuffIdx = 1:opts.shuff
                        y = [naiveSh.atom1{outIdx}(shuffIdx,t),  shuff2.atom1{outIdx}(shuffIdx,t),  shuff4.atom2{outIdx}(shuffIdx,t)];
                        p = polyfit(x_extrap, y, 2);
                        correctedSh.atom1{outIdx}(shuffIdx,t)  =  correctedSh.atom1{outIdx}(shuffIdx,t) + (p(3)/xtrp);
                        y = [naiveSh.atom2{outIdx}(shuffIdx,t),  shuff2.atom2{outIdx}(shuffIdx,t),  shuff4.atom2{outIdx}(shuffIdx,t)];
                        p = polyfit(x_extrap, y, 2);
                        correctedSh.atom2{outIdx}(shuffIdx,t)  =  correctedSh.atom2{outIdx}(shuffIdx,t) + (p(3)/xtrp);
                    end
                end
            elseif strcmp(corr,'le_shuffSub')
                for t = 1:nTimepoints
                    p = polyfit(x_extrap(1:2), [atom1_naive{outIdx}(t), part2.atom1{outIdx}(t)], 1);
                    correctedQE.atom1{outIdx}(t) = correctedQE.atom1{outIdx}(t) +(p(2)/xtrp);
                    p = polyfit(x_extrap(1:2), [atom2_naive{outIdx}(t), part2.atom2{outIdx}(t)], 1);
                    correctedQE.atom2{outIdx}(t) = correctedQE.atom2{outIdx}(t) +(p(2)/xtrp);
                    for shuffIdx = 1:opts.shuff
                        y = [naiveSh.atom1{outIdx}(shuffIdx,t), shuff2.atom1{outIdx}(shuffIdx,t)];
                        p = polyfit(x_extrap(1:2), y, 1);
                        correctedSh.atom1{outIdx}(shuffIdx,t)  =  correctedSh.atom1{outIdx}(shuffIdx,t) + (p(2)/xtrp);
                        y = [naiveSh.atom2{outIdx}(shuffIdx,t), shuff2.atom2{outIdx}(shuffIdx,t)];
                        p = polyfit(x_extrap(1:2), y, 1);
                        correctedSh.atom2{outIdx}(shuffIdx,t)  =  correctedSh.atom2{outIdx}(shuffIdx,t) + (p(2)/xtrp);
                    end
                end
            end
         end
    end
    for outIdx = 1:length(outputs)
        corrected_atom.atom1{outIdx} = correctedQE.atom1{outIdx} - mean(correctedSh.atom1{outIdx},1);
        corrected_atom.atom2{outIdx} = correctedQE.atom2{outIdx} - mean(correctedSh.atom2{outIdx},1);
        corrected_v{outIdx} = min(corrected_atom.atom1{outIdx},corrected_atom.atom2{outIdx});
        naive_v_all{outIdx} = min(atom1_naive{outIdx},atom2_naive{outIdx});
    end
    naive_v = naive_v_all;
end



