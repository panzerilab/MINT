function [corrected_v, naive_v, shuff_all] = extrapolation(inputs, outputs, corr, corefunc, varargin)
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
shuff_all = cell(nTimepoints, length(outputs));
corrected_v = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
xtrp = opts.xtrp;

if strcmp(func2str(corefunc), 'PID')
    naive_opts = opts;
    naive_opts.bias = 'naive';
    nSources = length(inputs)-1;

    if ~isfield(opts, 'pid_constrained')
        if nSources == 2
            opts.pid_constrained = true;
        else
            opts.pid_constrained = false;
        end
    end
    PID_naive = feval(corefunc, inputs, outputs, naive_opts);
    PID_corrected = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
    if opts.pid_constrained
        I1_naive  = cell2mat(MI({inputs{1}, inputs{end}}, {'I(A;B)'}, naive_opts));
        I2_naive  = cell2mat(MI({inputs{2}, inputs{end}}, {'I(A;B)'}, naive_opts));
        I12_naive = cell2mat(MI({cat(1, inputs{1}, inputs{2}), inputs{end}}, {'I(A;B)'}, naive_opts));
        I1_corrected = zeros(1, nTimepoints);
        I2_corrected = zeros(1, nTimepoints);
        I12_corrected = zeros(1, nTimepoints);
    end
    nTrials = size(inputs{1},length(size(inputs{1})));
    if opts.parallel == false || ((opts.parallel == true) && (opts.shuff > 0))
        for i=1:opts.xtrp
            randidx = randperm(nTrials, nTrials);
            npartition = [1 2 4];
            part2 =  repmat({zeros(1,nTimepoints)}, 1, length(outputs));
            part4 =  repmat({zeros(1,nTimepoints)}, 1, length(outputs));
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
                    inputs_p = partition(inputs_s, npartition(np), pidx,0);

                    if strcmp(corr, 'qe_shuffSub') || strcmp(corr, 'le_shuffSub')
                        PID_value_tmp = shuffle_subtraction(inputs_p, outputs, corefunc, naive_opts);
                        if opts.pid_constrained
                            I1_tmp = shuffle_subtraction({inputs_p{1}, inputs_p{end}}, {'I(A;B)'}, @MI, naive_opts);
                            I2_tmp = shuffle_subtraction({inputs_p{2}, inputs_p{end}}, {'I(A;B)'}, @MI, naive_opts);
                            I12_tmp = shuffle_subtraction({cat(1, inputs_p{1}, inputs_p{2}), inputs_p{end}}, {'I(A;B)'}, @MI, naive_opts);
                        end
                    else
                        PID_value_tmp = feval(corefunc, inputs_p, outputs, naive_opts);
                        if opts.pid_constrained
                            I1_tmp = cell2mat(MI({inputs_p{1}, inputs_p{end}}, {'I(A;B)'}, naive_opts));
                            I2_tmp = cell2mat(MI({inputs_p{2}, inputs_p{end}}, {'I(A;B)'}, naive_opts));
                            I12_tmp = cell2mat(MI({cat(1, inputs_p{1}, inputs_p{2}), inputs_p{end}}, {'I(A;B)'}, naive_opts));
                        end
                    end
                    for outIdx = 1:length(outputs)
                        if npartition(np)==2
                            part2{outIdx} = part2{outIdx} + PID_value_tmp{outIdx}/2;
                        elseif npartition(np)==4
                            part4{outIdx} = part4{outIdx} + PID_value_tmp{outIdx}/4;
                        end
                    end
                    if opts.pid_constrained
                        if npartition(np)==2
                            I1_2 = I1_2 + I1_tmp/2;
                            I2_2 = I2_2 + I2_tmp/2;
                            I12_2 = I12_2 + I12_tmp/2;
                        elseif npartition(np)==4
                            I1_4 = I1_4 + I1_tmp/4;
                            I2_4 = I2_4 + I2_tmp/4;
                            I12_4 = I12_4 + I12_tmp/4;
                        end
                    end

                end
            end
            x_extrap = npartition./nTrials;
            for outIdx = 1:length(outputs)
                if strcmp(opts.bias,'qe') ||strcmp(opts.bias,'qe_shuffSub')
                    for t = 1:nTimepoints
                        y = [PID_naive{outIdx}(t), part2{outIdx}(t), part4{outIdx}(t)];
                        p = polyfit(x_extrap, y, 2); 
                        PID_corrected{outIdx}(1,t)  =  PID_corrected{outIdx}(t) + p(3)/xtrp;
                    end
                elseif strcmp(opts.bias,'le')||strcmp(opts.bias,'le_shuffSub')
                    for t = 1:nTimepoints
                        y = [PID_naive{outIdx}(t), part2{outIdx}(t)];
                        p = polyfit(x_extrap, y, 1); 
                        PID_corrected{outIdx}(1,t) = PID_corrected{outIdx}(t) + p(2)/xtrp;
                    end
                end
            end
            if opts.pid_constrained
                if strcmp(opts.bias,'qe') ||strcmp(opts.bias,'qe_shuffSub')
                    for t = 1:nTimepoints
                        y = [I1_naive(t), I1_2(t), I1_4(t)];
                        p = polyfit(x_extrap, y, 2);
                        I1_corrected(t) = I1_corrected(t) + p(3)/xtrp;
                        y = [I2_naive(t), I2_2(t), I2_4(t)];
                        p = polyfit(x_extrap, y, 2);
                        I2_corrected(t) = I2_corrected(t) + p(3)/xtrp;                      
                        y = [I12_naive(t), I12_2(t), I12_4(t)];
                        p = polyfit(x_extrap, y, 2);
                        I12_corrected(t) = I12_corrected(t) + p(3)/xtrp;
                    end
                elseif strcmp(opts.bias,'le')||strcmp(opts.bias,'le_shuffSub')
                    for t = 1:nTimepoints
                        y = [I1_naive(t), I1_2(t)];
                        p = polyfit(x_extrap, y, 1);
                        I1_corrected(t) = I1_corrected(t) +  p(2)/xtrp;
                        y = [I2_naive(t), I2_2(t)];
                        p = polyfit(x_extrap, y, 1);
                        I2_corrected(t) = I2_corrected(t) +  p(2)/xtrp;
                        y = [I12_naive(t), I12_2(t)];
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
    atom1_corr = repmat({0}, 1, length(outputs));
    atom2_corr = repmat({0}, 1, length(outputs));
    naive_opts = opts;
    naive_opts.bias = 'naive';
    naive_opts.recall = true;
    [~,~,~,atom1_naive, atom2_naive] = feval(corefunc, inputs, outputs, naive_opts);
    nTrials = size(inputs{1},length(size(inputs{1})));
    if opts.parallel == false || ((opts.parallel == true) && (opts.shuff > 0))
        for i=1:opts.xtrp
            randidx = randperm(nTrials, nTrials);
            npartition = [1 2 4];
            atom1_2 =  repmat({0}, 1, length(outputs));
            atom1_4 =  repmat({0}, 1, length(outputs));
            atom2_2 =  repmat({0}, 1, length(outputs));
            atom2_4 =  repmat({0}, 1, length(outputs));

            inputs_s = inputs;
            idx = repmat({':'}, 1, ndims(inputs{1}));
            idx{end} = randidx;
            for var = 1:length(inputs)
                inputs_s{var} = inputs{var}(idx{:});
            end
            for np=2:length(npartition)
                for pidx = 1:npartition(np)
                    inputs_p = partition(inputs_s, npartition(np), pidx,0);
                    if strcmp(corr, 'qe_shuffSub') || strcmp(corr, 'le_shuffSub')
                        [~, ~, ~, addOut] = shuffle_subtraction(inputs_p, outputs, corefunc, naive_opts);
                        for outIdx = 1:length(outputs)
                            atom1_tmp{outIdx} = addOut{1,outIdx}(end);
                            atom2_tmp{outIdx} = addOut{2,outIdx}(end);
                        end
                    else
                        [~,~,~,atom1_tmp, atom2_tmp] = feval(corefunc, inputs_p,outputs, naive_opts);
                    end
                    for outIdx = 1:length(outputs)
                        if npartition(np)==2
                            atom1_2{outIdx} = atom1_2{outIdx}+ atom1_tmp{outIdx}/2;
                            atom2_2{outIdx} = atom2_2{outIdx} + atom2_tmp{outIdx}/2;
                        elseif npartition(np)==4
                            atom1_4{outIdx} = atom1_4{outIdx} + atom1_tmp{outIdx}/4;
                            atom2_4{outIdx} = atom2_4{outIdx} + atom2_tmp{outIdx}/4;
                        end
                    end
                end
            end
            x_extrap = npartition./nTrials;
            if strcmp(func2str(corefunc), 'cFIT')
                for outIdx = 1:length(outputs)
                    for idx = 1:length(atom1_naive)
                        if strcmp(opts.bias,'qe') ||strcmp(opts.bias,'qe_shuffSub')
                            atom1_corr_tmp = polyfit(x_extrap, [atom1_naive{outIdx}(idx), atom1_2{outIdx}(idx), atom1_4{outIdx}(idx)], 2);
                            atom2_corr_tmp = polyfit(x_extrap, [atom2_naive{outIdx}(idx), atom2_2{outIdx}(idx), atom2_4{outIdx}(idx)], 2);
                            atom1_corr{outIdx}(idx) = atom1_corr{outIdx}(idx) + (atom1_corr_tmp(3)/xtrp);
                            atom2_corr{outIdx}(idx) = atom2_corr{outIdx}(idx) + (atom2_corr_tmp(3)/xtrp);
                        elseif strcmp(opts.bias,'le')||strcmp(opts.bias,'le_shuffSub')
                            atom1_corr_tmp = polyfit(x_extrap(1:2), [atom1_naive{outIdx}(idx), atom1_2{outIdx}(idx)], 1);
                            atom2_corr_tmp = polyfit(x_extrap(1:2), [atom2_naive{outIdx}(idx), atom2_2{outIdx}(idx)], 1);
                            atom1_corr{outIdx}(idx) = atom1_corr{outIdx}(idx) + (atom1_corr_tmp(2)/xtrp);
                            atom2_corr{outIdx}(idx) = atom2_corr{outIdx}(idx) + (atom2_corr_tmp(2)/xtrp);
                        end
                    end
                end
            else
                for outIdx = 1:length(outputs)
                    if strcmp(opts.bias,'qe') ||strcmp(opts.bias,'qe_shuffSub')
                        atom1_corr_tmp = polyfit(x_extrap, [atom1_naive{outIdx}, atom1_2{outIdx}, atom1_4{outIdx}], 2);
                        atom2_corr_tmp = polyfit(x_extrap, [atom2_naive{outIdx}, atom2_2{outIdx}, atom2_4{outIdx}], 2);
                        atom1_corr{outIdx} = atom1_corr{outIdx} + (atom1_corr_tmp(3)/xtrp);
                        atom2_corr{outIdx} = atom2_corr{outIdx} + (atom2_corr_tmp(3)/xtrp);
                    elseif strcmp(opts.bias,'le')||strcmp(opts.bias,'le_shuffSub')
                        atom1_corr_tmp = polyfit(x_extrap(1:2), [atom1_naive{outIdx}, atom1_2{outIdx}], 1);
                        atom2_corr_tmp = polyfit(x_extrap(1:2), [atom2_naive{outIdx}, atom2_2{outIdx}], 1);
                        atom1_corr{outIdx} = atom1_corr{outIdx} + (atom1_corr_tmp(2)/xtrp);
                        atom2_corr{outIdx} = atom2_corr{outIdx} + (atom2_corr_tmp(2)/xtrp);
                    end
                end
            end
        end
    end
    for outIdx = 1:length(outputs)
        if strcmp(func2str(corefunc), 'cFIT')
            atom1 = min(atom1_corr{outIdx});
            atom2 = min(atom2_corr{outIdx});
            corrected_v{outIdx} = atom1-atom2;
            atom1_naive = min(atom1_naive{outIdx});
            atom2_naive = min(atom2_naive{outIdx});
            naive_v{outIdx} = atom1_naive-atom2_naive;
        else
            corrected_v{outIdx} = min(atom1_corr{outIdx},atom2_corr{outIdx});
            naive_v{outIdx} = min(atom1_naive{outIdx}, atom2_naive{outIdx});
        end
    end
elseif strcmp(func2str(corefunc), 'II')
    atom1_corr = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
    atom2_corr = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
    naive_opts = opts;
    naive_opts.bias = 'naive';
    naive_opts.recall = true;
    [~,~,~,atom1_naive, atom2_naive] = feval(corefunc, inputs, outputs, naive_opts);
    nTrials = size(inputs{1},length(size(inputs{1})));
    if opts.parallel == false || ((opts.parallel == true) && (opts.shuff > 0))
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
                    inputs_p = partition(inputs_s, npartition(np), pidx,0);
                    if strcmp(corr, 'qe_shuffSub') || strcmp(corr, 'le_shuffSub')
                        [~, ~, ~, addOut] = shuffle_subtraction(inputs_p, outputs, corefunc, naive_opts);
                        for outIdx = 1:length(outputs)
                            atom1_tmp{outIdx} = addOut{1,outIdx}(end);
                            atom2_tmp{outIdx} = addOut{2,outIdx}(end);
                        end
                    else
                        [~,~,~,atom1_tmp, atom2_tmp] = feval(corefunc, inputs_p,outputs, naive_opts);
                    end
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
                if strcmp(opts.bias,'qe') ||strcmp(opts.bias,'qe_shuffSub')
                     for t = 1:nTimepoints
                        y = [atom1_naive{outIdx}(t), atom1_2{outIdx}(t), atom1_4{outIdx}(t)];
                        p = polyfit(x_extrap, y, 2); 
                        atom1_corr{outIdx}(1,t)  =  atom1_corr{outIdx}(t) + p(3)/xtrp;
                        
                        y = [atom2_naive{outIdx}(t), atom2_2{outIdx}(t), atom2_4{outIdx}(t)];
                        p = polyfit(x_extrap, y, 2);
                        atom2_corr{outIdx}(1,t)  =  atom2_corr{outIdx}(t) + p(3)/xtrp;
                     end
                elseif strcmp(opts.bias,'le')||strcmp(opts.bias,'le_shuffSub')
                    for t = 1:nTimepoints
                        y = [atom1_naive{outIdx}(t), atom1_2{outIdx}(t)];
                        p = polyfit(x_extrap, y, 1); 
                        atom1_corr{outIdx}(1,t)  =  atom1_corr{outIdx}(t) + p(2)/xtrp;
                        
                        y = [atom2_naive{outIdx}(t), atom2_2{outIdx}(t)];
                        p = polyfit(x_extrap, y, 1);
                        atom2_corr{outIdx}(1,t)  =  atom2_corr{outIdx}(t) + p(2)/xtrp;
                    end
                end
            end
        end
    end
    for outIdx = 1:length(outputs)
        if strcmp(func2str(corefunc), 'cFIT')
            atom1 = min(atom1_corr{outIdx});
            atom2 = min(atom2_corr{outIdx});
            corrected_v{outIdx} = atom1-atom2;
            atom1_naive = min(atom1_naive{outIdx});
            atom2_naive = min(atom2_naive{outIdx});
            naive_v{outIdx} = atom1_naive-atom2_naive;
        else
            corrected_v{outIdx} = min(atom1_corr{outIdx},atom2_corr{outIdx});
            naive_v{outIdx} = min(atom1_naive{outIdx}, atom2_naive{outIdx});
        end

    end
elseif strcmp(func2str(corefunc), 'H')
    naive_opts = opts;
    naive_opts.bias = 'naive';
    naive_opts.computeNulldist = 'false';
    naive_v = feval(corefunc, inputs, outputs, naive_opts);
    nTrials = size(inputs{1},length(size(inputs{1})));
    if opts.parallel == false || ((opts.parallel == true) && (opts.shuff > 0))
        for i=1:opts.xtrp
            randidx = randperm(nTrials, nTrials);
            npartition = [1 2 4];
            part2 = repmat({zeros(1,nTimepoints)}, 1, length(outputs));
            part4 = repmat({zeros(1,nTimepoints)}, 1, length(outputs));

            inputs_s = inputs;
            idx = repmat({':'}, 1, ndims(inputs{1}));
            idx{end} = randidx;
            for var = 1:length(inputs)
                inputs_s{var} = inputs{var}(idx{:});
            end
            for np=2:length(npartition)
                for pidx = 1:npartition(np)
                    inputs_p = partition(inputs_s, npartition(np), pidx, 0);                   
                    if strcmp(corr, 'qe_shuffSub') || strcmp(corr, 'le_shuffSub')
                       shuffle_subtraction(inputs_p, outputs, corefunc, naive_opts);
                    else 
                        value_tmp = feval(corefunc, inputs_p, outputs, naive_opts);
                    end 
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
                        p = polyfit(x_extrap, [naive_v{outIdx}(1,t), part2{outIdx}(1,t), part4{outIdx}(1,t)], 2);
                        corrected_v{outIdx}(1,t) = corrected_v{outIdx}(1,t) +(p(3)/xtrp);
                    elseif strcmp(opts.bias,'le')||strcmp(opts.bias,'le_shuffSub')
                        p = polyfit(x_extrap(1:2), [naive_v{outIdx}(1,t), part2{outIdx}(1,t)], 1);
                        corrected_v{outIdx}(1,t) = corrected_v{outIdx}(1,t) +(p(2)/xtrp);
                    end
                end 
            end
        end
    end
end
end






