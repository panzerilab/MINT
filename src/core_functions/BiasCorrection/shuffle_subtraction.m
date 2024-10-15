function [corrected_v, naive_v, shuff_all, addOut] = shuffle_subtraction(inputs, outputs, corefunc, varargin)

opts = varargin{1};
if nargin > 5
    inputs_nD = varargin{2};
end
shuff_all = cell(1, length(outputs));
corrected_v = repmat({0}, 1, length(outputs));
naive_v = repmat({0}, 1, length(outputs));
addOut = 0;
naive_opts = opts;
naive_opts.bias = 'naive';

if strcmp(func2str(corefunc), 'PID')
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
        I1_shuff_all = zeros(1, opts.shuff);
        I2_shuff_all = zeros(1, opts.shuff);
        I12_shuff_all = zeros(1, opts.shuff);
    end
    PID_corrected = repmat({0}, 1, length(outputs));
    for sIdx = 1:opts.shuff
        inputs_sh = inputs;
        for var = 1:length(inputs)
            inputs_sh{var} = (shuffle_core(0, inputs_sh{var}', 0, [1, 0]))';
        end
        shuff_v  = feval(corefunc, inputs_sh, outputs, naive_opts);
        for outIdx = 1:length(outputs)
            shuff_all{outIdx} = [shuff_all{outIdx}, shuff_v{outIdx}];
        end
        if opts.pid_constrained
            I1_shuff_all(sIdx) = cell2mat(MI({inputs_sh{1}, inputs_sh{end}}, {'I(A;B)'}, naive_opts));
            I2_shuff_all(sIdx) = cell2mat(MI({inputs_sh{2}, inputs_sh{end}}, {'I(A;B)'}, naive_opts));
            I12_shuff_all(sIdx) = cell2mat(MI({cat(1, inputs_sh{1}, inputs_sh{2}), inputs_sh{end}}, {'I(A;B)'}, naive_opts));
        end
    end

    if opts.pid_constrained
        I1_corrected = I1_naive-mean(I1_shuff_all);
        I2_corrected = I2_naive-mean(I2_shuff_all);
        I12_corrected = I12_naive-mean(I12_shuff_all);
    end
    for outIdx = 1:length(outputs)
        PID_corrected{outIdx} = PID_naive{outIdx} - mean(shuff_all{outIdx});
        naive_v{outIdx} = PID_naive{outIdx};
    end
    if opts.pid_constrained
        if ~isfield(opts, 'chosen_atom')
            opts.chosen_atom = 'Red';
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
                syn_naive  = I12_naive-I1_naive-I2_naive+red;
                unq1_naive  = I1_naive-red;
                unq2_naive  = I2_naive-red;
            case 'Unq1'
                unq1 = PID_corrected{pos};
                red = I1_corrected-unq1;
                syn =  I12_corrected-I2_corrected-unq1;
                unq2 = I2_corrected-I1_corrected+unq1;
                %naive
                unq1_naive = PID_naive{pos};
                red_naive = I1_naive-unq1;
                syn_naive =  I12_naive-I2_naive-unq1;
                unq2_naive = I2_naive-I1_naive+unq1;
            case 'Unq2'
                unq2 = PID_corrected{pos};
                red = I2_corrected-unq2;
                syn =  I12_corrected-I1_corrected-unq2;
                unq1 = I1_corrected-I2_corrected+unq2;
                %naive
                unq2_naive = PID_naive{pos};
                red_naive = I2_naive-unq2;
                syn_naive =  I12_naive-I1_naive-unq2;
                unq1_naive = I1_naive-I2_naive+unq2;
            case 'Syn'
                syn = PID_corrected{pos};
                red = I1_corrected+I2_corrected-I12_corrected+syn;
                unq1 = I12_corrected-I2_corrected-syn;
                unq2 = I12_corrected-I1_corrected-syn;
                %naive
                syn_naive  = PID_naive{pos};
                red_naive  = I1_naive+I2_naive-I12_naive+syn;
                unq1_naive  = I12_naive-I2_naive-syn;
                unq2_naive  = I12_naive-I1_naive-syn;
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
elseif strcmp(func2str(corefunc), 'FIT') || strcmp(func2str(corefunc), 'cFIT') || strcmp(func2str(corefunc), 'II')
    atom1_corrected = repmat({0}, 1, length(outputs));
    atom2_corrected = repmat({0}, 1, length(outputs));
    [FIT_naive,~,~,atom1_naive, atom2_naive] = feval(corefunc, inputs, outputs, naive_opts);
    addOut = cell(2,length(outputs));
    for sIdx = 1:opts.shuff
        inputs_sh = inputs;
        for var = 1:length(inputs)
            inputs_sh{var} = (shuffle_core(0, inputs_sh{var}', 0, [1, 0]))';
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
        corrected_v{outIdx} = min(atom1_corrected{outIdx},atom2_corrected{outIdx});
        naive_v{outIdx} = FIT_naive{outIdx};
    end



elseif strcmp(func2str(corefunc), 'II')

elseif strcmp(func2str(corefunc), 'MI')
    [corrected_v, naive_v, shuff_all] =  MI(inputs, outputs, opts);
elseif strcmp(func2str(corefunc), 'H')
    possibleOutputs = { 'H(A|B)', 'Hind(A|B)', 'Hsh(A|B)', 'Chi(A)','Hind(A)'};
    naive_v = feval(corefunc, inputs, outputs, naive_opts);
    for sIdx = 1:opts.shuff
        inputs_sh = inputs;
        for var = 1:length(inputs)
            inputs_sh{var} = (shuffle_core(0, inputs_sh{var}', 0, [1, 0]))';
        end
        shuff_v  = feval(corefunc, inputs_sh, outputs, naive_opts);
        for outIdx = 1:length(outputs)
            shuff_all{outIdx} = [shuff_all{outIdx}, shuff_v{outIdx}];
        end
    end
    for outIdx = 1:length(outputs)
        if ismember(outputs{outIdx}, possibleOutputs)
            corrected_v{outIdx} = {naive_v{outIdx}, shuff_all{outIdx}};
        else
            corrected_v{outIdx} = naive_v{outIdx};
        end
    end
end

