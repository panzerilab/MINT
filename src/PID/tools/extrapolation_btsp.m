function [corrected, naive] = extrapolation_btsp(core_function, inputs, opts, btspflag)

if nargin < 3
    msg = 'not enough input arguments.';
    error('extrapolation:notEnoughInput', msg);
end

sanity_check(inputs);


xtrp = opts.xtrp;
pdfopts.method='none';
ntrials = opts.n_trials;
prob_matrix = pdf(inputs, pdfopts);
naive = feval(core_function, prob_matrix, opts);
corrected = zeros(size(naive));
if opts.parallel == 0 || ((opts.parallel == 1) && (btspflag ==1))
    for i=1:xtrp
        ri = randperm(ntrials, ntrials);
        npartition = [1 2 4];

        part2 = 0;
        part4 = 0;

        inputs_s = inputs(:, ri);
        % prob_matrix = pdf(inputs_s, pdfopts);
        % prob_matrix = prob_matrix(2:end, 2:end, 2:end);
        % naive = feval(core_function, prob_matrix, opts);

        for np=2:length(npartition)
            for pidx = 1:npartition(np)
                ridx = randsample(ntrials, int32(ntrials/npartition(np)));
                sources_p = inputs_s(:,ridx);%partition_X(inputs_s, ntrials, npartition(np), pidx); %partition(inputs_s(1,:), inputs_s(2,:), inputs_s(3,:), npartition(np), pidx,1);%
                [prob_matrix] = pdf(sources_p, pdfopts);
                % prob_matrix = prob_matrix(2:end, 2:end, 2:end);

                if npartition(np)==2
                    part2 = part2 + feval(core_function, prob_matrix, opts)/2;
                elseif npartition(np)==4
                    part4 = part4 + feval(core_function, prob_matrix, opts)/4;
                end
            end
        end
        x_extrapolation = npartition./ntrials;
        
        if length(naive) == 1
            if strcmp(opts.bias,'qe')
                p = polyfit( x_extrapolation, [naive, part2, part4], 2);
                corrected = corrected + p(3)/xtrp;
            elseif strcmp(opts.bias,'le')
                p = polyfit( x_extrapolation(1:2), [naive, part2], 1);
                corrected = corrected + p(2)/xtrp;
            end
        else
            for p_fit = 1:length(naive)
                if strcmp(opts.bias,'qe')
                    p = polyfit(x_extrapolation, [naive(p_fit), part2(p_fit), part4(p_fit)], 2);
                    corrected(p_fit) = corrected(p_fit) + p(3)/xtrp;
                elseif strcmp(opts.bias,'le')
                    p = polyfit(x_extrapolation(1:2), [naive(p_fit), part2(p_fit)], 1);
                    corrected(p_fit) = corrected(p_fit) + p(2)/xtrp;
                end
            end
        end
    end
else
    parforArg = Inf;
    parfor (draw = 1:xtrp, parforArg)
        ri = randperm(ntrials, ntrials);
        npartition = [1 2 4];

        part2 = 0;
        part4 = 0;

        inputs_s = inputs(:, ri);
        prob_matrix = pdf(inputs_s, pdfopts);
        % prob_matrix = prob_matrix(2:end, 2:end, 2:end);
        naive = feval(core_function, prob_matrix, opts);

        for np=2:length(npartition)
            for pidx = 1:npartition(np)
                sources_p = partition_X(inputs, ntrials, npartition(np), pidx);
                [prob_matrix] = pdf(sources_p, pdfopts);
                % prob_matrix = prob_matrix(2:end, 2:end, 2:end);
                if npartition(np)==2
                    part2 = part2 + feval(core_function, prob_matrix, opts)/2;
                elseif npartition(np)==4
                    part4 = part4 + feval(core_function, prob_matrix, opts)/4;
                end
            end
        end
        x_extrapolation = npartition./ntrials;
        if length(naive) == 1
            if strcmp(opts.bias,'qe')
                p = polyfit( x_extrapolation, [naive, part2, part4], 2);
                corrected = corrected + p(3)/xtrp;
            elseif strcmp(opts.bias,'le')
                p = polyfit( x_extrapolation, [naive, part2, part4], 1);
                corrected = corrected + p(2)/xtrp;
            end
        else
            for p_fit = 1:length(naive)
                if strcmp(opts.bias,'qe')
                    p = polyfit(x_extrapolation, [naive(p_fit), part2(p_fit), part4(p_fit)], 2);
                    corrected(p_fit) = naive(p_fit) + p(3)/xtrp;
                elseif strcmp(opts.bias,'le')
                    p = polyfit(x_extrapolation, [naive(p_fit), part2(p_fit), part4(p_fit)], 1);
                    corrected(p_fit) = naive(p_fit) + p(2)/xtrp;
                end
            end
        end
    end
end 