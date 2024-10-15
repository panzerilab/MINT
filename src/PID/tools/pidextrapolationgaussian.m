function [corrected, naive] = pidextrapolationgaussian(X1, X2, Y, opts, btspflag)
%%% *function pid_v = pidimmi(pdf_dirty)*
%%%
%%% ### Description
%%% Compute Partial Information Decomposition (PID) values from a given probability distribution (pdf_dirty).
%%%
%%% Partial Information Decomposition separates the information that multiple sources provide about a target into unique, redundant, and synergistic components. 
%%% This function takes a 3D probability distribution (joint pdf of three variables) and calculates the PID values.
%%%
%%% ### Inputs:
%%% - pdf_dirty: A 3-dimensional array representing the joint probability density function of three variables, typically denoted as X, Y, and Z. 
%%%
%%% ### Outputs:
%%% - pid_v: A vector containing four values representing the components of PID:
%%% - Redundant information (common information both sources provide about the target)
%%% - Unique information from the first variable (information only the first source provides about the target)
%%% - Unique information from the second variable (information only the second source provides about the target)
%%% - Synergistic information (information that is only available when both sources are considered together)


if nargin < 3
    msg = 'not enough input arguments.';
    error('extrapolation:notEnoughInput', msg);
end

% sanity_check(inputs);

infopts = [];
infopts.bias = 'naive';
infopts.method = 'gs';
infopts.bin_methodX = 'none'; %opts.bin_methodX; 
infopts.bin_methodY = 'none'; %opts.bin_methodY;
% infopts.n_binsX     = opts.n_binsX; 
% infopts.n_binsY     = opts.n_binsY;
infopts.verbose = false;
infopts.btsp= 0;

if nargin < 3
    msg = 'not enough input arguments.';
    error('extrapolation:notEnoughInput', msg);
end

% sanity_check(inputs);

xtrp = opts.xtrp;
pdfopts.method='none';
ntrials = size(X1,2);

% inputs =[inputs(2,:); inputs(3,:); inputs(1,:)];

I12 = cell2mat(information([X1;X2], Y, infopts, {'I'})); %mutualInformationXYZ(p);
I1  = cell2mat(information(X1, Y, infopts, {'I'}));
I2  = cell2mat(information(X2, Y, infopts, {'I'}));

red = min([I1, I2]);
u1 = I1 - red;
u2 = I2 - red;
syn = I12 - u1 - u2 - red;
naive = [red u1 u2 syn];
corrected = zeros([xtrp, size(naive)]);
if opts.parallel == 0 || ((opts.parallel == 1) && (btspflag ==1))
    for i=1:xtrp
        ri = randperm(ntrials, ntrials);
        npartition = [1 2 4];

        part2 = 0;
        part4 = 0;

        X1_s = X1(:, ri);
        X2_s = X2(:, ri);
        Y_s = Y(:, ri);
        % inputs_s = inputs(:, ri);
        % prob_matrix = pdf(inputs_s, pdfopts);
        % prob_matrix = prob_matrix(2:end, 2:end, 2:end);
        % naive = feval(core_function, prob_matrix, opts);

        for np=2:length(npartition)
            for pidx = 1:npartition(np)

                X1_p = partition_X(X1_s, ntrials, npartition(np), pidx);%partition(inputs_s(1,:), inputs_s(2,:), inputs_s(3,:), npartition(np), pidx,0);%
                X2_p = partition_X(X2_s, ntrials, npartition(np), pidx);
                Y_p = partition_X(Y_s, ntrials, npartition(np), pidx);
                I12 = cell2mat(information([X1_p; X2_p], Y_p, infopts, {'I'})); %mutualInformationXYZ(p);
                I1  = cell2mat(information(X1_p, Y_p, infopts, {'I'}));
                I2  = cell2mat(information(X2_p, Y_p, infopts, {'I'}));
                
                red = min([I1, I2]);
                u1 = I1 - red;
                u2 = I2 - red;
                syn = I12 - u1 - u2 - red;
                
                pid_v = [red u1 u2 syn];

                if npartition(np)==2
                    part2 = part2 + pid_v/2;
                elseif npartition(np)==4
                    part4 = part4 + pid_v/4;
                end
            end
        end
        x_extrapolation = npartition./ntrials;
        
        if length(naive) == 1
            if strcmp(opts.bias,'qe')
                p = polyfit( x_extrapolation, [naive, part2, part4], 2);
                corrected(i) =p(3);
            elseif strcmp(opts.bias,'le')
                p = polyfit( x_extrapolation(1:2), [naive, part2], 1);
                corrected(i) =  p(2);
            end
        else
            for p_fit = 1:length(naive)
                if strcmp(opts.bias,'qe')
                    p = polyfit(x_extrapolation, [naive(p_fit), part2(p_fit), part4(p_fit)], 2);
                    corrected(i,p_fit) = p(3);
                elseif strcmp(opts.bias,'le')
                    p = polyfit(x_extrapolation(1:2), [naive(p_fit), part2(p_fit)], 1);
                    corrected(i,p_fit) = p(2);
                end
            end
        end
    end
else
  
end 
corrected = mean(corrected,1);
corrected = corrected(1,:);

I12 = cell2mat(information([X1;X2], Y, infopts, {'I'})); %mutualInformationXYZ(p);
I1  = cell2mat(information(X1, Y, infopts, {'I'}));
I2  = cell2mat(information(X2, Y, infopts, {'I'}));

red = corrected(1);
u1 = I1 - red;
u2 = I2 - red;
syn = I12 - u1 - u2 - red;
corrected = [red u1 u2 syn];

end