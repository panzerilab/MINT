function [p_src, p_crs, n_S, n_R, n_C, n_singleton_dims] = build_p(S, R, C)
%%% *function [p_src, p_crs, n_S, n_R, n_C, n_singleton_dims] = build_p(S, R, C)*
%%%
%%% ### Description
%%% The function estimates the probability distribution p(s,r,c) from the 3D histogram of the input (S, R, C) occurrences.
%%%
%%% ### Inputs:
%%% - *S*: must be a one dimensional array of *n_trials X 1* elements representing the discrete value of the stimulus presented in each trial.
%%% - *R*: must be a one dimensional array of *n_trials X 1* elements representing the response at each trial, here we suppose the response has been binned already.
%%% - *C*: must be a one dimensional array of *n_trials X 1* elements representing the discrete value of the choice made by the subject in each trial.
%%%
%%% ### Outputs:
%%% - *p_src*: joint probability p(s,r,c).
%%% - *p_crs*: joint probability p(c,r,s).
%%% - *n_S*: number of stimuli.
%%% - *n_R*: number of responses.
%%% - *n_C*: number of choices.
%%% - *n_singleton_dims*: number of singleton dimensions.


% Initialize variables and compute the number of unique values for S, R, and C
% Create flags to store data for each dimension for performance

if nargin < 3
    msg = 'not enough input arguments.';
    error('buildp:notEnoughInput', msg);
end

variables = {{S, 'S'}, {R, 'R'}, {C, 'C'}};
nanVariableName = '';
for i = 1:length(variables)
    var = variables{i}{1};
    name = variables{i}{2};
    if any(isnan(var), 'all')
        nanVariableName = name;
        break;
    end
end

if ~isempty(nanVariableName)
    msg = sprintf('%s contains NaNs. Aborting.', nanVariableName);
    error('buildp:NaNInput', msg);
end

R = reshape(R,length(R),1);
R_discrete_values=unique(R);
S_values=unique(S);
C_values=unique(C);

N_trials=length(R(:,1));
n_R=numel(R_discrete_values);
n_S=numel(S_values);
n_C=numel(C_values);

% Preallocate arrays for storing data for each dimension
R_stored=false(N_trials,n_R);
S_stored=false(N_trials,n_S);

% Initialize the joint probability distribution p(s, r, c)
p_src=zeros(n_S,n_R,n_C);
n_singleton_dims = 3-length(size(p_src))+sum(size(p_src)==1);

% Iterate through choice values (c), response values (r), and stimulus values (s) to compute the joint probabilities.
for cc=1:n_C
    for rr=1:n_R
        if cc == 1
            % store for performance
            R_stored(:,rr)=(R==R_discrete_values(rr));
        end
        for ss=1:n_S
            if rr == 1
                % store for performance
                S_stored(:,ss)=(S==S_values(ss));
            end
            % Calculate the joint probability p(s, r, c)
            p_src(ss,rr,cc)=sum(R_stored(:,rr).*S_stored(:,ss).*(C==C_values(cc)),'all');
        end
    end
end
% Normalize the joint probability distribution
p_src=p_src/sum(p_src(:));
% Compute p(c, r, s) by permuting p(s, r, c)
p_crs=permute(p_src,[3 2 1]);

end