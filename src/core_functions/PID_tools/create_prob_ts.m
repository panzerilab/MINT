function p_ts = create_prob_ts(p_distr, sources)
%%% function p_ts = create_prob_ts(p_distr, sources)
%%%
%%% ### Description
%%% Function to create a bivariate prob distribution between target and sources. If sources is just one variable, it's the same as marginalizing
%%% the other variables. If, for example, the target is X with size 2 and the sources are two variables Y and Z with number of values 3 and 4, the
%%% resulting matrix is (|X|, |Y|*|Z|)
%%%
%%% ### Inputs:
%%% - p_distr: A multidimensional array representing the original probability distribution.
%%% - sources: A vector indicating the indices of the dimensions in p_distr that should be treated as sources.
%%%
%%% ### Outputs:
%%% - p_ts: A two-dimensional matrix representing the joint probability distribution between
%%%         the target and source variables.

if nargin < 2
    msg = 'not enough input arguments.';
    error('createprobts:notEnoughInput', msg);
end

if  any(isnan(p_distr), 'all')
    msg = 'p_distr contains NaNs. Aborting.';
    error('createprobts:NaNInput', msg);
end


if length(sources)+1~=ndims(p_distr)
    p_ts = squeeze(sum(p_distr,setdiff(1:ndims(p_distr), [sources ndims(p_distr)])));
    new_dims = circshift(1:ndims(p_ts),1);
    p_ts = permute(p_ts, new_dims);
else
    p_ts = p_distr;
end

dim_ranges={};
size_array = size(p_ts);
for i=1:ndims(p_ts)
    dim_ranges{i}=1:size_array(i);
end

if ndims(p_ts)==3
    [source1,source2,target]=ndgrid(dim_ranges{1}, dim_ranges{2}, dim_ranges{3});
    table_ts = array2table([target(:),source1(:),source2(:),p_ts(:)], 'VariableNames', {'target', 'source_1', 'source_2', 'prob'});
    table_ts.sourcesingle = table_ts.source_1 + (table_ts.source_2 - 1)*size(unique(table_ts.source_1),1);
    reshaped_p_ts = zeros(size(unique(table_ts.target),1), size(unique(table_ts.sourcesingle),1));
    reshaped_p_ts(sub2ind(size(reshaped_p_ts), table_ts.target, table_ts.sourcesingle)) = table_ts.prob;
    p_ts = reshaped_p_ts;
end

