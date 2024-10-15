function partitioned_data = partition(inputs, nparts, partidx)

shape_data = size(inputs{1});
nVars   = length(inputs);
nTrials = shape_data(end);
if partidx ==1
    bin_edges = 1:round(nTrials/nparts);
elseif partidx ==nparts
    bin_edges = (partidx-1)*round(nTrials/nparts)+1:nTrials;
else
    bin_edges = (partidx-1)*round(nTrials/nparts)+1:partidx*round(nTrials/nparts);
end

partitioned_data = cell(nVars,1);
if length(shape_data)==2
    for varidx =1:nVars
        inputidx = inputs{varidx};
        partitioned_data{varidx} = inputidx(:,bin_edges);
    end
elseif length(shape_data)==3
    for varidx =1:nVars
        inputidx = inputs{varidx};
        partitioned_data{varidx} = inputidx(:,:,bin_edges);
    end
end

end