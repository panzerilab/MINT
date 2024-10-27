function partitioned_data = partition(inputs, nparts, partidx,conservePs)
% partition - Partition input data into specified subsets based on trial indices or unique stimuli.
%
% This function partitions the provided input data into specified subsets. It can operate in two modes: 
% either dividing trials evenly across partitions or grouping data based on unique stimulus values 
% while maintaining the proportion of samples within each partition. The choice of partitioning method is controlled 
% by the `conservePs` parameter, which indicates whether to conserve the proportions of stimuli.
%
% Inputs:
%   - inputs: A cell array containing data arrays to be partitioned. The size and dimensions of each 
%             array can vary, but typically include data for different variables across trials.
%   - nparts: A scalar value specifying the number of partitions to create from the input data.
%   - partidx: A scalar value representing the index of the current partition to extract.
%   - conservePs: A boolean flag (1 or 0) indicating whether to conserve proportions of unique stimuli 
%                 when partitioning the data. If set to 1, the function will group data based on 
%                 unique stimulus values; if 0, it will split the data based on trial indices.
%
% Outputs:
%   - partitioned_data: A cell array containing the partitioned datasets. The size of this array matches 
%                       the number of input variables, with each cell containing the corresponding 
%                       partitioned data for that variable.
%
% Method:
% The function checks the dimensions of the input data to determine the appropriate partitioning strategy. 
% If `conservePs` is set to 1, the function will identify unique stimulus values and partition the data 
% accordingly, ensuring that each partition maintains a proportionate representation of the stimuli. 
% If `conservePs` is set to 0, the function will simply split the input data based on the specified number 
% of partitions and the trial indices.
%
% Example:
% To partition data from two input variables across 4 partitions while conserving stimulus proportions:
%
%     X1 = rand(100, 10);   % Example input variable 1
%     X2 = rand(100, 10);   % Example input variable 2
%     Y = randi([1, 3], 100, 1); % Example stimulus variable
%
%     partitioned_data = partition({X1, X2, Y}, 4, 1, 1);
%
% This call will return the first partition of the input data arrays, ensuring that the partition maintains 
% the proportions of unique stimuli present in the original dataset.

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

shape_data = size(inputs{1});

nVars   = length(inputs);
if length(shape_data)==2
    nTrials = shape_data(end);
elseif length(shape_data)==3
    nTimpoints = shape_data(2);
    nTrials = shape_data(end);
end

part= [];
if conservePs==1 && length(shape_data)==2
    X1 = inputs{1};
    X2 = inputs{2};
    Y = inputs{end};
    [uniqStim,~,~] = unique(Y);
    for stim=1:length(uniqStim)
        mask= (Y==uniqStim(stim));
        Ytotstim =  Y(mask);
        X1totstim = X1(mask);
        X2totstim = X2(mask);
        if partidx ==1
            bin_edges = 1:round(length(Ytotstim)/nparts);
        elseif partidx ==nparts
            bin_edges = (partidx-1)*round(length(Ytotstim)/nparts)+1:length(Ytotstim);
        else
            bin_edges = (partidx-1)*round(length(Ytotstim)/nparts)+1:partidx*round(length(Ytotstim)/nparts);
        end
        if length(shape_data)==2
            part = [part [X1totstim(1,bin_edges); X2totstim(1,bin_edges); Ytotstim(1,bin_edges)]];
        elseif length(shape_data)==3
            part = [part [X1totstim(1,:,bin_edges); X2totstim(1,:,bin_edges); Ytotstim(1,:,bin_edges)]];
        end
       
    end

    partitioned_data = cell(nVars,1);
    if length(shape_data) == 2
            partitioned_data{1} = part(1,:);
            partitioned_data{2} = part(2,:);
            partitioned_data{3} = part(3,:);
        elseif length(shape_data) == 3
            partitioned_data{1} = part(1,:,:);
            partitioned_data{2} = part(2,:,:);
            partitioned_data{3} = part(3,:,:);
    end
    
else
    if partidx ==1
        bin_edges = 1:round(nTrials/nparts);
    elseif partidx ==nparts
        bin_edges = (partidx-1)*round(nTrials/nparts)+1:nTrials;
    else
        bin_edges = (partidx-1)*round(nTrials/nparts)+1:partidx*round(nTrials/nparts);
    end


    partitioned_data = cell(1,nVars);
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
end