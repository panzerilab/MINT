function [ L, num, values ] = createClusters_v2(infQuant, clusterThreshold)
% createClusters creates clusters based of sample above the information
% cluster-forming threshold
%   
% infQuant should be a rectangular matrix of values of the information
% quantity

% take only non-zero values beacuse sometimes we don't compute the whole
% matrix and our quantities shouldn't be exactly zero even if they are very
% low
allValues = infQuant(:);

% just a corner case when all values are smaller than the clusterThreshold
if max(allValues) < clusterThreshold
    L = zeros(size(infQuant));
    num = 0;
    values = [];
else
    % create a bit mask of all those that are above the threshold
    minValid = clusterThreshold;
    valid = infQuant >= minValid;
    
    % create clusters
    [L, num] = bwlabeln(valid,4); % you can set the second parameter to 4 if you want clusters to be connected only with the 4 closest neighbours

    % sum the t-values for them
    values = zeros(num, 1);
    for i = 1:num
        values(i) = sum(sum(infQuant(L == i)));
    end
end

