function [ bitMask ] = clusterStat_v2( infQuant, infQuantSh, significanceThreshold, clusterPercentilThreshold)
%%
% clusterStat computes cluster based statistic of the given infQuant
%   
% infQuant should be a rectangular matrix of values of the information
% quantity
%
% infQuantSh should be a 3D matrix having the first two dimensions equal to
% the size of infQuant and the last one equal to number of shuffling of
% stimulus that was used
% 
% it returns a bitmask of the same size as infQuant having 1s for
% significant values and 0s for not significant values.
%
% It is an implementation of the procedure described here:
% http://www.sciencedirect.com/science/article/pii/S0165027007001707?via%3Dihub
%
% Example:
% mask = clusterStat(infQuant, allrDFIQeSh, 0.05, 0.975);
% computes cluster statistic based on clusters composed of values higher
% than 97.5th percentile. We consider as significant only those that are
% lower than 5th in a rank test.
%

threshold = significanceThreshold;

allVals = infQuantSh(:);
clusterThresh = prctile(allVals,clusterPercentilThreshold*100);

% create clusters of interest from the given 2d matrix
[L, num, clusterVals] = createClusters_v2(infQuant, clusterThresh);

% count, how many times there is a cluster in the shuffled values
% with higher value than the previously computed clusters of interest
higherClusterValCount = zeros(size(clusterVals));
for s = 1:size(infQuantSh, 3)
    [~, ~, surrogateClusterVals] = createClusters_v2(squeeze(infQuantSh(:, :, s)), clusterThresh);
    if isempty(surrogateClusterVals)
        continue
    end
    higherClusterValCount(clusterVals < max(surrogateClusterVals)) = higherClusterValCount(clusterVals < max(surrogateClusterVals)) + 1;
    %maxclus(s) = max(surrogateClusterVals);
end

% put all the significant clusters into one bit mask
bitMask = zeros(size(infQuant));
for c = 1:num
    if (higherClusterValCount(c) / size(infQuantSh, 3) < threshold)
        bitMask(L == c) = 1;
    end
end

bitMask = logical(bitMask);

