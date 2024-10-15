function data_out = btsp_shuff_groupLevel(data, null_samples)
% btsp_shuff_groupLevel performs a bootstrap shuffling procedure at the group level to create a null distribution of data. 
% This function is designed for multidimensional arrays (3D to 6D) representing repeated measures 
% 
% Inputs:
% - data: A multidimensional array containing the data to analyze, structured as (repetitions, dimensions, shuffles) or higher dimensions.
%         first dimension has to be nSubj/nReps and last dimension nShuff
% - null_samples: An integer indicating the number of null samples (shuffles) to generate.
% 
% Output:
% - data_out: A multidimensional array of shuffled data averages, with the shape depending on the input dimensions minus the repetitions.
%
% The function generates a null distribution by reshuffling data across repetitions and calculating the mean for each shuffle, 
% allowing for the assessment of group-level interference. This facilitates statistical testing by comparing observed effects 
% against the variability expected under the null hypothesis.

% This method is developed by Hamed Nili
% (https://scholar.google.co.uk/citations?user=QBgOje0AAAAJ&hl=en)


nReps = size(data,1);
nDim = numel(size(data));

size_surrogate = size(data);
nShuff = size_surrogate(end);
size_surrogate(end) = [];

size_out = size(data);
size_out(1) = [];

data_out = nan(size_out);

for bIdx = 1:null_samples
    randIdxs = randi(nShuff,1,nReps); 
    tmpSurrogate = nan(size_surrogate);
    if nDim == 3
        for repIdx = 1:nReps
            tmpSurrogate(repIdx,:) = data(repIdx,:,randIdxs(repIdx));
        end
        data_out(:,bIdx) = squeeze(mean(tmpSurrogate,1));
    elseif nDim == 4
        for repIdx = 1:nReps
            tmpSurrogate(repIdx,:,:) = data(repIdx,:,:,randIdxs(repIdx));
        end
        data_out(:,:,bIdx) = squeeze(mean(tmpSurrogate,1));
    elseif nDim == 5
        for repIdx = 1:nReps
            tmpSurrogate(repIdx,:,:,:) = data(repIdx,:,:,:,randIdxs(repIdx));
        end
        data_out(:,:,:,bIdx) = squeeze(mean(tmpSurrogate,1));
    elseif nDim == 6
        for repIdx = 1:nReps
            tmpSurrogate(repIdx,:,:,:,:) = data(repIdx,:,:,:,:,randIdxs(repIdx));
        end
        data_out(:,:,:,:,bIdx) = squeeze(mean(tmpSurrogate,1));
    end
end

end

