function data_out = btstrp_shuff(data_in,n_boot)
% Data in should have dimensions 1 = nSubj/nReps and last = nShuff

nReps = size(data_in,1);
nDim = numel(size(data_in));

size_surrogate = size(data_in);% size of single surrogate dataset
nShuff = size_surrogate(end);
size_surrogate(end) = [];

size_out = size(data_in);
size_out(1) = [];

data_out = nan(size_out);

for bIdx = 1:n_boot
    randIdxs = randi(nShuff,1,nReps); % idxs of shufflings to take from each repetition
    tmpSurrogate = nan(size_surrogate);
    if nDim == 3
        for repIdx = 1:nReps
            tmpSurrogate(repIdx,:) = data_in(repIdx,:,randIdxs(repIdx));
        end
        data_out(:,bIdx) = squeeze(mean(tmpSurrogate,1));
    elseif nDim == 4
        for repIdx = 1:nReps
            tmpSurrogate(repIdx,:,:) = data_in(repIdx,:,:,randIdxs(repIdx));
        end
        data_out(:,:,bIdx) = squeeze(mean(tmpSurrogate,1));
    elseif nDim == 5
        for repIdx = 1:nReps
            tmpSurrogate(repIdx,:,:,:) = data_in(repIdx,:,:,:,randIdxs(repIdx));
        end
        data_out(:,:,:,bIdx) = squeeze(mean(tmpSurrogate,1));
    elseif nDim == 6
        for repIdx = 1:nReps
            tmpSurrogate(repIdx,:,:,:,:) = data_in(repIdx,:,:,:,:,randIdxs(repIdx));
        end
        data_out(:,:,:,:,bIdx) = squeeze(mean(tmpSurrogate,1));
    end
end

end