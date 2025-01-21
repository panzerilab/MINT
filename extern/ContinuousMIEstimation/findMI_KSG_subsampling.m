function [MIs] = findMI_KSG_subsampling(X,Y, kvalue, listSplitSizes, do_plot,pathSave)

% function [MIs] = findMI_KSG_subsampling(X,Y, kvalue, listSplitSizes, do_plot)
% 
% uses kraskov-grassberger to estimate I(X;Y)
%
% X,Y - continuous variables - should each be given as matrices, so that
% each row corresponds to a single data point, and each column to a
% different dimension. These *must* be the same length as each other, but
% do not need to have the same dimensionality.
%
% kvalue - scalar, positive integer, number of nearest neighbors to count 
% when using the KSG estimator.
% 
% listSplitSizes - 1-d vector of required partitionings of the data when 
% looking at the dependence of the mutual information on the data set size. 
% For example, listSplitSizes = [1,2,3] means that we want to estimate the
% information from the whole data, then from two non-overlapping subsets of 
% the data, each of size N/2, and then from three non-obverlapping subsets 
% each of size N/3. The first element of listSplitSizes should always be 1.
% 
% do_plot -- a binary scalar variable. If 1, plots showing dependence of the 
% mutual information on the sample size will automatically be generated.
% 
% MIs - a cell array, with length(listSplitSizes) rows. Each row contains 
% mutual information data from a different subsample (see description of 
% listSplitSizes). There are two columns - the first entry contain the
% SplitSize, the number of partitions of the data used for that row; the 
% second entry contain the estimated mutual informations in bits for each of 
% SplitSize partitions of the data.
%
%-----------------------------------------------------------------------------------------
% Copyright 2018 Caroline Holmes, Ilya Nemenman
%-----------------------------------------------------------------------------------------
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
% You should receive a copy of the GNU General Public License
% along with this program.  See also <http://www.gnu.org/licenses/>.
%-----------------------------------------------------------------------------------------
% Please reference
% 
% Holmes, C.M. & Nemenman, I.  Estimation of mutual information for
% real-valued data with error bars and controlled bias. 
% Submitted, 2018.
%
% A. Kraskov, H. Stogbauer, and P. Grassberger,
% Estimating mutual information.
% Phys. Rev. E 69 (6) 066138, 2004
%
% in your published research.

if size(X,2) > size(X,1)
    X = X';
end
if size(Y,2) > size(Y,1)
    Y = Y';
end


n = length(listSplitSizes); 
MIs = cell(n,2);

N = length(X);

means = zeros(n,1);
errorBars = zeros(n,1);

for i = 1:n
    MIs{i,1} = listSplitSizes(i);
end

for i = 1:n
    
    %randomly permute the data to assign it at random to one of partitions
    a = randperm(length(X)); 
    nSplits = listSplitSizes(i); %number of partitions the data is split into
    l = round(linspace(0,length(X),nSplits + 1));
    
    MI_T = zeros(nSplits,1); %a temporary variable to hold the different 
                %mutual information values for each data partition
                
    for j = 1: nSplits %move through each of the subsets of equal size
        xT = X(a(l(j) + 1:l(j+1)),:); % the 'X' values that belong to this subset
        yT = Y(a(l(j) + 1:l(j+1)),:); % the 'Y' values that belong to this subset
        MI_T(j) = MIxnyn_matlab(xT,yT,kvalue,pathSave); %call the KSG MI estimator function
    end
    MIs{i,2} = MI_T./log(2); %we need this log(2) to convert from nats to bits.
    means(i) = mean(MIs{i,2}); %mean value of MIs of partitions of this size
    errorBars(i) = std(MIs{i,2}); %std.dev. of MIs of partitions of this size
end

if do_plot == 1
    figure
    errorbar(listSplitSizes./N, means, errorBars)
    xlabel('1/N')
    ylabel('Mutual Information, bits')
    title(horzcat('Dependence of Mutual Information on the data set size for k = ', num2str(kvalue)))
end