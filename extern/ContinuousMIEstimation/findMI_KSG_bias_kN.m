function [means, errorBars] = findMI_KSG_bias_kN(X,Y,listOfKs, listSplitSizes, do_plot, do_plot_stddev, do_plot_subsampling, pathSave)

% function [means, errorBars] = findMI_KSG_bias_kN(X,Y,listOfKs, listSplitSizes, do_plot, do_plot_stddev, do_plot_subsampling)
%
% uses KSG k-nearest neighbors mutual information estimator at various values
% of k, allowing to select the best k
% 
% X,Y - continuous variables - should each be given as matrices, so that
% each row corresponds to a single data point, and each column to a
% different dimension. These *must* be the same length as each other, but
% do not need to have the same dimensionality.
% 
% listOfKs - vector of integers, the list of different values of k in the 
% k-nearest neighbors KSG estimator to try. 
% 
% listSplitSizes - list of paritioning of the data to be used for errror bars
% calculations; see detailed description in findMI_KSG_subsampling.m 
% 
% means - values of estimated I(X;Y), in bits, at each k value in listOfKs.
% Vector of the same dimension as listOfKs.
% 
% errorBars - errors for the estimated values of I(X;Y) at each k, in bits.
% Vector of the same dimension as listOfKs.
% 
% do_plot - binary scalar. If 1, plots of dependence of MI estimated values 
% on k will automatically be generated.
% 
% do_plot_stddev, do_plot_subsampling -- binary scalars, are passed directly 
% to findMI_KSG_stddev and findMI_KSG_subsampling and control generation of
% figures in these function; see documentation for these function
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

means = zeros(length(listOfKs),1); %the estimated value of the mutual 
                    %information at each k, in bits
errorBars = zeros(length(listOfKs),1); %the errorbars at each k, in bits
N = length(X);      %We do not check that X and Y are of the same size; 
        % this is done by MIxnyn

for i = 1:length(listOfKs) %for each k value, find the values of the information and the error estimates
    [MIs] = findMI_KSG_subsampling(X,Y, listOfKs(i), listSplitSizes, do_plot_subsampling,pathSave);
    means(i) = mean(MIs{1,2});
    errorBars(i) = findMI_KSG_stddev(MIs, N, do_plot_stddev); 
end


%plot the obtained values to see the relationship between k and N on the 
%one hand and the estimates of the information on the other 

if do_plot == 1
    figure
    errorbar(listOfKs, means, errorBars)
    xlabel('k')
    ylabel('mutual information, bits')
    title('Mutual Information different k values')
end
