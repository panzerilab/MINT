function [listOfKs, mutualInformationEstimates, errorBars, mutualInformationAnalytical] = gaussianExample(correlation, numPoints)

% function [listOfKs, mutualInformationEstimates, errorBars, mutualInformationAnalytical] = gaussianExample(correlation, numPoints)
%    
% This function performs mutual information estimation for correlated Gaussian 
% data. It generates figures similar to Figs. 2, 3 in (Holmes and Nemenman, 2018)
% it illustrated how error bars are calculated, and how mutual information
% estimates depends on k and on N for these data. Since for bivariate
% Gaussians we can calculate the mutual information analytically, this
% example allows us to compare the estimates with true information values. 
% 
% This function will iterates across values of 'k' and plots an equivalant
% to Fig. 2 (estimation of error bars) at each k. It also outputs the 
% dependence of the estimated mutual information on k.
%
% correlation - this is the correlation of the bivariate gaussian
%
% numPoints - the number of points in the generated data
%
% listOfKs - this is just a list of k values used
%
% mutualInformationEstimates - the estimated value of the mutual
% information at each k
%
% errorBars - the estimated error bars on the mutual information estimates
% at each k
%
% mutualInformationAnalytical - the analytical value of the mutual
% information between X and Y, for comparison with the estimates
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


% Compile the C code:
mex MIxnyn.C

% Set path for temporary file:
pathToSave = pwd;


muX = 0;         %X mean of the bivariate gaussian
muY = 0;         %Y mean of the bivariate gaussian
sigmaX = 1;      %X standard deviation   
sigmaY = 1;      %Y standard deviation

%generate a correlated X and Y
x1 = normrnd(0, 1, numPoints, 1);
x2 = normrnd(0, 1, numPoints, 1);
x3 = correlation .* x1 + (1 - correlation^2)^.5 .* x2;
X = muX + x1 * sigmaX;
Y = muY + x3 * sigmaY;


disp(horzcat('correlation = ', num2str(correlation)))


%calculate actual expected mutual information between X and Y
mutualInformationAnalytical = - 1/2* log2(1-correlation^2);
disp(horzcat('Analytical value of the mutual information = ', num2str(mutualInformationAnalytical),' bits'))


%now we will do the estimates. This will generate quite a lot of plots.
listOfKs = [1 2 3 4 5 7 10 15 20]; %values of k in k-nearest neighbors to try
listSplitSizes = [1 2 3 4 5 6 7 8 9 10]; %list of partitionings of data to try
[values, errorBars] = findMI_KSG_bias_kN(X,Y,listOfKs, listSplitSizes, 1, 1, 1, pathToSave); %runs the actual estimates. the last three entries determine whether or not various plots are generated during the estimation process.

disp(horzcat('The values of k used are ', num2str(listOfKs)))
disp(horzcat('Estimated values of mutual information at different k"s are ', num2str(values'), ' and standard deviations are ', num2str(errorBars')))

mutualInformationEstimates = values;

end

 
