function [ks, unreparamaterizedEstimates, unreparamaterizedErrorBars, ...
    reparamaterizedEstimates, reparamaterizedErrorBars, mutualInformationAnalytical] = heavyTailExample()
% 
% function [ks, unreparamaterizedEstimates, unreparamaterizedErrorBars, ...
%    reparamaterizedEstimates, reparamaterizedErrorBars, mutualInformationAnalytical] = heavyTailExample()
%
% This is an example of the effects of reparamaterization, using an extreme
% case with very heavy-tailed data. This is similar to Ref. [1],
% Fig. 4. Note that the form of the k and N dependence may depend on 
% on exactly which sample points are generated -- in long-tailed data, a
% few outliers can have a dramatic effect.
%
% Outputs:
%   ks - this is the list of k values used
%
%   unreparamaterizedEstimates - this is a vector of the estimated mutual
%   information at the values of k, for the full data set size.
%
%   unreparamaterizedErrorBars - this is a vector of the estimated error bars
%   for the mutual information estimates stored in unreparamaterizedEstimates
%
%   reparamaterizedEstimates - this is equivalent to
%   unreparamaterizedEstimates, but estimates were done after
%   reparamaterizing the data
%
%   reparamaterizedErrorBars - this is equivalent to
%   unreparamaterizedErrorBars, but for reparamaterized data
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

% Compile C code
mex MIxnyn.C

% Set path for temp files
pathToSave = pwd;

%We want data where we can easily analytically calculate the mutual information,
%so we will take a reparamaterized bivariate gaussian. Here, X is
%log-normal distributed, while Y is even more heavy-tailed
correlation = 0.5;
N = 10000;

x1 = normrnd(0, 1, N, 1);
x2 = normrnd(0, 1, N, 1);
x3 = correlation .* x1 + (1 - correlation^2)^.5 .* x2;
X = [exp(x1)];
Y = [exp(exp(x3))]; 
mutualInformationAnalytical = log2((1/(1-correlation^2)).^.5);





ks = [1 2 3 4 5 20]; %this are the values of k we're using
listSplitSizes = [1:10]; %this is a standard choice


[values, errorBars] = findMI_KSG_bias_kN(X,Y,ks, listSplitSizes, 1, 1, 1, pathToSave); %run MI estimation

figure;
plot(ks, values, errorbars)
xlabel('K')
ylabel('MI')

%%
% alternatively, we could call findMI_KSG_subsampling, which will give us the results from the subsampling individually for each K. 

%for each k, we need a list of values of the mutual information and error
%bars for those values.
means = zeros(length(ks), length(listSplitSizes)); 
stds = zeros(length(ks), length(listSplitSizes));

%instead of calling findMI_KSG_bias_kN.m, we can directly call
%findMI_KSG_subsampling.
for i = 1:length(ks)
    [MIs] = findMI_KSG_subsampling(X,Y, ks(i), listSplitSizes, 0, pathToSave);
    for j = 1:length(listSplitSizes)
        means(i,j) = mean(MIs{j,2});
        stds(i,j) = std(MIs{j,2});
    end
end

figure
hold on
errorbar(listSplitSizes + 0.1,means(1,:),stds(1,:),'b') %the 0.1 offset for the x-axis is just for visualization purposes
errorbar(listSplitSizes,means(2,:),stds(2,:),'g')
errorbar(listSplitSizes - 0.1,means(3,:),stds(3,:),'r') 
xlabel('Number of divisions of the data')
ylabel('Estimated Mutual Information')
legend(horzcat('k = ', num2str(ks(1))), horzcat('k = ', num2str(ks(2))),horzcat('k = ', num2str(ks(3))))
title('N-dependence, unreparameterized data')
plot([0 11], [mutualInformationAnalytical mutualInformationAnalytical], 'k')
xlim([0 11])

unreparamaterizedEstimates = means(:,1);
unreparamaterizedErrorBars = stds(:,1);

%Reparameterize the data to compare the effects of the reparameterization.
[X2] = reparameterize_data(X);
[Y2] = reparameterize_data(Y);



%and now we can generate a similar plot, using the reparameterized data
for i = 1:length(ks)
    [MIs] = findMI_KSG_subsampling(X2,Y2, ks(i), listSplitSizes, 0,pathToSave);
    for j = 1:length(listSplitSizes)
        means(i,j) = mean(MIs{j,2});
        stds(i,j) = std(MIs{j,2});
    end
end

figure
hold on
errorbar(listSplitSizes + 0.1,means(1,:),stds(1,:),'b')
errorbar(listSplitSizes,means(2,:),stds(2,:),'g')
errorbar(listSplitSizes - 0.1,means(3,:),stds(3,:),'r')
xlabel('Number of divisions of the data')
ylabel('Estimated Mutual Information')
legend(horzcat('k = ', num2str(ks(1))), horzcat('k = ', num2str(ks(2))),horzcat('k = ', num2str(ks(3))))
title('N-dependence, reparameterized data')
plot([0 11], [mutualInformationAnalytical mutualInformationAnalytical], 'k')
xlim([0 11])

reparamaterizedEstimates = means(:,1);
reparamaterizedErrorBars = stds(:,1);

end
