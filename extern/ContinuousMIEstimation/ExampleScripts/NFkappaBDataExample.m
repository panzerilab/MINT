function [ks, unreparamaterizedEstimates, unreparamaterizedErrorBars, reparamaterizedEstimates, ...
    reparamaterizedErrorBars] = NFkappaBDataExample()
%function [ks, unreparamaterizedEstimates, unreparamaterizedErrorBars, reparamaterizedEstimates, ...
%    reparamaterizedErrorBars] = NFkappaBDataExample()
%
% Runs through some different analyses of the NF-kB data - comparing k and N
% dependence of mutual information estimates, and showing the effects of
% reparameterization. The data set is described in more detail in Ref. [1].
% Briefly, the data describe the joint activity of two transcription factors 
% NF-kappaB and p-ATF-2 measured in 335 individual wildtype mouse fibroblast
% cells 30 min after exposure to the tumor necrosis factor (TNF) ligand at 
% the concentration of 1.3 ng/mL. The two transcription factors are 
% activated downstream of the same TNF receptor, and hence their activity 
% is correlated. The mutual information between these two sets quantifies 
% this relation. 
%
% Returned values:
%   ks - this is the list of k values used
%
%   unreparamaterizedEstimates - this is a vector of the estimated mutual
%   information at the values of k, for the full data set size.
%
%   unreparamaterizedErrorBars - this is a vector of the estimated error bars
%   for the mutual information estimates stored in unreparamaterizedEstimates
%
%   reparamaterizedEstimates - this is equivalent to unreparamaterizedEstimates,
%   but estimates were done after reparamaterizing the data into marginal
%   standard normals.
%
%   reparamaterizedErrorBars - this is equivalent to 
%   unreparamaterizedErrorBars, but for the reparamaterized data.
%
%-----------------------------------------------------------------------------------------
% Copyright 2018, 2019 Caroline Holmes, Ilya Nemenman
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

load('NFkappaBData.mat')

% Compile the C code:
mex MIxnyn.C

% Set path for temporary file:
pathToSave = pwd;

%First, let's generate some basic plots, showing the dependence of the
%mutual information estimate on N and k
listOfKs = [1 2 3 4 5 7 10 15 20];
listSplitSizes = [1:10];
[values errorbars] = findMI_KSG_bias_kN(X,Y,listOfKs, listSplitSizes, 1, 0, 1, pathToSave);





%Instead, we could make a single plot similar to figure 7 in Holmes &
%Nemenman.
ks = [1:10 15 20]; %to put it all on one plot, we need a shorter list of ks.

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
ylim([0 1])

unreparamaterizedEstimates = means(:,1);
unreparamaterizedErrorBars = stds(:,1);

%Reparameterize the data to compare the effects of the reparameterization.
[X2] = reparameterize_data(X);
[Y2] = reparameterize_data(Y);



%and now we can generate a similar plot, using the reparameterized data
for i = 1:length(ks)
    [MIs] = findMI_KSG_subsampling(X2,Y2, ks(i), listSplitSizes, 0, pathToSave);
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
ylim([0 1])

reparamaterizedEstimates = means(:,1);
reparamaterizedErrorBars = stds(:,1);

end
