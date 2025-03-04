function [errorEstimate] = findMI_KSG_stddev(MIs, N, do_plot)

% function [errorEstimate] = findMI_KSG_stddev(MIs, N, do_plot)
% 
% MIs - the output from findMI_KSG_subsampling. Cell array; see detailed 
% description in findMI_KSG_subsampling, values in bits.
% 
% N - scalar, total size of the data set: how many data points were there 
% in the original X and Y variables, for which we are estimating I(X;Y)?
% 
% do_plot -- binary scalar. If 1, plots for estimation of MI error bars will 
% automatically be generated.
% 
% errorEstimate - a scalar, the estimate of the standard deviation of the KSG
% mutual information estimator, in bits, at the full data set size N
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

listSplitSizes = MIs(:,1);  % the list of partition sizes
MIs = MIs(:,2); %now MIs only contains the mutual information data, and not
            %the number of partitions
            
n = length(MIs); % number of different partition sizes used to partition the
            %data

Nsamples = zeros(n, 1); 
for i = 1:n
    Nsamples(i) = floor(N./(listSplitSizes{i,1})); %how many data points were actually in each subset for each partition size
end

Xvals = 1./Nsamples; %because we want to look at dependence on 1/N, not N. 
xActual = Xvals(1);  %'actual' meaning with using the entire available data set
Xvals = Xvals(2:end); %We only want to do linear fits for where subsets of 
            %the data areused (where there is a nonzero variance)

listVariances = zeros(n - 1,1); %n-1 because there is no variance in the first row, which we assume is where the entire data set was used
for i = 2:n
    listVariances(i - 1) = var(MIs{i,1});  
end


%now need to fit line with slope = 1 to the above: we want the line to be
%of the form that log(variance) = log(Xvals) + b, and we only want to find
%b

%we are only using this line as a visual guide , so we will do this without
%worrying about how to properly weigh the points.
b = mean(log10(listVariances) -log10(Xvals));

%create the line we'll plot
Xs = linspace(xActual./2,max(Xvals), 100);
Ys = Xs .* (10^b); %and then we'll plot everything in log space
 

%The following estimated the error for the actual mutual information. We
%won't do this by projecting along the fitted line, but rather using a
%chi^2 distribution as described in the paper, although the two answers
%shouldn't be too different.

k = cell2mat(listSplitSizes(2:end)); %because we can't use the first entry here (listSplitSizes = 1) to help us predict the actual variance. 
variancePredicted = sum((k-1)./k.*listVariances)./sum((k-1)); %The estimated variance of the mutual information at the full N
Sml=variancePredicted*N;
varS = 2*Sml^2/sum((k-1)); %Estimated variance of the variance of the estimate of the mutual information at full N
stdvar = sqrt(varS/N^2); %the error bars on our estimate of the variance


if do_plot == 1   
    figure
    hold on
    plot(Xvals, listVariances, 'x')
    plot(Xs,Ys)

    plot(xActual, variancePredicted, 'ks','MarkerFaceColor','k') %draws the projected point
    plot([xActual xActual], [variancePredicted - stdvar, variancePredicted + stdvar],'k-')

    ax = gca;
    ax.XScale = 'log';
    ax.YScale = 'log';
    xlim([min(Xs) max(Xs)]);
    
    xlabel('1/N')
    ylabel('Variance of MI, in bits')
    
    title('Estimating error bars')
    
    legend('Variances from splitting data', '1/N scaling fit to variances from split data','Estimated variance for full data size')
end

errorEstimate = variancePredicted^.5;
