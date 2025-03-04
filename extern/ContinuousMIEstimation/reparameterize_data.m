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
% Phys. Rev. E, 2019.
%
% A. Kraskov, H. Stogbauer, and P. Grassberger,
% Estimating mutual information.
% Phys. Rev. E 69 (6) 066138, 2004
%
% in your published research.


function [transformedValues] = reparameterize_data(listOfValues)

%This reparamaterizes our data to a gaussian. listOfValues can represent
%multi-dimensional data, and the reparameterization will be done on each
%dimension independently. If reparamaterization on only some dimensions is
%desired, only those dimensions should be input to this function.

%This function assumes that the listOfValues is a matrix, where the larger
%dimension counts through samples of the distribution, and the smaller
%dimension counts through dimensions of the data.

%This reparamaterization does not change the ordering of the points -
%instead, it just changes the distances along the axis so that you have a
%gaussian distribution. 

% However, in order to remove a potential artifact associated with data
% points with the same value, we reorder the data randomly before
% reparamaterizing, and then return it to its original ordering.
if size(listOfValues,2) > size(listOfValues,1)
    listOfValues = listOfValues';
end



if min(size(listOfValues)) == 1

    %reorder the data:


    ordering = randperm(length(listOfValues));
    newListOfValues = listOfValues(ordering);


    [A,I] = sort(newListOfValues); 

    listPositions = zeros(length(listOfValues),1);

    for i = 1:length(listOfValues)
        listPositions(I(i)) = i;
    end

    %reassign to be in range where can use erfinv: need to be between -1 and 1.

    listPositions = listPositions-mean(listPositions);
    listPositions = 2.*listPositions./(length(listPositions));

    randOrdering_transformedValues = erfinv(listPositions);
    
    transformedValues = zeros(size(listOfValues));
    %return to the original ordering
    for i = 1:length(listOfValues)
        transformedValues(ordering(i)) = randOrdering_transformedValues(i);
    end

else
    % for the case where we have multiple dimensions. Essentially, this just
    % loops through the above behaviour, operating on each dimension
    % independently.
    
    x2 = zeros(size(listOfValues));
    input = listOfValues;
    for iiii = 1:min(size(input)) 
        listOfValues = input(:,iiii);
        ordering = randperm(length(listOfValues));
        newListOfValues = listOfValues(ordering);


        [A,I] = sort(newListOfValues); 

        listPositions = zeros(length(listOfValues),1);

        for i = 1:length(listOfValues)
            listPositions(I(i)) = i;
        end

        %reassign to be in range where can use erfinv: need to be between -1 and 1.


        listPositions = listPositions-mean(listPositions);
        listPositions = 2.*listPositions./(length(listPositions));

        randOrdering_transformedValues = erfinv(listPositions);

        transformedValues = zeros(size(listOfValues));
        %return to the original ordering
        for i = 1:length(listOfValues)
            transformedValues(ordering(i)) = randOrdering_transformedValues(i);
        end
        x2(:,iiii) = transformedValues;
    end
    transformedValues = x2;
    
end
end
