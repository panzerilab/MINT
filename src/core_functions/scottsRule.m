function num_bins = scottsRule(data)   
%%% Description:
%%% Calculates the number of bins using Scott's rule for histogram bin width.
%%%
%%% Inputs:
%%% - *data*: Input data vector.
%%%
%%% Outputs:
%%% - *num_bins*: Optimal number of bins calculated using Scott's rule.
%%%
%%% Source:
%%% This implementation is based on Scott's rule for histogram bin width.
%%% Reference: Scott, D. W. (1979). On optimal and data-based histograms. Biometrika, 66(3), 605-610.
%%% DOI: https://doi.org/10.1093/biomet/66.3.605

    N = numel(data);
    h = 3.49 * std(data) * N^(-1/3);
    data_range = range(data);
    if data_range == 0 || h == 0
        num_bins = 1;
    else
        num_bins = ceil(data_range / h);
    end

end
