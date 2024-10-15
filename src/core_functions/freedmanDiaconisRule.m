function num_bins = freedmanDiaconisRule(data)
%%% Description:
%%% Calculates the optimal number of bins using the Freedman-Diaconis rule for histogram bin width.
%%%
%%% Inputs:
%%% - *data*: Input data vector.
%%%
%%% Outputs:
%%% - *num_bins*: Optimal number of bins calculated using the Freedman-Diaconis rule.
%%%
%%% Source:
%%% This implementation is based on the Freedman-Diaconis rule for histogram bin width.
%%% Reference: Freedman, D. and Diaconis, P. (1981). On the histogram as a density estimator: L2 theory.
%%% DOI: https://doi.org/10.1007/BF01025868

    N = numel(data);
    IQR = iqr(data); % Interquartile range
    h = 2 * IQR * N^(-1/3); % Bin width based on Freedman-Diaconis rule
    data_range = range(data);
    
    if data_range == 0 || h == 0
        num_bins = 1;
    else
        num_bins = ceil(data_range / h);
    end

end