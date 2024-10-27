%%% pdf - Compute the Joint Probability Distribution of input data.
%%%
%%% Input:
%%%   data: *nDims X nTrials X time* input matrix
%%%   opts: Binning options
%%%       - opts.nb: Scalar with number of bins/Vector with the number of bins for each dimension
%%%       - opts.method: Binning methods for each dimension
%%%               Available options:
%%%                   - 'eqpop': Equipopulated binning
%%%                   - 'eqspace': Equispaced binning (default if no values is given in opts.method)
%%%                   - 'ceqspace': Centered equispaced binning
%%%                   - 'gseqspace': Gaussian equispaced binning
%%%                   - 'none': No binning, use raw data values as bins
%%%
%%% Output: probdist: Joint Probability Distribution of input parameters
