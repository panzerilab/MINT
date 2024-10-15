function pid_v = pidimmigaussian(X, Y)

if nargin < 1
    msg = 'not enough input arguments.';
    error('pidimin:notEnoughInput', msg);
end

infopts = [];
infopts.bias = 'naive';
infopts.method = 'gs';
infopts.bin_methodX = 'none'; %opts.bin_methodX; 
infopts.bin_methodY = 'none'; %opts.bin_methodY;
% infopts.n_binsX     = opts.n_binsX; 
% infopts.n_binsY     = opts.n_binsY;
infopts.verbose = false;
infopts.btsp= 0;

if nargin < 3
    msg = 'not enough input arguments.';
    error('extrapolation:notEnoughInput', msg);
end

I12 = cell2mat(information(vertcat(X{:}), Y, infopts, {'I'})); %mutualInformationXYZ(p);
I1  = cell2mat(information(X{1}, Y, infopts, {'I'}));
I2  = cell2mat(information(X{2}, Y, infopts, {'I'}));

red = min([I1, I2]);
u1 = I1 - red;
u2 = I2 - red;
syn = I12 - u1 - u2 - red;
pid_v = [red u1 u2 syn];
end