% This script allows to easily compile the MEX functions in 
% the MI - Continuous Method.
% Once the compiler has been setup properly (see help for function MEX in
% Matlab) simply run the script. Remember that, in order for the make file
% to work properly, the toolbox folders need to be included in the Matlab
% path: see the installation instruction and function STARTUP_infolab.

% Toolbox directory:
toolboxDir = 'C:\Users\glorenz\Documents\PhD\MINT\src\MI\ContinuousMethods';

disp('Compiling files for NPC:');

% Compiling NPC_DenseNaiveMat and NPC_DenseNaiveMatwin
outDir   = fullfile(toolboxDir, 'NPC');
if 0
	compileOpts = {'-g', '-largeArrayDims'};
else
	compileOpts = {'-O', '-largeArrayDims'};
end

if ispc
    c_file = fullfile(outDir,'NPC_DenseNaiveMatwin.c');
else
    c_file = fullfile(outDir,'NPC_DenseNaiveMat.c');
end

disp(['   Compiling: ' c_file]);

mex(compileOpts{:}, '-outdir', outDir, c_file);
