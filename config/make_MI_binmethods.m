% This script allows to easily compile the MEX functions in the MI - BinnedMethods
% module. Once the compiler has been setup properly (see help for function MEX in
% Matlab) simply run the script. Remember that, in order for the make file
% to work properly, the toolbox folders need to be included in the Matlab
% path: see the installation instruction and function STARTUP_infolab.

% Toolbox directory:
toolboxDir = 'C:\Users\glorenz\Documents\PhD\MINT\src\MI\BinnedMethods';
if 0
	compileOpts = {'-g', '-largeArrayDims'};
else
	compileOpts = {'-O', '-largeArrayDims'};
end

disp('Compiling files for NIT:');
% Compiling PARTITION_X -----------------------------------------------
outDir   = fullfile(toolboxDir, 'Extrapolation', 'Helpers');
filePath = fullfile(outDir, 'partition_X.c');

disp(['   Compiling: ' filePath]);

mex(compileOpts{:}, '-outdir', outDir, filePath);

% Compiling SHUFFLE_R_ACROSS_TRIALS -----------------------------------
outDir   = fullfile(toolboxDir, 'Extrapolation', 'Helpers');
filePath = fullfile(outDir, 'shuffle_X_across_trials.c');

disp(['   Compiling: ' filePath]);

mex(compileOpts{:}, '-outdir', outDir, filePath);

% Compiling SHUFFLE_X_ACROSS_CELLS ------------------------------------
outDir   = fullfile(toolboxDir, 'Methods', 'Helpers');
filePath = fullfile(outDir, 'shuffle_X_across_cells.c');

disp(['   Compiling: ' filePath]);

mex(compileOpts{:}, '-outdir', outDir, filePath);

% Compiling DIRECT_METHOD ---------------------------------------------
outDir = fullfile(toolboxDir, 'Methods', 'DirectMethod');
file1Path = fullfile(outDir, 'direct_method.c');
file2Path = fullfile(toolboxDir, 'BiasCorrection', 'DirectMethod', 'panzeri_and_treves_96.c');

disp(['   Compiling: ' file1Path]);
disp(['              ' file2Path]);

mex(compileOpts{:}, '-outdir', outDir, file1Path, file2Path);
    
disp('File compiling successful.');
