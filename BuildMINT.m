% This MATLAB script configures and initializes the MINT project environment, including building components, 
% updating MATLAB paths, and ensuring all necessary configurations are applied at MATLAB startup. 
% Execute this script from the MINT project root directory to set up the environment.

function BuildMINT()
    scriptdir = fileparts(mfilename('fullpath'));
    % Define the path to the Ecos directory (external dependencies)
    ecosdir = fullfile(scriptdir, 'extern', 'ecos-matlab-master','bin');
    % Define the path to the config directory which holds configuration files
    configdir = fullfile(scriptdir, 'config');
    
    % Uncomment the following lines if you want to build and test MI and tools
    % MI_BuildAndTest(configdir, 0);  % Build MI components
    % tools_BuildAndTest(configdir);  % Build tools components
    
    % Build and test PID components (called when you execute this function)
    PID_BuildAndTest(configdir, ecosdir);
    compile_KSGmex_file(scriptdir, fullfile(scriptdir, 'extern', 'ContinuousMIEstimation','bin'));
    % Add the necessary paths to MATLAB startup by modifying the startupMINT.m
    addReqToStartup(scriptdir);
    
    % Print success message
    disp('_________________________________________________________________________');
    disp('Installation Successful.');
    disp('Congratulations - You are now ready to use the MINT Toolbox. Happy coding!');
end

function find_and_replace_in(templatefile, newfile, findstr, replacestr)
    % This function replaces placeholder strings in the template file with actual values 
    % and saves the modified content to a new file.
    
    replace = containers.Map;  % Create a map to store the find and replace pairs
    
    % Ensure findstr and replacestr are cell arrays for consistency
    if ~iscell(findstr)
        findstr = {findstr};
    end
    if ~iscell(replacestr)
        replacestr = {replacestr};
    end
    
    % Determine the minimum length between findstr and replacestr to avoid index issues
    minLen = min(length(findstr), length(replacestr));
    
    % Fill the replace map with find-replace pairs
    for i = 1:minLen
        replace(findstr{i}) = replacestr{i};
    end
    
    % Open the template file for reading and new file for writing
    infile = fopen(templatefile, 'r');
    outfile = fopen(newfile, 'w');
    
    % Check if the template file and new file opened successfully
    if infile == -1
        error(['Could not open file ' templatefile '.']);
    end
    if outfile == -1
        error(['Could not create file ' newfile '.']);
    end
    
    % Read through the template file line by line, replace placeholders, and write to new file
    while ~feof(infile)
        line = fgetl(infile);
        keysReplace = keys(replace);
        for i = 1:length(keysReplace)
            line = strrep(line, keysReplace{i}, replace(keysReplace{i}));
        end
        fprintf(outfile, '%s\n', line);
    end
    
    % Close the files after processing
    fclose(infile);
    fclose(outfile);
end

function add_to_startup(startupfile)
    % This function appends the content of the provided startupfile to the user's startupMINT.m script.
    % The startupMINT.m file is responsible for setting up necessary paths when MATLAB starts.
    
    % Define the path to the user's startupMINT.m file
    targetFilePath = fullfile(userpath, 'startupMINT.m');
    
    % If the file does not exist, create it
    if exist(targetFilePath, 'file') ~= 2
        fid = fopen(targetFilePath, 'w');
        fclose(fid);
        disp(['Created ' targetFilePath]);
    end 
    
    % Read the current content of the startupMINT.m file and the new content from the startupfile
    targetFileContent = fileread(targetFilePath);
    sourceFileContent = fileread(startupfile);
    
    % If the content is not already in the file, append it
    if ~contains(targetFileContent, sourceFileContent)
        stfileout = fopen(targetFilePath, 'a');
        fwrite(stfileout, sourceFileContent);
        fclose(stfileout);
        disp(['Added content of ' startupfile ' to startupMINT.m']);
    else
        disp(['Content of ' startupfile ' already exists in startupMINT.m']);
    end
end

function addReqToStartup(configdir)
    % This function adds the required paths for the MINT project to MATLAB's startupMINT.m and startup.m files.
    
    % Define the path to the startupMINT.m file
    startupMINTPath = fullfile(userpath, 'startupMINT.m');  
    
    % Read the current content of the startupMINT.m file
    stfileoutread_mint = fopen(startupMINTPath, 'r');
    stfileoutlines_mint = textscan(stfileoutread_mint, '%s', 'Delimiter', '\n');
    fclose(stfileoutread_mint);
    
    % If the end marker is not found, add the necessary paths
    if ~any(strcmp(stfileoutlines_mint{1}, '% finish adding MINT requirements'))
        stfileout_mint = fopen(startupMINTPath, 'a');
        fprintf(stfileout_mint, '%% Add all config Folders to MATLAB path.\n');
        fprintf(stfileout_mint, 'addpath(genpath("%s"));\n', configdir);
        fprintf(stfileout_mint, '%% finish adding MINT requirements\n');
        fclose(stfileout_mint);
    end

    % Define the path to the general startup.m file
    startupPath = fullfile(userpath, 'startup.m');
    
    % If the startup.m file doesn't exist, create it
    if ~exist(startupPath, 'file')
        fid = fopen(startupPath, 'w');
        fclose(fid);
    end
    
    % Read the current content of the startup.m file
    stfileoutread_startup = fopen(startupPath, 'r');
    stfileoutlines_startup = textscan(stfileoutread_startup, '%s', 'Delimiter', '\n');
    fclose(stfileoutread_startup);
    
    % If the startupMINT command isn't found, add it to the startup.m file
    if ~any(strcmp(stfileoutlines_startup{1}, 'startupMINT'))
        stfileout_startup = fopen(startupPath, 'a');
        fprintf(stfileout_startup, '%% execute requirements of the MINT toolbox\n');
        fprintf(stfileout_startup, 'startupMINT\n');
        fprintf(stfileout_startup, '%% finish adding MINT requirements\n');
        fclose(stfileout_startup);
    end
end

function PID_BuildAndTest(configdir, ecosdir)
    % This function installs the ECOS solver (via the makemex script),
    % generates configuration files for the PID module, and adds the generated startup file 
    % to MATLAB's startup sequence.

    disp('Installing Ecos');
    
    % Define the path to the makemex script for Ecos installation
    makemex_path = fullfile(ecosdir, 'makemex.m');
    
    % Run the makemex script to install Ecos if it's available
    if exist(makemex_path, 'file') == 2
        run(makemex_path);
    else
        error('The makemex script was not found in the specified directory.');
    end  
    
    disp('Generating PID config files');
    
    % Define paths for the PID startup template and output startup file
    startupfiletemplate = fullfile(configdir, 'startup_PID_TEMPLATE.m');
    startupfile = fullfile(configdir, 'startup_PID.m');
    
    % Define the source directory for the PID module
    moduleFolder = fullfile('src','PID');
    scriptdir = fileparts(mfilename('fullpath'));
    srcdir = fullfile(scriptdir, moduleFolder);
    
    % Replace the placeholder in the template with the actual PID root directory
    find_and_replace_in(startupfiletemplate, startupfile, '$PID_MATLAB_ROOT', srcdir);
    
    % Add the new startup file to MATLAB's startup sequence
    add_to_startup(startupfile);
end

function compile_KSGmex_file(scriptdir, outputdir)
    % This function compiles the Mxnyn.C file using the mex command and places the output in the specified folder
    
    % Define the path to the Mxnyn.C file
    mexfile = fullfile(scriptdir, 'extern', 'ContinuousMIEstimation', 'MIxnynmint.C');
    
    % Check if the Mxnyn.C file exists
    if exist(mexfile, 'file') == 2
        % Ensure the output directory exists
        if ~exist(outputdir, 'dir')
            mkdir(outputdir);  % Create the output directory if it doesn't exist
        end
        
        % If the file exists, compile the file and specify the output folder
        disp('Compiling Mxnynmint.C using mex...');
        mex('-outdir', outputdir, mexfile);  % Compile the C file and specify output directory
    else
        % If the file does not exist, display an error message
        error('The Mxnynmint.C file was not found in the script directory.');
    end
end
