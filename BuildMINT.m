% This MATLAB script configures and initializes the MINT project environment, including building components, 
% updating MATLAB paths, and ensuring all necessary configurations are applied at MATLAB startup. Execute from the 
% MINT project root directory.

function BuildMINT()
    ecosdir = fullfile(pwd, 'extern', 'ecos-matlab-master','bin');
    configdir = fullfile(pwd, 'config');
    % MI_BuildAndTest(configdir, 0);
    % tools_BuildAndTest(configdir)
    PID_BuildAndTest(configdir, ecosdir)
    addReqToStartup(configdir);
    disp('_________________________________________________________________________');
    disp('Installation Successfull.');
    disp('Congratulations - You are now ready to use the MINT Toolbox. Happy coding!');
end


function find_and_replace_in(templatefile, newfile, findstr, replacestr)
    replace = containers.Map;
    if ~iscell(findstr)
        findstr = {findstr};
    end
    if ~iscell(replacestr)
        replacestr = {replacestr};
    end
    minLen = min(length(findstr), length(replacestr));
    for i = 1:minLen
        replace(findstr{i}) = replacestr{i};
    end
    infile = fopen(templatefile, 'r');
    outfile = fopen(newfile, 'w');
    if infile == -1
        error(['Could not open file ' templatefile '.']);
    end
    if outfile == -1
        error(['Could not create file ' newfile '.']);
    end
    while ~feof(infile)
        line = fgetl(infile);
        keysReplace = keys(replace);
        for i = 1:length(keysReplace)
            line = strrep(line, keysReplace{i}, replace(keysReplace{i}));
        end
        fprintf(outfile, '%s\n', line);
    end
    fclose(infile);
    fclose(outfile);
end

function add_to_startup(startupfile)
    targetFilePath = fullfile(userpath, 'startupMINT.m');
    if exist(targetFilePath, 'file') ~= 2
        fid = fopen(targetFilePath, 'w');
        fclose(fid);
        disp(['Created ' targetFilePath]);
    end 
    targetFileContent = fileread(targetFilePath);
    sourceFileContent = fileread(startupfile);
    if ~contains(targetFileContent, sourceFileContent)
        stfileout = fopen(targetFilePath, 'a');
        fwrite(stfileout, sourceFileContent);
        fclose(stfileout);
        disp(['Added content of ' startupfile ' to startupMINT.m']);
    else
        disp(['Content of ' startupfile ' already exists in startupMINT.m']);
    end
end

function MI_BuildAndTest(configdir, debugflag)
    disp('Generating MI config files');
    startupfiletemplate = fullfile(configdir, 'startup_MI_TEMPLATE.m');
    startupfile = fullfile(configdir, 'startup_MI.m');
    moduleFolder = fullfile('src','MI');
    srcdir = fullfile(pwd, moduleFolder);
    find_and_replace_in(startupfiletemplate, startupfile, '$MI_ROOT', srcdir);
    add_to_startup(startupfile);
end

function tools_BuildAndTest(configdir)
    disp('Generating Tools config files');
    startupfiletemplate = fullfile(configdir, 'startup_tools_TEMPLATE.m');
    startupfile = fullfile(configdir, 'startup_tools.m');
    moduleFolder = fullfile('src','tools');
    srcdir = fullfile(pwd, moduleFolder);
    find_and_replace_in(startupfiletemplate, startupfile, '$TOOLS_SRC_ROOT', srcdir);
    add_to_startup(startupfile);
end

function addReqToStartup(configdir)
    startupMINTPath = fullfile(userpath, 'startupMINT.m');  
    stfileoutread_mint = fopen(startupMINTPath, 'r');
    stfileoutlines_mint = textscan(stfileoutread_mint, '%s', 'Delimiter', '\n');
    fclose(stfileoutread_mint);
    if ~any(strcmp(stfileoutlines_mint{1}, ['% finish adding MINT requirements']))
        stfileout_mint = fopen(startupMINTPath, 'a');
        fprintf(stfileout_mint, '%% Add all config Folders to MATLAB path.\n');
        fprintf(stfileout_mint, 'addpath(genpath("%s"));\n', configdir);
        fprintf(stfileout_mint, '%% finish adding MINT requirements\n');
        fclose(stfileout_mint);
    end
    startupPath = fullfile(userpath, 'startup.m');
    if ~exist(startupPath, 'file')
            fid = fopen(startupPath, 'w');
            fclose(fid);
    end
    stfileoutread_startup = fopen(startupPath, 'r');
    stfileoutlines_startup = textscan(stfileoutread_startup, '%s', 'Delimiter', '\n');
    fclose(stfileoutread_startup);
    if ~any(strcmp(stfileoutlines_startup{1}, 'startupMINT'))
        stfileout_startup = fopen(startupPath, 'a');
        fprintf(stfileout_startup, '%% execute requirements of the MINT toolbox\n');
        fprintf(stfileout_startup, 'startupMINT\n');
        fprintf(stfileout_startup, '%% finish adding MINT requirements\n');
        fclose(stfileout_startup);
    end
end

function PID_BuildAndTest(configdir, ecosdir)
    disp('Installing Ecos');
    makemex_path = fullfile(ecosdir, 'makemex.m');
    if exist(makemex_path, 'file') == 2
        run(makemex_path);
    else
        error('The makemex script was not found in the specified directory.');
    end  
    disp('Generating Tools config files');
    startupfiletemplate = fullfile(configdir, 'startup_PID_TEMPLATE.m');
    startupfile = fullfile(configdir, 'startup_PID.m');
    moduleFolder = fullfile('src','PID');
    srcdir = fullfile(pwd, moduleFolder);
    find_and_replace_in(startupfiletemplate, startupfile, '$PID_MATLAB_ROOT', srcdir);
    add_to_startup(startupfile);
end


