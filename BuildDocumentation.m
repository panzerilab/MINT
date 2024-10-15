function BuildDocumentation()
    main();
end

function main()
    
    docsdirname = 'Documentation';
    MINTrootdir = pwd;
    srcdir = fullfile(MINTrootdir, 'src');
    docsdir = fullfile(MINTrootdir, docsdirname);
   
    
    if exist(docsdir, 'dir')
        rmdir(docsdir, 's');
    end
         
    dirs = {srcdir};
    markdown_files = {};
    ignorelist = {};
    
    while ~isempty(dirs)
        dirname = dirs{end};
        dirs(end) = [];
        docsignore = fullfile(dirname, '.docsignore');
    
        if exist(docsignore, 'file')
            ignorelist = textscan(docsignore, '%s', 'delimiter', '\n');
            ignorelist = cellfun(@(x) strrep(x, '*', '\\\*'), ignorelist, 'UniformOutput', false);
            ignoreregexps = strjoin(string(ignorelist), '|');
        end
    
        files = dir(dirname);
        files = files(~startsWith({files.name}, '.'));
    
        for entry = files'
            if entry.isdir && ~strcmp(entry.name, docsdirname)
                dirs{end+1} = fullfile(dirname, entry.name);
            elseif endsWith(entry.name, '.md') && ...
                    ((isempty(ignorelist) || isempty(regexp(entry.name, ignoreregexps, 'once'))) || ...
                    (strcmp(entry.name, 'Readme.m') && endsWith(entry.name, '.m')))
                relativepath = fullfile(relpath(entry.folder, srcdir), entry.name);
                absolutedocpath = fullfile(docsdir, relativepath);
                if ~exist(fileparts(absolutedocpath), 'dir')
                    mkdir(fileparts(absolutedocpath));
                end
                copyfile(fullfile(entry.folder, entry.name), absolutedocpath);
                markdown_files{end+1} = absolutedocpath;
            else
                for ex = {'.m'}
                    if endsWith(entry.name, ex) && ...
                            (isempty(ignorelist) || isempty(regexp(entry.name, ignoreregexps, 'once')))
                        rgx = ['\' ex '$'];
                        relativepath = fullfile(relpath(entry.folder, srcdir), regexprep(entry.name, '\.m$', '.md'));
                        absolutedocpath = fullfile(docsdir, relativepath);
                        disp(['Processing file: ' entry.name]);
                        createfile = createDocsFromSrc(fullfile(entry.folder, entry.name), absolutedocpath, '%%%');
                        if createfile
                            markdown_files{end+1} = absolutedocpath;
                        end
                    end
                end
            end
        end
        ignorelist = {};
    end
    
    rootfiles = dir(MINTrootdir);
    for entry = rootfiles'
        if entry.isdir
            continue;
        end
        if endsWith(entry.name, '.md')
            if strcmp(entry.name, 'README.md')
                copyfile(fullfile(MINTrootdir, entry.name), fullfile(docsdir, 'home.md'));
            else
                copyfile(fullfile(MINTrootdir, entry.name), fullfile(docsdir, entry.name));
            end
        end
    end 
    buildSidebar(docsdir);
    
end

function buildSidebar(rootDocsPath)
    sidebarFile = fopen(fullfile(rootDocsPath, '_Sidebar.md'), 'w');
    rgx = '\.md$';
    [allFiles, ~, ~] = dirwalk(rootDocsPath);

    for i = 1:numel(allFiles)
        file = allFiles{i};
        [~, name, ext] = fileparts(file);
        level = count(file, filesep)-7;
        lastTwoChars = file(end-1:end);
    
        if level == 0 && strcmp(lastTwoChars, 'md') && ~strcmp(name, '_Sidebar')
            fprintf(sidebarFile, '%s* [%s](.%s)\n', repmat('  ', 1, level), name, ['/' name]);
        end
    end

    for i = 1:numel(allFiles)
        file = allFiles{i};
        [~, name, ext] = fileparts(file);
        level = count(file, filesep)-7;
        lastTwoChars = file(end-1:end);
        [~, nameFolder, ~] = fileparts(fileparts(file));
    
        if level > 0
            if strcmp(lastTwoChars, '/.')
                fprintf(sidebarFile, '%s* %s\n', repmat('  ', 1, level-1), nameFolder);
            elseif strcmp(lastTwoChars, 'md')
                parts = strsplit(file, '/');
                extractedPath = strjoin(parts(end - level:end-1), '/');
                fprintf(sidebarFile, '%s* [%s](%s/%s)\n', repmat('  ', 1, level), name, extractedPath, name);
            end
        end
    end
    fclose(sidebarFile);
end

function createfile = createDocsFromSrc(srcFilePath, docsFilePath, identifier)
    warning('off', 'all');    
    srcFile = fopen(srcFilePath, 'r');
    docsFile = [];
    mkdir(fileparts(docsFilePath));       
    docsFile = fopen(docsFilePath, 'w', 'n', 'UTF-8');
    createfile = true;
    
    if srcFile == -1
        error('Error opening source file.');
    end
    
    docsFile = fopen(docsFilePath, 'w');
    if docsFile == -1
        fclose(srcFile);
        error('Error opening target file.');
    end
    while ~feof(srcFile)
        line = fgetl(srcFile);
        if startsWith(line, identifier)
            fprintf(docsFile, '%s\n', line);
        end
    end
    
    fclose(srcFile);
    fclose(docsFile);
    warning('on', 'all');
end

function rpath = relpath(path, basepath) 
    commonPrefix = commonprefix(path, basepath);
    rpath = strtrim(path(length(commonPrefix)+1:end));
end

function common = commonprefix(str1, str2)
    len = min(length(str1), length(str2));
    common = '';
    for i = 1:len
        if str1(i) == str2(i)
            common = [common, str1(i)];
        else
            break;
        end
    end
end

function [allFiles, allDirs, allPaths] = dirwalk(directory)
    allFiles = {};
    allDirs = {};
    allPaths = {};
    entries = dir(directory);
    for entry = entries'
        if entry.isdir && ~strcmp(entry.name, '.') && ~strcmp(entry.name, '..')
            subdir = fullfile(directory, entry.name);
            [subFiles, subDirs, subPaths] = dirwalk(subdir);
            allFiles = [allFiles; subFiles];
            allDirs = [allDirs; subDirs];
            allPaths = [allPaths; subPaths];
        else            
            allFiles = [allFiles; fullfile(directory, entry.name)];
            allDirs = [allDirs; entry.isdir];
            allPaths = [allPaths; fullfile(directory, entry.name)];
        end
    end
end
