function combine_runs_runorder_file(exp, us, ...
    original_runtype, original_runs, varargin)

% Creates run order file for combined runs
%
% 2017-01-10: Created, Sam NH

global root_directory;

I.overwrite = true;
I.combined_runtype = [original_runtype '_combined'];
I.combined_runs = 1:length(original_runs);
I.create_para_files = false;
I = parse_optInputs_keyvalue(varargin, I);

%% Combine surface files

runorders_directory = [root_directory '/' exp '/data/brain/runorders'];

combined_runorder_file = [runorders_directory '/' ...
    I.combined_runtype '_us' num2str(us) '.txt'];

if ~exist(combined_runorder_file, 'file') || I.overwrite
    
    % open file for writing
    fid = fopen(combined_runorder_file, 'w');
    
    % runs and scan ids for original runtype
    [runs, ~, scanid] = read_runs(exp, us, original_runtype);
    
    for i = 1:length(original_runs)
        
        % scan id for combined run
        xi = original_runs{i}(1) == runs;
        combined_scanid = scanid(xi);
        clear xi;
        
        % write file
        fprintf(fid, '%d %d %d\n', I.combined_runs(i), -1, combined_scanid);

    end
    
    % close file
    fclose(fid);
    
end
