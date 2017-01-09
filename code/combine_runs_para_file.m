function combined_para_files = combine_runs_para_file(exp, us, ...
    original_para_prefix, original_runs, varargin)

% Combines para files across sets of runs and for each group creates a
% new run of a new combined runtype.
% 
% 2017-01-08: Created, Sam NH

global root_directory;

% optional input arguments
I.overwrite = false;
I.combined_para_prefix = [original_para_prefix '_combined'];
I.combined_runs = 1:length(original_runs);
I.list_all_runs_in_fname = false;
I = parse_optInputs_keyvalue(varargin, I);

% if prefix for the original and combined is the same, then list all runs
if strcmp(I.combined_para_prefix, original_para_prefix)
    I.list_all_runs_in_fname = true;
end

combined_para_files = cell(1,length(original_runs));
for i = 1:length(original_runs)
    
    % new combined para file
    if I.list_all_runs_in_fname
        combined_para_files{i} = [root_directory '/' exp '/data/para/usub' num2str(us) ...
            '/' I.combined_para_prefix  '_r' sprintf(original_runs{i}) '.par'];
    else
        combined_para_files{i} = [root_directory '/' exp '/data/para/usub' num2str(us) ...
            '/' I.combined_para_prefix  '_r' num2str(I.combined_runs(i)) '.par'];
    end
    
    % create if it doesn't exist
    if ~exist(combined_para_files{i},'file') || I.overwrite
        
        fprintf('Combining runs %s\n', sprintf('%d', original_runs{i}));
        drawnow;
        
        % open file to write to
        fid = fopen(combined_para_files{i}, 'w');
        run_onset = 0;
        for j = 1:length(original_runs{i})
            
            % read the original para file
            original_para_file = [root_directory '/' exp '/data/para/usub' num2str(us) ...
                '/' original_para_prefix  '_r' num2str(original_runs{i}(j)) '.par'];
            P = read_para(original_para_file, varargin);
            
            % write each line of the new para file
            for k = 1:length(P.onsets)
                fprintf(fid, '%8.2f%5d%8.2f%5d %s\n', ...
                    P.onsets(k) + run_onset, P.condition_indices(k), P.durs(k), 1, P.conds{k});
            end
            
            % add the duration of this run to the onset time of the nexxt
            [TR, ~, nTR, ~] = read_functional_scan_parameters(...
                exp,us,original_para_prefix,original_runs{i}(j));
            run_onset = run_onset + TR * nTR;

        end
        
        % close file
        fclose(fid);
    end    
end

