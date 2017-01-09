function combine_runs_scan_parameters(exp, us, ...
    original_runtype, original_runs, varargin)

% Creates scan parameter file for combined runs
% 
% 2017-01-09: Created, Sam NH

global root_directory;

I.overwrite = false;
I.combined_runtype = [original_runtype '_combined'];
I.combined_runs = 1:length(original_runs);
I.list_all_runs_in_fname = false;
I.create_para_files = false;
I = parse_optInputs_keyvalue(varargin, I);

% list all runs in file name if the original and combined runtype are the same
if strcmp(I.combined_runtype, original_runtype)
    I.list_all_runs_in_fname = true;
end

%% Combine surface files

scan_parameter_directory = [root_directory '/' exp '/data/brain/scan-parameters'];

combined_scan_parameter_files = cell(1,length(original_runs));
for i = 1:length(original_runs)
    
    % new directory
    if I.list_all_runs_in_fname
        combined_scan_parameter_files{i} = [scan_parameter_directory '/' ...
            I.combined_runtype '_us' num2str(us) '_r' sprintf(original_runs{i}, '%d') '.txt'];
    else
         combined_scan_parameter_files{i} = [scan_parameter_directory '/' ...
            I.combined_runtype '_us' num2str(us) '_r' num2str(I.combined_runs(i)) '.txt'];
    end
    
    if ~exist(combined_scan_parameter_files{i}, 'file') || I.overwrite
        
        fprintf('Combining runs %s\n', sprintf('%d', original_runs{i})); drawnow;
        
        % combine runs
        % determine TR, TA, and nTR
        for j = 1:length(original_runs{i})
            [TR, TA, nTR, ~] = read_functional_scan_parameters(exp, us, ...
                original_runtype, original_runs{i}(j));
            
            % update parameters
            if j == 1
                P.TR = TR;
                P.TA = TA;
                P.nTR = int32(nTR);
            else
                P.nTR = P.nTR + int32(nTR);
                assert(P.TR == TR);
                assert(P.TA == TA);
            end
        end
        clear X;
        
        % write file
        fid = fopen(combined_scan_parameter_files{i}, 'w');
        fprintf(fid, 'TR %.6f\nTA %.6f\nnTR %d\nnDISDAQS NaN\n', ...
            P.TR, P.TA, P.nTR);
        fclose(fid);    
    end
end
