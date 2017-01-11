function combined_grid_files = combine_runs_surf_grid(exp, us, original_runtype, ...
    original_runs, fwhm, grid_roi, grid_spacing_mm, varargin)

% Combines surface grid files across sets of runs and for each group creates a
% new run of a new combined runtype.
% 
% 2017-01-08: Created, Sam NH

global root_directory;

% optional input arguments
I.overwrite = false;
I.detrend = 0;
I.combined_runtype = [original_runtype '_combined'];
I.combined_runs = 1:length(original_runs);
I.list_all_runs_in_fname = false;
I.create_para_files = false;
I.keyboard = false;
I = parse_optInputs_keyvalue(varargin, I);

% list all runs in file name if the original and combined runtype are the same
if strcmp(I.combined_runtype, original_runtype)
    I.list_all_runs_in_fname = true;
end

if I.keyboard
    keyboard;
end

%% Combine surface files

preprocessing_directory = [root_directory '/' exp ...
    '/analysis/preprocess/usub' num2str(us)];

combined_grid_files = cell(1,length(original_runs));
for i = 1:length(original_runs)
    
    % new directory
    if I.list_all_runs_in_fname
        combined_preprocessing_directory = [preprocessing_directory '/' ...
            I.combined_runtype '_r' sprintf(original_runs{i}, '%d') '/myfsaverage'];
    else
         combined_preprocessing_directory = [preprocessing_directory '/' ...
            I.combined_runtype '_r' num2str(I.combined_runs(i)) '/myfsaverage'];
    end
    
    % create if not already present
    if ~exist(combined_preprocessing_directory, 'dir')
        mkdir(combined_preprocessing_directory);
    end
    
    % new file
    combined_grid_files{i} = [combined_preprocessing_directory '/smooth-' ...
        num2str(fwhm) 'mm_grid-' num2str(grid_spacing_mm) 'mm_' grid_roi '.mat'];
    
    if ~exist(combined_grid_files{i}, 'file') || I.overwrite
        
        fprintf('Combining runs %s\n', sprintf('%d', original_runs{i})); drawnow;
        
        % concatenate runs
        D = [];
        run_index = [];
        for j = 1:length(original_runs{i})
            original_grid_file = [preprocessing_directory '/' ...
                original_runtype '_r' num2str(original_runs{i}(j)) '/myfsaverage'...
                '/smooth-' num2str(fwhm) 'mm_grid-' num2str(grid_spacing_mm) 'mm' ...
                '_' grid_roi '.mat'];
            
            % load grid file and convert to a matrix
            load(original_grid_file, 'G');
            X = grid2matrix(G);
            D = [D; X]; %#ok<AGROW>
            run_index = [run_index; j * ones(size(X,1),1)]; %#ok<AGROW>
        end
        
        % optionally detrend
        if I.detrend >= 0
            
            D_detrended = nan(size(D));
            for j = 1:length(original_runs{i})
                xi = run_index == j;
                                
                % detrending regression matrix
                X = (1:sum(xi))'; % time vector
                X = X-mean(X); % demean time vector
                X = X * ones(1, I.detrend+1); % replicate to matrix
                for z = 1:I.detrend+1
                    X(:,z) = X(:,z).^(z-1);
                end
                
                % detrend
                D_detrended(xi,:) = D(xi,:) ...
                    - X * (pinv(X) * D(xi,:));
            end
            
            % add the global mean back
            D = bsxfun(@plus, D_detrended, nanmean(D,1));
            
        end
        
        % convert back to grid
        G = matrix2grid_v2(D, G);
        
        % save grid file
        save(combined_grid_files{i}, 'G');
    
    end
end

