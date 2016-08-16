function mcflirt(exp,us,runtype,r,varargin)

% function mcflirt(exp,us,runtype,r,varargin)
% 
% motion corrects data using FSL's mcflirt algorithm

global root_directory

% fsl version for this experiment
[fsl_version, ~] = read_version_info(exp,varargin{:});

% functional volume to which the run is motion-corrected
[~, ~, nTR] = read_functional_scan_parameters(exp,us,runtype,r,varargin{:});
reference_frame = floor(nTR/2);

% raw data, un-motion corrected
input_file = [root_directory '/' exp '/data/brain/nifti/usub' num2str(us) '/' runtype '_r' num2str(r) '.nii.gz'];

% check nifti file exists
if ~exist(input_file, 'file')
    error('Input file below doesn''t exist: %s\n', input_file);
end

% preprocessing directory
preprocessing_directory = [root_directory '/' exp '/analysis/preprocess/usub' num2str(us) '/' runtype '_r' num2str(r)];
if ~exist(preprocessing_directory,'dir');
    mkdir(preprocessing_directory);
end

% volume used for motion-correction
example_func_file = [preprocessing_directory '/example_func.nii.gz'];

% motion corrected file
mcflirt_file = [preprocessing_directory '/motcorr.nii.gz'];

% create reference volume
if ~exist(example_func_file, 'file') || optInputs(varargin, 'overwrite');
    unix_fsl(fsl_version, ['fslroi ' input_file ' ' example_func_file ' ' num2str(reference_frame) ' 1']);
end

% motion correct
if ~exist(mcflirt_file,'file') || optInputs(varargin, 'overwrite');
    unix_fsl(fsl_version, ['mcflirt -in ' input_file ' -out ' strrep(mcflirt_file,'.nii.gz','') ' -reffile ' example_func_file ' -mats -plots -rmsrel -rmsabs']);
end

% plot rms motion, see motion_summary for more extensive plots
if optInputs(varargin, 'plot')
    fid = fopen([preprocessing_directory '/motcorr_rel.rms']);
    x = textscan(fid,'%f'); fclose(fid);
    if optInputs(varargin, 'monkey')
      relrms = x{1}(:)/2.8571;
    else
      relrms = x{1}(:);
    end
    fprintf('\n\n%s, run %d\n',runtype,r);
    fprintf('Mean: %.3f\n',mean(relrms(:)));
    fprintf('Time points > 0.5 mm: %d\n',sum(relrms(:)>0.5));
    fprintf('Time points > 1 mm: %d\n',sum(relrms(:)>1));
    figure;
    plot(relrms);
    ylim([0 1]);
    xlabel('Volumes'); ylabel('Time-Point-to-Time-Point Motion (mm)');
    title(sprintf('%s, run %d',strrep(runtype,'_',' '),r));
    box off;
    export_fig([preprocessing_directory '/relrms.pdf'],'-pdf','-nocrop','-transparent');
end