function [PCs, exvar] = whitematter_PCs(exp, us, runtype, r, volname, regtype, varargin)

% Calculates principal components from white matter using a masked computed by
% Freesurfer.

global root_directory;

if nargin < 5
    volname = 'motcorr';
end

if nargin < 6
    regtype = 'bbreg';
end

[~, freesurfer_version] = read_version_info(exp);

% preprocessing directory for the functional data
preprocessing_directory = [root_directory '/' exp '/analysis/preprocess' ...
    '/usub' num2str(us) '/' runtype '_r' num2str(r)];

% create white matter mask from segmentation
freesurfer_segmentation = [root_directory '/freesurfer/us' num2str(us) '/mri/aseg.mgz'];
whitematter_mask_anatomical = [root_directory '/freesurfer/us' num2str(us) '/mri/wm_bin.mgz'];
if ~exist(whitematter_mask_anatomical,'file') || optInputs(varargin,'overwrite')
    unix_freesurfer_version(freesurfer_version, ...
        ['mri_binarize --i ' freesurfer_segmentation ...
        ' --o ' whitematter_mask_anatomical ' --match 41 2']);
end

% register to functional volume
exfunc = [preprocessing_directory '/example_func.nii.gz'];
whitematter_mask_functional = [preprocessing_directory '/wm_bin.nii.gz'];
reg_exfunc2highres = [preprocessing_directory '/' ...
    'reg_' regtype '/example_func2highres.dat'];
if ~exist(whitematter_mask_functional,'file') || optInputs(varargin,'overwrite')
    unix_freesurfer_version(freesurfer_version, ...
        ['mri_vol2vol' ...
        ' --inv ', ...
        ' --mov ' exfunc ...
        ' --targ ' whitematter_mask_anatomical ...
        ' --reg ' reg_exfunc2highres ...
        ' --o ' whitematter_mask_functional, ...
        ' --trilin']);
end

% compute the PCs
whitematter_PCs_mat_file = [preprocessing_directory '/wm_PCs.mat'];
if ~exist(whitematter_PCs_mat_file, 'file')
    
    % read in functional volume and white-matter
    func = MRIread([preprocessing_directory '/' volname '.nii.gz']);
    wm = MRIread(whitematter_mask_functional);
    
    % unwrap space dimensions and transpose
    % -> time x all voxels
    dims = size(func.vol);
    D = reshape(func.vol, [prod(dims(1:3)), dims(4)])';
    W = reshape(wm.vol, [prod(dims(1:3)),1])';
    
    % select white-matter voxels
    D = D(:,W > 0.99);
    
    % demean
    D = D - repmat(mean(D,1), size(D,1), 1);
    
    % time x PC matrix
    [PCs,S,~] = svd(D, 'econ');
    exvar = cumsum(diag(S).^2 / sum(diag(S).^2));
    save(whitematter_PCs_mat_file, 'PCs', 'exvar');
    
else
    
    load(whitematter_PCs_mat_file, 'PCs', 'exvar');
    
end

if optInputs(varargin, 'plot')
    
    % cumulative scree plot
    figure;
    plot(exvar,'k-o');
    hold on;
    xlabel('PC #'); ylabel('Cumulative Explained Variance'); ylim([0 1]);
    box off;
    export_fig([preprocessing_directory '/wm_exvar.pdf'],'-pdf','-transparent','-nocrop');
    
    % timecourses for first N PCs
    figure;
    n_PCs = 5;
    legendlabels = cell(1,n_PCs);
    cols = colormap(sprintf('jet(%d)',n_PCs))*0.8;
    for i = 1:n_PCs
        plot(PCs(:,i),'Color',cols(i,:));
        if i == 1
            hold on;
        end
        legendlabels{i} = sprintf('PC %d',i);
    end
    legend(legendlabels{:}); xlabel('TRs'); ylabel('Response');
    box off;
    export_fig([preprocessing_directory 'wm_top' num2str(n_PCs) '.pdf']...
        ,'-pdf','-transparent','-nocrop');
    
end



