

root_directory = '/mindhive/nklab/u/svnh';

grid_roi = 'hand-stp-stg';
grid_spacing_mm = 1.5;

%%
% 
us = 123;
loc.exp = 'naturalsound-new-pipeline';
loc.runtype = 'main_v3';
loc.analysis_name = 'auditory-ICA-components';
loc.n_perms = 1000;

% us = 360;
% loc.exp = 'naturalsound-ica-localizer-v1';
% loc.runtype = 'ica_localizer';
% loc.analysis_name = 'auditory-ICA-components';
% loc.n_perms = 1000;

loc.fwhm = 3;
loc.runs = read_runs(loc.exp, us, loc.runtype);

loc_betas = zeros(6,6309);
loc_all_betas = nan(6,6309,length(loc.runs));
for i = 1:length(loc.runs)
    fname = [root_directory '/' loc.exp '/analysis/glm/' loc.analysis_name ...
        '/fsaverage_smooth-' num2str(loc.fwhm) 'mm_grid-' num2str(grid_spacing_mm) ...
        'mm_' grid_roi '/usub' num2str(us) ...
        '/r' num2str(loc.runs(i)) '_' num2str(loc.n_perms) 'perms.mat'];
    load(fname, 'beta_one_per_regressor');
    loc_betas = loc_betas + beta_one_per_regressor/length(loc.runs);
    loc_all_betas(:,:,i) = beta_one_per_regressor;
end

fname = [root_directory '/' loc.exp '/analysis/glm/' loc.analysis_name ...
    '/fsaverage_smooth-' num2str(loc.fwhm) 'mm_grid-' num2str(grid_spacing_mm) ...
    'mm_' grid_roi '/usub' num2str(us) ...
    '/allruns_' num2str(loc.n_perms) 'perms.mat'];

load(fname, 'logP_permtest');

xi = any(logP_permtest > 2);

%% test

test.exp = 'music-quilting-fMRI';
test.analysis_name = 'intact-vs-quilt-31ms';
test.n_perms = 1000;
test.runs = 1:10;
test.fwhm = 3;

test_betas = zeros(16,6309);
for i = 1:length(test.runs);
    fname = [root_directory '/' test.exp '/analysis/glm/' test.analysis_name ...
        '/fsaverage_smooth-' num2str(test.fwhm) 'mm_grid-' num2str(grid_spacing_mm) ...
        'mm_' grid_roi '/usub' num2str(us) ...
        '/r' num2str(test.runs(i)) '_' num2str(test.n_perms) 'perms.mat'];
    load(fname, 'beta_one_per_regressor');
    test_betas = test_betas + beta_one_per_regressor/length(test.runs);
end

%%

xi = all(~isnan(test_betas)) & all(~isnan(loc_betas));% & any(logP_permtest > 1.3);
betas = test_betas(:,xi) * pinv(loc_betas(:,xi)); 

