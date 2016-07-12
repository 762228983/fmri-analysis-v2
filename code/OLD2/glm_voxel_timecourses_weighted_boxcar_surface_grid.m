function [component_weights, component_names] = glm_voxel_timecourses_weighted_boxcar_surface_grid(grid_file, para_file, weight_file, varargin)
% Runs first-level analysis based on components
% 
% Last Modified by Sam Norman-Haignere on 4/7/2015



% estimate component weights
[component_weights, component_names] = glm_voxel_timecourses_weighted_boxcar(data_matrix, TR, para_file, weight_file);
